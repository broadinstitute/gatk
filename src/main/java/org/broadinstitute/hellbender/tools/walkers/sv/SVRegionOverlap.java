package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.stratify.SVStratificationEngine;
import org.broadinstitute.hellbender.tools.sv.stratify.SVStratificationEngineArgumentsCollection;
import org.broadinstitute.hellbender.utils.*;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Annotate records from a structural variant (SV) VCF with overlap metrics against one or more interval tracks. For each
 * intervals track with name "MYTRACK", the following INFO fields are added:
 * <ul>
 *     <li>NUM_END_OVERLAPS_MYTRACK: number of ends of the variant that overlap the track [integer].</li>
 *     <li>OVERLAP_FRAC_MYTRACK: fraction of the variant that spans the track [float], or NaN for interchromosomal
 *     and SVs with 0 length on the reference.</li>
 * </ul>
 *
 * All track names are converted to upper-case. Note that -L/-XL and associated arguments maintain usual behavior
 * of filtering variants by location and are not used for annotation.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         SV VCF
 *     </li>
 *     <li>
 *         One or more interval tracks (interval file and name)
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         SV VCF with overlap annotations
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk SVRegionOverlap \
 *       --track-intervals intervals1.bed \
 *       --track-name TRACK1 \
 *       --track-intervals intervals2.bed \
 *       --track-name TRACK2 \
 *       -V variants.vcf.gz \
 *       -O annotated.vcf.gz
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */

@CommandLineProgramProperties(
        summary = "Annotates structural variants with overlap metrics against one or more interval tracks",
        oneLineSummary = "Annotates structural variants with overlap metrics against one or more interval tracks",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@BetaFeature
@DocumentedFeature
public final class SVRegionOverlap extends VariantWalker {

    public static final String REGIONS_FILE_LONG_NAME = SVStratificationEngineArgumentsCollection.TRACK_INTERVAL_FILE_LONG_NAME;
    public static final String REGIONS_NAME_LONG_NAME = SVStratificationEngineArgumentsCollection.TRACK_NAME_LONG_NAME;
    public static final String REGIONS_SET_RULE_LONG_NAME = "region-set-rule";
    public static final String REGIONS_MERGING_RULE_LONG_NAME = "region-merging-rule";
    public static final String REGION_PADDING_LONG_NAME = "region-padding";
    public static final String SUPPRESS_ENDPOINT_COUNTS_LONG_NAME = "suppress-endpoint-counts";
    public static final String SUPPRESS_OVERLAP_FRACTION_LONG_NAME = "suppress-overlap-fraction";

    @Argument(
            doc = "Output VCF",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private GATKPath outputFile;

    @Argument(
            doc = "Track interval files, may be specified multiple times",
            fullName = REGIONS_FILE_LONG_NAME
    )
    private List<GATKPath> regionPaths;

    @Argument(
            doc = "Track names, may be specified multiple times. All values must be unique after converting to " +
                    "upper-case and must correspond with the input order of --" + REGIONS_FILE_LONG_NAME,
            fullName = REGIONS_NAME_LONG_NAME
    )
    private List<String> regionNames;

    @Argument(
            doc = "Track interval set rule (applies to all tracks)",
            fullName = REGIONS_SET_RULE_LONG_NAME,
            optional=true
    )
    private IntervalSetRule intervalSetRule = IntervalSetRule.UNION;

    @Argument(
            doc = "Track interval merging rule (applies to all tracks)",
            fullName = REGIONS_MERGING_RULE_LONG_NAME,
            optional=true
    )
    private IntervalMergingRule intervalMergingRule = IntervalMergingRule.OVERLAPPING_ONLY;

    @Argument(
            doc = "Track padding (in bp, applies to all tracks)",
            fullName = REGION_PADDING_LONG_NAME,
            optional=true
    )
    private int regionPadding = 0;

    @Argument(
            doc = "Suppress endpoint counts annotation",
            fullName = SUPPRESS_ENDPOINT_COUNTS_LONG_NAME,
            optional = true
    )
    private boolean suppressEndpointCounts = false;

    @Argument(
            doc = "Suppress overlap fraction annotation",
            fullName = SUPPRESS_OVERLAP_FRACTION_LONG_NAME,
            optional = true
    )
    private boolean suppressOverlapFraction = false;

    private SAMSequenceDictionary dictionary;
    private List<String> formattedRegionNames;
    private SVStratificationEngine engine;
    private VariantContextWriter writer;

    @Override
    public void onTraversalStart() {
        // Dictionary defined by input vcf
        dictionary = getSequenceDictionaryForDrivingVariants();

        Utils.validateArg(!(suppressOverlapFraction && suppressEndpointCounts), "Cannot use both --" +
                SUPPRESS_ENDPOINT_COUNTS_LONG_NAME + " and --" + SUPPRESS_OVERLAP_FRACTION_LONG_NAME);

        // Load interval sets
        Utils.validateArg(regionPaths.size() == regionNames.size(),
                "Number of --" + REGIONS_NAME_LONG_NAME + " and --" + REGIONS_FILE_LONG_NAME + " arguments must be equal");
        Utils.validateArg(dictionary != null, "Sequence dictionary not found in variants header");
        formattedRegionNames = regionNames.stream().map(String::toUpperCase).collect(Collectors.toList());
        Utils.validateArg(new HashSet<>(formattedRegionNames).size() == formattedRegionNames.size(), "Found duplicate region names (not case-sensitive)");

        final Iterator<String> nameIterator = formattedRegionNames.iterator();
        final Iterator<GATKPath> pathIterator = regionPaths.iterator();
        final GenomeLocParser genomeLocParser = new GenomeLocParser(dictionary);
        engine = new SVStratificationEngine(dictionary);
        while (nameIterator.hasNext() && pathIterator.hasNext()) {
            final String name = nameIterator.next();
            final GATKPath intervalsPath = pathIterator.next();
            final GenomeLocSortedSet genomeLocs = IntervalUtils.loadIntervals(Collections.singletonList(intervalsPath.toString()), intervalSetRule, intervalMergingRule, regionPadding, genomeLocParser);
            final List<Locatable> intervals = Collections.unmodifiableList(genomeLocs.toList());
            engine.addTrack(name, intervals);
            engine.addStratification(name, null, null, null, Collections.singleton(name));
        }

        // Initialize output
        writer = createVCFWriter(outputFile);
        writeVCFHeader();
    }

    @Override
    public Object onTraversalSuccess() {
        writer.close();
        return null;
    }

    private void writeVCFHeader() {
        final VCFHeader header = new VCFHeader(getHeaderForVariants());
        for (final String name : formattedRegionNames) {
            if (!suppressEndpointCounts) {
                header.addMetaDataLine(new VCFInfoHeaderLine(getFieldName(GATKSVVCFConstants.NUM_END_OVERLAPS_INFO_BASE, name), 1, VCFHeaderLineType.Integer, "Number of variant endpoints overlapping region " + name));
            }
            if (!suppressOverlapFraction) {
                header.addMetaDataLine(new VCFInfoHeaderLine(getFieldName(GATKSVVCFConstants.OVERLAP_FRACTION_INFO_BASE, name), 1, VCFHeaderLineType.Float, "Fraction overlap of region " + name));
            }
        }
        writer.writeHeader(header);
    }

    private static String getFieldName(final String fieldBaseName, final String regionSetName) {
        return fieldBaseName + regionSetName;
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext,
                      final ReferenceContext referenceContext, final FeatureContext featureContext) {
        final SVCallRecord record = SVCallRecordUtils.create(variant, dictionary);
        final VariantContextBuilder builder = new VariantContextBuilder(variant);
        for (final SVStratificationEngine.Stratum stratum : engine.getStrata()) {
            if (!suppressEndpointCounts) {
                builder.attribute(getFieldName(GATKSVVCFConstants.NUM_END_OVERLAPS_INFO_BASE, stratum.getName()), stratum.countBreakpointOverlaps(record));
            }
            if (!suppressOverlapFraction) {
                builder.attribute(getFieldName(GATKSVVCFConstants.OVERLAP_FRACTION_INFO_BASE, stratum.getName()), stratum.trackOverlapFraction(record));
            }
        }
        writer.add(builder.make());
    }
}
