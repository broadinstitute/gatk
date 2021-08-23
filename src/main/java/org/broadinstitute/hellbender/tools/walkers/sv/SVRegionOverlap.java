package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
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
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.IntervalOverlapCalculator;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.utils.*;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Annotate records from a structural variant (SV) VCF with overlap metrics against one or more interval sets.
 *
 * Note that -L/-XL arguments maintain normal behavior of filtering variants by location and will not be used for annotation.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         SV VCF
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
 *       --sequence-dictionary ref.dict \
 *       --region-file my_regions.bed \
 *       --region-name MY_REGIONS \
 *       -V variants.vcf.gz \
 *       -O annotated.vcf.gz
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */

@CommandLineProgramProperties(
        summary = "Annotates structural variants with overlap metrics against one or more interval sets",
        oneLineSummary = "Annotates structural variants with overlap metrics against one or more interval sets",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@BetaFeature
@DocumentedFeature
public final class SVRegionOverlap extends VariantWalker {
    public static final String REGIONS_FILE_LONG_NAME = "region-file";
    public static final String REGIONS_NAME_LONG_NAME = "region-name";
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
            doc = "Region interval files, may be specified multiple times",
            fullName = REGIONS_FILE_LONG_NAME
    )
    private List<GATKPath> regionPaths;

    @Argument(
            doc = "Region names. All values must be unique after converting to upper-case and must correspond with the " +
                    "input order of --" + REGIONS_FILE_LONG_NAME,
            fullName = REGIONS_NAME_LONG_NAME
    )
    private List<String> regionNames;

    @Argument(
            doc = "Region interval set rule",
            fullName = REGIONS_SET_RULE_LONG_NAME,
            optional=true
    )
    private IntervalSetRule intervalSetRule = IntervalSetRule.UNION;

    @Argument(
            doc = "Region interval merging rule",
            fullName = REGIONS_MERGING_RULE_LONG_NAME,
            optional=true
    )
    private IntervalMergingRule intervalMergingRule = IntervalMergingRule.OVERLAPPING_ONLY;

    @Argument(
            doc = "Region padding (bp)",
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
    private final Map<String, IntervalOverlapCalculator> intervalTreeMap = new HashMap<>();
    private VariantContextWriter writer;

    @Override
    public void onTraversalStart() {
        // Dictionary defined by input vcf
        dictionary = getMasterSequenceDictionary();

        Utils.validateArg(!(suppressOverlapFraction && suppressEndpointCounts), "Cannot use both --" +
                SUPPRESS_ENDPOINT_COUNTS_LONG_NAME + " and --" + SUPPRESS_OVERLAP_FRACTION_LONG_NAME);

        // Load interval sets
        Utils.validateArg(regionPaths.size() == regionNames.size(),
                "Number of --" + REGIONS_NAME_LONG_NAME + " and --" + REGIONS_FILE_LONG_NAME + " arguments must be equal");
        Utils.validateArg(dictionary != null, "Sequence dictionary not found in variants header");
        formattedRegionNames = regionNames.stream().map(String::toUpperCase).collect(Collectors.toList());
        Utils.validateArg(new HashSet<>(formattedRegionNames).size() == formattedRegionNames.size(), "Found duplicate region names (not case-sensitive)");
        for (int i = 0; i < regionPaths.size(); i++) {
            final IntervalOverlapCalculator calc = IntervalOverlapCalculator.create(
                    regionPaths.get(i),
                    dictionary,
                    intervalSetRule,
                    intervalMergingRule,
                    regionPadding
            );
            intervalTreeMap.put(formattedRegionNames.get(i), calc);
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

    public static String getFieldName(final String fieldBaseName, final String regionSetName) {
        return fieldBaseName + regionSetName;
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext,
                      final ReferenceContext referenceContext, final FeatureContext featureContext) {
        final SVCallRecord record = SVCallRecordUtils.create(variant);
        final VariantContextBuilder builder = new VariantContextBuilder(variant);
        for (final Map.Entry<String, IntervalOverlapCalculator> entry : intervalTreeMap.entrySet()) {
            if (!suppressEndpointCounts) {
                builder.attribute(getFieldName(GATKSVVCFConstants.NUM_END_OVERLAPS_INFO_BASE, entry.getKey()), entry.getValue().getEndpointOverlapCount(record));
            }
            if (!suppressOverlapFraction) {
                builder.attribute(getFieldName(GATKSVVCFConstants.OVERLAP_FRACTION_INFO_BASE, entry.getKey()), entry.getValue().getOverlapFraction(record));
            }
        }
        writer.add(builder.make());
    }
}
