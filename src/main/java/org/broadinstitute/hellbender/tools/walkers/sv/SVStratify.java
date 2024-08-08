package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.SVStatificationEngine;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.nio.file.Path;
import java.util.*;

/**
 * <p>Clusters structural variants based on coordinates, event type, and supporting algorithms. Primary use cases include:</p>
 * <ul>
 *     <li>
 *         Clustering SVs produced by multiple callers, based on interval overlap, breakpoint proximity, and sample overlap.
 *     </li>
 *     <li>
 *         Merging multiple SV VCFs with disjoint sets of samples and/or variants.
 *     </li>
 *     <li>
 *         Defragmentation of copy number variants produced with depth-based callers.
 *     </li>
 * </ul>
 *
 * <p>Clustering tasks can be accomplished using one of two algorithms. The SINGLE_LINKAGE algorithm produces clusters
 * for which all members cluster with <i>at least one</i> other member. The MAX_CLIQUE algorithm, however,
 * requires that all members cluster with <i>every other</i> member. The latter is in general non-polynomial in time and
 * space but implemented to minimize computations by traversing variants ordered by start position and efficiently
 * finalizing "active" clusters that are determined to be complete.</p>
 *
 * <p>The tool determines whether two given variants should cluster based following criteria:</p>
 * <ul>
 *     <li>
 *         Matching SV type. DEL and DUP are considered matching SV types if --enable-cnv is used and merged into
 *         a multi-allelic CNV type.
 *     </li>
 *     <li>
 *         Matching breakend strands (BND and INV only)
 *     </li>
 *     <li>
 *         Interval reciprocal overlap (inapplicable for BNDs).
 *     </li>
 *     <li>
 *         Distance between corresponding event breakends (breakend window).
 *     </li>
 *     <li>
 *         Sample reciprocal overlap, based on carrier status determined by available genotypes (GT fields). If
 *         no GT fields are called for a given variant, the tool attempts to find carriers based on copy number
 *         (CN field) and sample ploidy (as determined by the ECN FORMAT field).
 *     </li>
 * </ul>
 *
 * <p>For CNV defragmentation (DEFRAGMENT_CNV algorithm), the tool uses single-linkage clustering based
 * on the following linkage criteria:</p>
 * <ul>
 *     <li>
 *         Must be a DEL/DUP/CNV and only be supported by a depth algorithm.
 *     </li>
 *     <li>
 *         Matching SV type
 *     </li>
 *     <li>
 *         Overlapping after padding both sides of each variant by the fraction of event length specified by
 *         --defrag-padding-fraction.
 *     </li>
 *     <li>
 *         Sample overlap fraction (described above) specified by --defrag-sample-overlap.
 *     </li>
 * </ul>
 *
 * <p>Interval overlap, breakend window, and sample overlap parameters are defined for three combinations of event types
 * using the ALGORITHMS field, which describes the type of evidence that was used to call the variant:</p>
 * <ul>
 *     <li>
 *         Depth-only - both variants have solely "depth" ALGORITHMS
 *     </li>
 *     <li>
 *         PESR (paired-end/split-read) - both variants have at least one non-depth entry in ALGORITHMS
 *     </li>
 *     <li>
 *         Mixed - one variant is depth-only and the other is PESR
 *     </li>
 * </ul>
 *
 * <p>Users must supply one or more VCFs containing SVs with the following info fields:</p>
 *
 * <ul>
 *     <li>
 *         SVTYPE - event type (DEL, DUP, CNV, INS, INV, BND)
 *     </li>
 *     <li>
 *         SVLEN - variant length for INS only, if known
 *     </li>
 *     <li>
 *         STRANDS - breakend strands ("++", "+-", "-+", or "--") (BND and INV only)
 *     </li>
 *     <li>
 *         CHROM2 / END2 - mate coordinates (BND only)
 *     </li>
 *     <li>
 *         ALGORITHMS - list of supporting tools/algorithms. These names may be arbitrary, except for depth-based
 *         callers, which should be listed as "depth".
 *     </li>
 * </ul>
 *
 * <p>In addition, the following FORMAT fields must be defined:</p>
 *
 * <ul>
 *     <li>
 *         GT - genotype alleles
 *     </li>
 *     <li>
 *         ECN - expected copy number (i.e. ploidy)
 *     </li>
 *     <li>
 *         CN - genotype copy number (DEL/DUP/CNV only)
 *     </li>
 * </ul>
 *
 * <p>Note that for CNVs (DEL, DUP, multi-allelic CNV), GT alleles are set according to the CN/ECN fields. In
 * some cases, (e.g. diploid DUPs with CN 4), allele phasing cannot be determined unambiguously and GT is set with
 * no-call alleles.</p>
 *
 * <p>The tool generates a new VCF with clusters collapsed into single representative records. By default, a MEMBERS
 * field is generated that lists the input variant IDs contained in that record's cluster.</p>
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         One or more SV VCFs
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Clustered VCF
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk SVCluster \
 *       -V variants.vcf.gz \
 *       -O clustered.vcf.gz \
 *       --algorithm SINGLE_LINKAGE
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Stratifies structural variants by type, size, and context into mutually exclusive VCFs",
        oneLineSummary = "Stratifies structural variants by type, size, and context into mutually exclusive VCFs",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@BetaFeature
@DocumentedFeature
public final class SVStratify extends MultiVariantWalker {

    // Command-line arguments
    public static final String STRATIFY_CONFIG_FILE_LONG_NAME = "stratify-config";
    public static final String CONTEXT_NAME_FILE_LONG_NAME = "context-name";
    public static final String CONTEXT_INTERVAL_FILE_LONG_NAME = "context-intervals";
    public static final String OVERLAP_FRACTION_LONG_NAME = "overlap-fraction";
    public static final String NUM_BREAKPOINT_OVERLAPS_LONG_NAME = "num-breakpoint-overlaps";
    public static final String NUM_BREAKPOINT_INTERCHROM_OVERLAPS_LONG_NAME = "num-breakpoint-overlaps-interchromosomal";

    // Default output group suffix for unmatched records
    public static final String DEFAULT_STRATIFICATION = "non_matched";

    @Argument(
            doc = "Output directory",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private GATKPath outputPath;

    @Argument(
            doc = "Prefix for output filenames",
            fullName =  CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME
    )
    private String outputPrefix;

    /**
     * Expected format is tab-delimited and contains columns SVTYPE, MIN_SIZE, MAX_SIZE, CONTEXT. First line must be
     * a header with column names. Comment lines starting with {@link TableUtils#COMMENT_PREFIX} are ignored.
     */
    @Argument(
            doc = "Stratification configuration file (.tsv)",
            fullName = STRATIFY_CONFIG_FILE_LONG_NAME
    )
    private GATKPath configFile;

    @Argument(
            doc = "Context intervals file. Can be specified multiple times.",
            fullName = CONTEXT_INTERVAL_FILE_LONG_NAME,
            optional = true
    )
    private List<GATKPath> contextFileList;

    @Argument(
            doc = "Context names. Must be once for each --" + CONTEXT_INTERVAL_FILE_LONG_NAME,
            fullName = CONTEXT_NAME_FILE_LONG_NAME,
            optional = true
    )
    private List<String> contextNameList;

    @Argument(
            doc = "Minimum overlap fraction for contexts",
            minValue = 0,
            maxValue = 1,
            fullName = OVERLAP_FRACTION_LONG_NAME
    )
    private double overlapFraction = 0.5;

    @Argument(
            doc = "Minimum number of breakpoint overlaps for contexts",
            minValue = 0,
            maxValue = 2,
            fullName = NUM_BREAKPOINT_OVERLAPS_LONG_NAME
    )
    private int numBreakpointOverlaps = 0;

    @Argument(
            doc = "Minimum number of breakpoint overlaps for contexts for interchromosomal variants (e.g. BNDs)",
            minValue = 1,
            maxValue = 2,
            fullName = NUM_BREAKPOINT_INTERCHROM_OVERLAPS_LONG_NAME
    )
    private int numBreakpointOverlapsInterchrom = 1;

    protected SAMSequenceDictionary dictionary;
    protected Map<String, VariantContextWriter> writers;
    protected SVStatificationEngine engine;

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        dictionary = getMasterSequenceDictionary();
        Utils.validateArg(dictionary != null, "Reference dictionary is required; please specify with --" +
                StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME);
        loadStratificationConfig();
        logger.debug("Loaded stratification groups:");
        for (final SVStatificationEngine.SVStratification s : engine.getStratifications()) {
            logger.debug(s);
        }
        initializeWriters();
    }

    protected void createGroupWriter(final String name) {
        final String filename = outputPrefix + "." + name + ".vcf.gz";
        final Path path = outputPath.toPath().resolve(filename);
        final VariantContextWriter writer = createVCFWriter(path);
        writer.writeHeader(getHeaderForVariants());
        if (writers.containsKey(name)) {
            throw new GATKException.ShouldNeverReachHereException("Stratification name already exists: " + name);
        }
        writers.put(name, writer);
    }

    protected void initializeWriters() {
        writers = new HashMap<>();
        createGroupWriter(DEFAULT_STRATIFICATION);
        for (final SVStatificationEngine.SVStratification s : engine.getStratifications()) {
            createGroupWriter(s.getName());
        }
    }

    protected SVStatificationEngine loadStratificationConfig() {
        Utils.validateArg(contextNameList.size() == contextFileList.size(), "Arguments --" +
                CONTEXT_NAME_FILE_LONG_NAME + " and --" + CONTEXT_INTERVAL_FILE_LONG_NAME +
                " must be specified the same number of times.");
        final Map<String, List<Locatable>> map = new HashMap<>();
        final Iterator<String> nameIterator = contextNameList.iterator();
        final Iterator<GATKPath> pathIterator = contextFileList.iterator();
        final GenomeLocParser genomeLocParser = new GenomeLocParser(dictionary);
        while (nameIterator.hasNext() && pathIterator.hasNext()) {
            final String name = nameIterator.next();
            final GATKPath path = pathIterator.next();
            final GenomeLocSortedSet genomeLocs = IntervalUtils.loadIntervals(Collections.singletonList(path.toString()), IntervalSetRule.UNION, IntervalMergingRule.ALL, 0, genomeLocParser);
            final List<Locatable> intervals = Collections.unmodifiableList(genomeLocs.toList());
            if (map.containsKey(name)) {
                throw new UserException.BadInput("Duplicate context name was specified: " + name);
            }
            map.put(name, intervals);
        }
        return SVStatificationEngine.create(map, configFile);
    }


    @Override
    public void closeTool() {
        for (final VariantContextWriter writer : writers.values()) {
            writer.close();
        }
        super.closeTool();
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext,
                      final ReferenceContext referenceContext, final FeatureContext featureContext) {

        // Save a ton of compute by not copying genotypes into the new record
        final VariantContext variantNoGenotypes = new VariantContextBuilder(variant).genotypes(Collections.emptyList()).make();
        final SVCallRecord record = SVCallRecordUtils.create(variantNoGenotypes, dictionary);
        final SVStatificationEngine.SVStratification stratification = engine.getMatch(record, overlapFraction, numBreakpointOverlaps, numBreakpointOverlapsInterchrom);
        if (stratification == null) {
            writers.get(DEFAULT_STRATIFICATION).add(variant);
        } else {
            writers.get(stratification.getName()).add(variant);
        }
    }
}
