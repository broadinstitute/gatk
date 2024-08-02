package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
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
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVStratification;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.io.IOException;
import java.nio.file.Path;
import java.util.*;
import java.util.function.Function;

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

    // Default output group suffix for unmatched records
    public static final String DEFAULT_STRATIFICATION = "non_matched";

    // Configuration table column names
    public static final String SVTYPE_COLUMN = "SVTYPE";
    public static final String MIN_SIZE_COLUMN = "MIN_SIZE";
    public static final String MAX_SIZE_COLUMN = "MAX_SIZE";
    public static final String CONTEXT_COLUMN = "CONTEXT";

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
            fullName = CONTEXT_INTERVAL_FILE_LONG_NAME
    )
    private List<GATKPath> contextFileList;

    @Argument(
            doc = "Context names. Must be once for each --" + CONTEXT_INTERVAL_FILE_LONG_NAME,
            fullName = CONTEXT_NAME_FILE_LONG_NAME
    )
    private List<String> contextNameList;

    @Argument(
            doc = "Minimum overlap fraction for contexts",
            minValue = 0,
            maxValue = 1,
            fullName = OVERLAP_FRACTION_LONG_NAME
    )
    private double overlapFraction;

    @Argument(
            doc = "Minimum number of breakpoint overlaps for contexts",
            minValue = 0,
            maxValue = 2,
            fullName = NUM_BREAKPOINT_OVERLAPS_LONG_NAME
    )
    private int numBreakpointOverlaps;

    protected SAMSequenceDictionary dictionary;
    protected Map<String, List<Locatable>> contextMap;
    protected List<SVStratification> strats;
    protected Map<String, VariantContextWriter> writers;

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        dictionary = getSequenceDictionaryForDrivingVariants();
        Utils.validateArg(dictionary != null, "Reference dictionary is required in the input VCF");
        contextMap = loadContextIntervals();
        strats = loadConfigurationTable();
        logger.debug("Loaded stratification groups:");
        for (final SVStratification s : strats) {
            logger.debug(s);
        }
        initializeWriters();
    }

    protected void createGroupWriter(final String name) {
        final String filename = outputPrefix + "." + name + ".vcf.gz";
        final Path path = outputPath.toPath().resolve(filename);
        final VariantContextWriter writer = createVCFWriter(path);
        writer.setHeader(getHeaderForVariants());
        if (writers.containsKey(name)) {
            throw new GATKException.ShouldNeverReachHereException("Stratification name already exists: " + name);
        }
        writers.put(name, writer);
    }

    protected void initializeWriters() {
        createGroupWriter(DEFAULT_STRATIFICATION);
        for (final SVStratification s : strats) {
            createGroupWriter(s.getName());
        }
    }

    protected Map<String, List<Locatable>> loadContextIntervals() {
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
            map.put(name, intervals);
        }
        return map;
    }

    protected List<SVStratification> loadConfigurationTable() {
        final List<SVStratification> configList;
        try (final TableReader<SVStratification> tableReader = TableUtils.reader(configFile.toPath(), this::tableParser)) {
            // force into array list to ensure fast index lookup
            configList = new ArrayList<>(tableReader.toList());
        } catch (final IOException e) {
            throw new GATKException("IO error while reading config table", e);
        }
        for (int i = 0; i < configList.size(); i++) {
            for (int j = i + 1; j < configList.size(); j++) {
                if (!configList.get(i).isMutuallyExclusive(configList.get(j))) {
                    throw new UserException.BadInput("Configuration records " + configList.get(i) + " and "
                            + configList.get(j) + " are not mutually exclusive.");
                }
            }
        }
        return configList;
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
        for (final SVStratification s : strats) {
            if (s.matches(variant, overlapFraction, numBreakpointOverlaps, dictionary)) {
                writers.get(s.getName()).add(variant);
                return;
            }
        }
        writers.get(DEFAULT_STRATIFICATION).add(variant);
    }

    public Function<DataLine, SVStratification> tableParser(TableColumnCollection columns, Function<String, RuntimeException> exceptionFactory) {
        if (columns.columnCount() != 4) {
            throw exceptionFactory.apply("Expected 4 columns but found " + columns.columnCount());
        }
        if (!columns.contains(SVTYPE_COLUMN)) {
            throw exceptionFactory.apply("Missing column " + SVTYPE_COLUMN);
        }
        if (!columns.contains(MIN_SIZE_COLUMN)) {
            throw exceptionFactory.apply("Missing column " + MIN_SIZE_COLUMN);
        }
        if (!columns.contains(MAX_SIZE_COLUMN)) {
            throw exceptionFactory.apply("Missing column " + MAX_SIZE_COLUMN);
        }
        if (!columns.contains(CONTEXT_COLUMN)) {
            throw exceptionFactory.apply("Missing column " + CONTEXT_COLUMN);
        }
        return this::parseTableLine;
    }

    protected SVStratification parseTableLine(final DataLine dataLine) {
        final GATKSVVCFConstants.StructuralVariantAnnotationType svType = GATKSVVCFConstants.StructuralVariantAnnotationType.valueOf(dataLine.get(SVTYPE_COLUMN));
        Integer minSize = Integer.valueOf(dataLine.get(MIN_SIZE_COLUMN));
        Integer maxSize = Integer.valueOf(dataLine.get(MAX_SIZE_COLUMN));
        final String context = dataLine.get(CONTEXT_COLUMN);
        if (!contextMap.containsKey(context)) {
            throw new GATKException("Could not find context with name " + context);
        }
        return new SVStratification(svType, minSize, maxSize, context, contextMap.get(context));
    }
}
