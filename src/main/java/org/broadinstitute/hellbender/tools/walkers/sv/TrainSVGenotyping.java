package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.arrow.util.VisibleForTesting;
import org.apache.commons.io.IOUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.*;
import org.broadinstitute.hellbender.tools.sv.aggregation.DiscordantPairEvidenceAggregator;
import org.broadinstitute.hellbender.tools.sv.aggregation.PESREvidenceTester;
import org.broadinstitute.hellbender.tools.sv.aggregation.SplitReadEvidenceAggregator;
import org.broadinstitute.hellbender.tools.sv.cluster.PloidyTable;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterWalker;
import org.broadinstitute.hellbender.tools.sv.stratify.SVStratificationEngine;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.codecs.DepthEvidenceCodec;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.tsv.*;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.StreamSupport;

import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.CountingVariantFilter;

/**
 * <p>Trains SV genotyping models.</p>
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         SV VCF
 *     </li>
 *     <li>
 *         CNV training intervals
 *     </li>
 *     <li>
 *         Depth evidence file
 *     </li>
 *     <li>
 *         Discordant pairs evidence file
 *     </li>
 *     <li>
 *         Split reads evidence file
 *     </li>
 *     <li>
 *         Median binned read counts table
 *     </li>
 *     <li>
 *         Reference sequence dictionary
 *     </li>
 *     <li>
 *         Sample ploidy table
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         SV VCF (all records by default; if PE/SR training downsampling is activated because the number of
 *         eligible training records exceeds the configured cap, the VCF will contain only the retained training subset)
 *     </li>
 *     <li>
 *         Genotyping cutoff tables (RD, PE, SR)
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk TrainSVGenotyping \
 *       -V variants.vcf.gz \
 *       -O variants.genotyped.vcf.gz \
 *       --output-dir ./ \
 *       --output-name my_batch \
 *       --training-intervals train.bed \
 *       --median-coverage coverage.tsv \
 *       --rd-file batch.rd.txt.gz \
 *       --discordant-pairs-file batch.pe.txt.gz \
 *       --split-reads-file batch.sr.txt.gz \
 *       --sequence-dictionary /Homo_sapiens_assembly38.dict \
 *       --ploidy-table batch.ploidy_table.tsv
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        summary = "Trains SV genotyping models",
        oneLineSummary = "Trains SV genotyping models",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
public final class TrainSVGenotyping extends MultiplePassVariantWalker {

    public static final String DEPTH_EVIDENCE_FILE_PATH_LONG_NAME = "rd-file";
    public static final String MEDIAN_COUNTS_FILE_PATH_LONG_NAME = AggregateSVEvidence.MEDIAN_COVERAGE_LONG_NAME;
    public static final String TABLES_DIR_LONG_NAME = "output-dir";
    public static final String TABLES_NAME_LONG_NAME = "output-name";
    public static final String NUM_TRAINING_STATES_LONG_NAME = "n-training-states";
    public static final String TRAINING_INTERVALS_LONG_NAME = "training-intervals";
    public static final String PESR_EXCLUSION_INTERVALS_LONG_NAME = "pesr-exclusion-intervals";
    public static final String DEPTH_EXCLUSION_INTERVALS_LONG_NAME = "depth-exclusion-intervals";
    public static final String MIN_PE_QUALITY_LONG_NAME = "pe-quality";
    public static final String MIN_PESER_SIZE_LONG_NAME = "min-pesr-size";
    public static final String MIN_SR_QUALITY_LONG_NAME = "sr-quality";
    public static final String DISCORDANT_PAIR_QUERY_LOOKAHEAD_LONG_NAME = "pe-query-lookahead";
    public static final String SPLIT_READ_QUERY_LOOKAHEAD_LONG_NAME = "sr-query-lookahead";
    public static final String MAX_TRAINING_RECORDS_LONG_NAME = "max-training-records";

    @Argument(
            fullName = DEPTH_EVIDENCE_FILE_PATH_LONG_NAME,
            doc = "Indexed read counts file ending with " + DepthEvidenceCodec.FORMAT_SUFFIX + ".gz"
    )
    public GATKPath depthEvidenceFile;

    @Argument(
            fullName = MEDIAN_COUNTS_FILE_PATH_LONG_NAME,
            doc = "Median counts file"
    )
    public GATKPath medianFile;

    @Argument(
            doc = "Discordant pairs evidence file (indexed and ending in .pe.txt.gz)",
            fullName = AggregateSVEvidence.DISCORDANT_PAIRS_LONG_NAME,
            optional = true
    )
    private GATKPath discordantPairsFile;

    @Argument(
            doc = "Split reads evidence file (indexed and ending in .sr.txt.gz)",
            fullName = AggregateSVEvidence.SPLIT_READ_LONG_NAME,
            optional = true
    )
    private GATKPath splitReadsFile;

    /**
     * Expected format is tab-delimited and contains a header with the first column SAMPLE and remaining columns
     * contig names. Each row corresponds to a sample, with the sample ID in the first column and contig ploidy
     * integers in their respective columns.
     */
    @Argument(
            doc = "Sample ploidy table (.tsv)",
            fullName = SVClusterWalker.PLOIDY_TABLE_LONG_NAME
    )
    protected GATKPath ploidyTablePath;

    @Argument(
            doc = "Training site intervals file",
            fullName = TRAINING_INTERVALS_LONG_NAME,
            optional = true
    )
    protected GATKPath trainingIntervalsPath;

    @Argument(
            doc = "Intervals to exclude for RD evidence genotype training",
            fullName = DEPTH_EXCLUSION_INTERVALS_LONG_NAME,
            optional = true
    )
    protected GATKPath depthExclusionIntervalsPath;

    @Argument(
            doc = "Intervals to exclude for PE/SR evidence genotype training",
            fullName = PESR_EXCLUSION_INTERVALS_LONG_NAME,
            optional = true
    )
    protected GATKPath pesrExclusionIntervalsPath;

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output VCF"
    )
    public GATKPath outputVcf;

    @Argument(
            fullName = TABLES_DIR_LONG_NAME,
            doc = "Metric table output directory"
    )
    public GATKPath tablesDir;

    @Argument(
            fullName = TABLES_NAME_LONG_NAME,
            doc = "Metric table output file base name"
    )
    public String tableBaseName;

    @Argument(
            fullName = AggregateDepthEvidence.LARGE_VARIANT_SIZE_LONG_NAME,
            doc = "Large variant size",
            minValue = 0
    )
    public long largeVariantSize = 2500000L;

    @Argument(
            fullName = AggregateDepthEvidence.LARGE_VARIANT_POINTS_LONG_NAME,
            doc = "Large variant points",
            minValue = 1
    )
    public int largeVariantPoints = 500;

    @Argument(
            fullName = AggregateDepthEvidence.LARGE_VARIANT_WINDOW_LONG_NAME,
            doc = "Large variant window",
            minValue = 1
    )
    public int largeVariantWindow = 2000;

    @Argument(
            fullName = AggregateDepthEvidence.NUM_BINS_LONG_NAME,
            doc = "Number of bins to resample to",
            minValue = 1
    )
    public int numBins = 10;

    @Argument(
            fullName = MIN_PE_QUALITY_LONG_NAME,
            doc = "Discordant pair quality threshold, used for setting evidence count threshold",
            minValue = 1
    )
    public double minDiscordantPairQuality = 30;

    @Argument(
            fullName = MIN_PESER_SIZE_LONG_NAME,
            doc = "Discordant pair and split read training minimum size",
            minValue = 0
    )
    public int minPesrSize = 1000;

    @Argument(
            fullName = MIN_SR_QUALITY_LONG_NAME,
            doc = "Split read quality threshold, used for setting evidence count threshold",
            minValue = 1
    )
    public double minSplitReadQuality = 30;

        @Argument(
            fullName = DISCORDANT_PAIR_QUERY_LOOKAHEAD_LONG_NAME,
            doc = "Number of bases to prefetch after PE evidence query cache misses",
            minValue = 0,
            optional = true
        )
        public int discordantPairQueryLookahead = FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES;

        @Argument(
            fullName = SPLIT_READ_QUERY_LOOKAHEAD_LONG_NAME,
            doc = "Number of bases to prefetch after SR evidence query cache misses",
            minValue = 0,
            optional = true
        )
        public int splitReadQueryLookahead = FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES;

        @Argument(
            fullName = MAX_TRAINING_RECORDS_LONG_NAME,
            doc = "Maximum number of eligible PE/SR training records to use before deterministic every-nth downsampling is applied. Set to 0 to disable the cap.",
            minValue = 0,
            optional = true
        )
        public int maxTrainingRecords = 100_000;

    @Argument(
            fullName = AggregateDepthEvidence.MAX_QUALITY_LONG_NAME,
            doc = "Max quality score",
            minValue = 1
    )
    public int maxQual = 999;

    @Argument(
            fullName = NUM_TRAINING_STATES_LONG_NAME,
            doc = "Number of training copy states",
            minValue = 2
    )
    public int numTrainingStates = 5;

    /**
     * Paired-end window size downstream of the first mate and upstream of the second mate.
     */
    @Argument(
            doc = "Inner discordant pair window size (bp)",
            fullName = AggregateSVEvidence.PE_INNER_WINDOW_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private int innerWindow = 50;

    /**
     * Paired-end window size upstream of the first mate and downstream of the second mate.
     */
    @Argument(
            doc = "Outer discordant pair window size (bp)",
            fullName = AggregateSVEvidence.PE_OUTER_WINDOW_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private int outerWindow = 500;

    private SAMSequenceDictionary dictionary;
    private VariantContextWriter writer;
    private VCFHeader outputHeader;
    private Map<String, Double> sampleMedians;
    private PloidyTable ploidyTable;
    private DepthMatrixLoader loader;
    private DepthEvidenceGenotyper depthGenotyper;
    private List<String> masterSampleList;
    private FeatureDataSource<DepthEvidence> depthSource;
    private List<DepthEvidenceGenotyper.CopyStateStats> trainedCopyStateStats;

    private Map<String, DepthEvidenceGenotyper.DepthGenotypeResult> depthGenotypeResults;
    private Map<String, DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult> discordantPairGenotypeResults;
    private Map<String, SplitReadEvidenceGenotyper.SplitReadGenotypeResult> splitReadGenotypeResults;

    private FeatureDataSource<DiscordantPairEvidence> discordantPairSource;
    private DiscordantPairEvidenceAggregator discordantPairCollector;
    private DiscordantPairEvidenceGenotyper discordantPairGenotyper;
    private DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeParameters discordantPairParameters;

    private FeatureDataSource<SplitReadEvidence> splitReadStartSource;
    private FeatureDataSource<SplitReadEvidence> splitReadEndSource;
    private SplitReadEvidenceAggregator splitReadStartCollector;
    private SplitReadEvidenceAggregator splitReadEndCollector;
    private SplitReadEvidenceGenotyper splitReadGenotyper;
    private SplitReadEvidenceGenotyper.SplitReadGenotypeParameters splitReadParameters;
    private SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs splitReadFrequencyCutoffs;

    private OverlapDetector<SimpleInterval> depthExclusionIntervals;
    private SVStratificationEngine pesrExclusionEngine;
    private Set<String> selectedDiscordantPairTrainingRecordIds = Collections.emptySet();
    private Set<String> selectedSplitReadTrainingRecordIds = Collections.emptySet();
    private boolean emitTrainingSubsetOnly = false;

    private static final String PESR_EXCLUSION_STRATIFICATION = "pesrex";

    protected int numberOfPasses() {
        // Note: the actual traversal is managed by the overridden traverse() below,
        // which caches records and merges some passes. This value is retained for
        // compatibility with the parent contract but is not used to control traversal.
        return 8;
    }

    /**
     * Overrides the default multi-pass traversal to:
     * 1. Read the VCF only once (pass 0) and cache SVCallRecords in memory,
     *    then iterate the cached list for all subsequent passes. This eliminates
     *    7 redundant VCF decompressions and SVCallRecord re-creations.
     * 2. Merge old passes 3 (PE genotyping) and 4 (SR first pass) into a single
     *    traversal, since SR first pass can consume PE genotype results as they
     *    are produced within the same record iteration.
     */
    @Override
    public void traverse() {
        final CountingVariantFilter countingVariantFilter = makeVariantFilter();
        final CountingReadFilter readFilter = makeReadFilter();

        // Pass 0: Full VCF read — depth genotyping, PE overlap registration, and record caching
        logger.info("Starting pass 0 - depth genotyping (caching all records)");
        final List<SVCallRecord> cachedRecords = new ArrayList<>();
        StreamSupport.stream(getSpliteratorForDrivingVariants(), false)
                .filter(countingVariantFilter)
                .forEach(variant -> {
                    final SVCallRecord record = SVCallRecordUtils.create(variant, dictionary);
                    cachedRecords.add(record);
                    applyReadDepth(record);
                    if (discordantPairCollectionEnabled()) {
                        discordantPairGenotyper.registerVariantForOverlapCheck(record);
                    }
                    progressMeter.update(new SimpleInterval(variant));
                });
        logger.info("Finished pass 0 (" + cachedRecords.size() + " variants cached)");
        afterNthPass(0);
        initializeDiscordantPairTrainingSelection(cachedRecords);
        progressMeter.reset();

        // Pass 1: PE evidence collection (training sites only)
        logger.info("Starting pass 1 - discordant pair evidence collection");
        for (final SVCallRecord record : cachedRecords) {
            applyDiscordantPairFirstPass(record);
        }
        logger.info("Finished pass 1");
        afterNthPass(1);
        initializeSplitReadTrainingSelection(cachedRecords);
        progressMeter.reset();

        // Pass 2: PE parameter estimation
        logger.info("Starting pass 2 - discordant pair parameter estimation");
        for (final SVCallRecord record : cachedRecords) {
            applyDiscordantPairSecondPass(record);
        }
        logger.info("Finished pass 2");
        afterNthPass(2);
        progressMeter.reset();

        // Merged pass 3+4: PE genotyping + SR evidence collection in a single traversal.
        // This is safe because SR first pass only checks discordantPairGenotypeResults for
        // the current record, which is populated by applyDiscordantPairThirdPass just above.
        logger.info("Starting pass 3 - discordant pair genotyping + split read evidence collection (merged)");
        for (final SVCallRecord record : cachedRecords) {
            applyDiscordantPairThirdPass(record);
            applySplitReadFirstPass(record);
        }
        logger.info("Finished pass 3");
        afterNthPass(3);
        afterNthPass(4);
        // PE training data is no longer needed after SR first pass finalization
        if (discordantPairCollectionEnabled()) {
            discordantPairGenotyper.clearTrainingData();
        }
        progressMeter.reset();

        // Pass 5 (old pass 5): SR parameter estimation
        logger.info("Starting pass 4 - split read parameter estimation");
        for (final SVCallRecord record : cachedRecords) {
            applySplitReadSecondPass(record);
        }
        logger.info("Finished pass 4");
        afterNthPass(5);
        progressMeter.reset();

        // Pass 6 (old pass 6): SR genotyping
        logger.info("Starting pass 5 - split read genotyping");
        for (final SVCallRecord record : cachedRecords) {
            applySplitReadThirdPass(record);
        }
        logger.info("Finished pass 5");
        afterNthPass(6);
        progressMeter.reset();

        // Pass 7 (old pass 7): Write final genotypes
        logger.info("Starting pass 6 - writing genotypes");
        for (final SVCallRecord record : cachedRecords) {
            writeGenotypes(record);
        }
        logger.info("Finished writing genotypes");

        logger.info(countingVariantFilter.getSummaryLine());
        logger.info(readFilter.getSummaryLine());
    }

    private void initializeDiscordantPairTrainingSelection(final List<SVCallRecord> cachedRecords) {
        if (!discordantPairCollectionEnabled()) {
            selectedDiscordantPairTrainingRecordIds = Collections.emptySet();
            return;
        }

        final LinkedHashSet<String> eligibleIds = new LinkedHashSet<>();
        for (final SVCallRecord record : cachedRecords) {
            final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = depthGenotypeResults.get(record.getId());
            if (depthResult != null && discordantPairGenotyper.trainableRecord(record, depthResult, pesrExclusionEngine)) {
                eligibleIds.add(record.getId());
            }
        }
        selectedDiscordantPairTrainingRecordIds = downsampleEligibleRecords(eligibleIds, "PE");
    }

    private void initializeSplitReadTrainingSelection(final List<SVCallRecord> cachedRecords) {
        if (!splitReadCollectionEnabled()) {
            selectedSplitReadTrainingRecordIds = Collections.emptySet();
            return;
        }

        final LinkedHashSet<String> eligibleIds = new LinkedHashSet<>();
        for (final SVCallRecord record : cachedRecords) {
            if (selectedDiscordantPairTrainingRecordIds.contains(record.getId())
                    && splitReadGenotyper.trainableRecord(record, discordantPairGenotyper, pesrExclusionEngine)) {
                eligibleIds.add(record.getId());
            }
        }
        selectedSplitReadTrainingRecordIds = downsampleEligibleRecords(eligibleIds, "SR");
    }

    private Set<String> downsampleEligibleRecords(final LinkedHashSet<String> eligibleIds, final String label) {
        if (eligibleIds.isEmpty()) {
            logger.info("No eligible " + label + " training records found");
            return Collections.emptySet();
        }
        if (maxTrainingRecords == 0 || eligibleIds.size() <= maxTrainingRecords) {
            logger.info("Using all " + eligibleIds.size() + " eligible " + label + " training records");
            return eligibleIds;
        }

        final LinkedHashSet<String> selectedIds = selectEveryNthEligibleRecords(eligibleIds, maxTrainingRecords);
        final int stride = (int) Math.ceil(eligibleIds.size() / (double) maxTrainingRecords);
        emitTrainingSubsetOnly = true;
        logger.info("Downsampled " + label + " training records from " + eligibleIds.size() + " to " + selectedIds.size() + " using stride " + stride);
        return selectedIds;
    }

    @VisibleForTesting
    static LinkedHashSet<String> selectEveryNthEligibleRecords(final LinkedHashSet<String> eligibleIds, final int maxTrainingRecords) {
        Utils.nonNull(eligibleIds);
        Utils.validateArg(maxTrainingRecords > 0, "maxTrainingRecords must be greater than zero when downsampling");
        final int stride = (int) Math.ceil(eligibleIds.size() / (double) maxTrainingRecords);
        final LinkedHashSet<String> selectedIds = new LinkedHashSet<>();
        int index = 0;
        for (final String recordId : eligibleIds) {
            if (index % stride == 0) {
                selectedIds.add(recordId);
            }
            index++;
        }
        return selectedIds;
    }

    private boolean retainDiscordantPairTrainingRecord(final SVCallRecord record) {
        return !discordantPairCollectionEnabled() || selectedDiscordantPairTrainingRecordIds.contains(record.getId());
    }

    private boolean retainSplitReadTrainingRecord(final SVCallRecord record) {
        return !splitReadCollectionEnabled() || selectedSplitReadTrainingRecordIds.contains(record.getId());
    }

    private boolean shouldWriteRecord(final SVCallRecord record) {
        return !emitTrainingSubsetOnly || !discordantPairCollectionEnabled() || selectedDiscordantPairTrainingRecordIds.contains(record.getId());
    }

    private void initializeDiscordantPairCollection() {
        if (discordantPairCollectionEnabled()) {
            discordantPairSource = new FeatureDataSource<>(
                    discordantPairsFile.toString(),
                    "discordantPairsFile",
                    discordantPairQueryLookahead,
                    DiscordantPairEvidence.class,
                    cloudPrefetchBuffer,
                    cloudIndexPrefetchBuffer);
            discordantPairCollector = new DiscordantPairEvidenceAggregator(discordantPairSource, dictionary, innerWindow, outerWindow);
        }
    }

    private void initializeSplitReadCollection() {
        if (splitReadCollectionEnabled()) {
            splitReadStartSource = new FeatureDataSource<>(
                    splitReadsFile.toString(),
                    "splitReadsStartFile",
                    splitReadQueryLookahead,
                    SplitReadEvidence.class,
                    cloudPrefetchBuffer,
                    cloudIndexPrefetchBuffer);
            splitReadEndSource = new FeatureDataSource<>(
                    splitReadsFile.toString(),
                    "splitReadsEndFile",
                    splitReadQueryLookahead,
                    SplitReadEvidence.class,
                    cloudPrefetchBuffer,
                    cloudIndexPrefetchBuffer);
            splitReadStartCollector = new SplitReadEvidenceAggregator(splitReadStartSource, dictionary, 0, true);
            splitReadEndCollector = new SplitReadEvidenceAggregator(splitReadEndSource, dictionary, 0, false);
        }
    }

    @Override
    public void onTraversalStart() {
        depthSource = new FeatureDataSource<>(depthEvidenceFile.toString());
        sampleMedians = loadMedianSampleCoverageTable();
        dictionary = getBestAvailableSequenceDictionary();

        if (depthExclusionIntervalsPath != null) {
            final GenomeLocParser genomeLocParser = new GenomeLocParser(dictionary);
            final GenomeLocSortedSet genomeLocs = IntervalUtils.loadIntervals(Collections.singletonList(depthExclusionIntervalsPath.toString()), IntervalSetRule.UNION, IntervalMergingRule.OVERLAPPING_ONLY, 0, genomeLocParser);
            depthExclusionIntervals = OverlapDetector.create(genomeLocs.stream().map(SimpleInterval::new).toList());
        }
        if (pesrExclusionIntervalsPath != null) {
            pesrExclusionEngine = new SVStratificationEngine(dictionary);
            final GenomeLocParser genomeLocParser = new GenomeLocParser(dictionary);
            final GenomeLocSortedSet genomeLocs = IntervalUtils.loadIntervals(Collections.singletonList(pesrExclusionIntervalsPath.toString()), IntervalSetRule.UNION, IntervalMergingRule.OVERLAPPING_ONLY, 0, genomeLocParser);
            final List<Locatable> intervals = Collections.unmodifiableList(genomeLocs.toList());
            pesrExclusionEngine.addTrack(PESR_EXCLUSION_STRATIFICATION, intervals);
            pesrExclusionEngine.addStratification(PESR_EXCLUSION_STRATIFICATION, null, null, null, Collections.singleton(PESR_EXCLUSION_STRATIFICATION));
        }

        loader = new DepthMatrixLoader(depthSource, numBins, largeVariantSize, largeVariantPoints, largeVariantWindow, depthExclusionIntervals, dictionary);
        writer = createVCFWriter(outputVcf);
        outputHeader = createHeader(getHeaderForVariants());
        writer.writeHeader(outputHeader);
        masterSampleList = outputHeader.getSampleNamesInOrder();
        depthGenotyper = new DepthEvidenceGenotyper(null, masterSampleList, maxQual, dictionary);
        trainCopyNumberSites();

        depthGenotypeResults = new HashMap<>();
        discordantPairGenotypeResults = new HashMap<>();
        splitReadGenotypeResults = new HashMap<>();

        if (splitReadCollectionEnabled()) {
            Utils.validate(discordantPairCollectionEnabled(), "Discordant pairs file must be provided for split read training");
        }

        initializeDiscordantPairCollection();
        initializeSplitReadCollection();
        discordantPairGenotyper = new DiscordantPairEvidenceGenotyper(sampleMedians, minDiscordantPairQuality, minPesrSize, PESREvidenceTester.DEPTH_BASIS, maxQual);
        splitReadGenotyper = new SplitReadEvidenceGenotyper(sampleMedians, masterSampleList.size(), minSplitReadQuality, minPesrSize, PESREvidenceTester.DEPTH_BASIS, maxQual);
    }

    private void trainCopyNumberSites() {
        final GenomeLocParser genomeLocParser = new GenomeLocParser(dictionary);
        final GenomeLocSortedSet trainingLocs = IntervalUtils.loadIntervals(Collections.singletonList(trainingIntervalsPath.toString()), IntervalSetRule.UNION, IntervalMergingRule.OVERLAPPING_ONLY, 0, genomeLocParser);
        logger.info("Training on " + trainingLocs.size() + " CNV sites");
        final List<DepthEvidenceGenotyper.DepthGenotypeResult> genotypeResults = new ArrayList<>();
        for (final GenomeLoc genomeLoc : trainingLocs) {
            final DepthMatrix depthMatrix = loader.load(new SimpleInterval(genomeLoc), sampleMedians);
            genotypeResults.add(depthGenotyper.genotype(depthMatrix));
        }
        trainedCopyStateStats = depthGenotyper.train(genotypeResults, numTrainingStates);
        try (final TableWriter<DepthEvidenceGenotyper.CopyStateStats> tableWriter = TableUtils.writer(getTablePath(".rd_geno_params.tsv").toPath(), DepthEvidenceGenotyper.DepthTableParser.CUTOFFS_COLUMNS, new DepthEvidenceGenotyper.DepthTableParser()::composeCutoffsLine)) {
            tableWriter.writeAllRecords(trainedCopyStateStats);
        } catch (IOException e) {
            throw new GATKException("Error while writing RD cutoffs table", e);
        }
        logger.info("Training completed");
    }

    @VisibleForTesting
    protected Map<String, Double> loadMedianSampleCoverageTable() {
        final List<String> lines = IOUtils.readLines(BucketUtils.openFile(medianFile.toString()), Charset.defaultCharset());
        Utils.validate(lines.size() >= 2, "Median coverage file must contain at least two lines");
        final String[] samples = lines.get(0).split("\t");
        final String[] values = lines.get(1).split("\t");
        Utils.validate(samples.length == values.length,
                "Median file's first two lines must have the same number of columns");
        final Map<String, Double> sampleMedians = new HashMap<>();
        try {
            for (int i = 0; i < samples.length; i++) {
                sampleMedians.put(samples[i], Double.valueOf(values[i]));
            }
        } catch (final NumberFormatException nfe) {
            throw new UserException.BadInput(nfe.getMessage());
        }
        final List<String> vcfSamples = getHeaderForVariants().getSampleNamesInOrder();
        Utils.validate(sampleMedians.keySet().containsAll(vcfSamples), "Median counts table does not contain all samples in the VCF");
        return sampleMedians;
    }

    @Override
    public void nthPassApply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext, final int n) {
        final SVCallRecord record = SVCallRecordUtils.create(variant, dictionary);
        if (n == 0) {
            applyReadDepth(record);
            if (discordantPairCollectionEnabled()) {
                discordantPairGenotyper.registerVariantForOverlapCheck(record);
            }
        } else if (n == 1) {
            applyDiscordantPairFirstPass(record);
        } else if (n == 2) {
            applyDiscordantPairSecondPass(record);
        } else if (n == 3) {
            applyDiscordantPairThirdPass(record);
        } else if (n == 4) {
            applySplitReadFirstPass(record);
        } else if (n == 5) {
            applySplitReadSecondPass(record);
        } else if (n == 6) {
            applySplitReadThirdPass(record);
        } else if (n == 7) {
            writeGenotypes(record);
        } else {
            throw new GATKException("Unexpected number of passes: " + n);
        }
    }

    public void writeGenotypes(final SVCallRecord record) {
        if (!shouldWriteRecord(record)) {
            return;
        }
        final ArrayList<Genotype> newGenotypeList = new ArrayList<>(masterSampleList.size());
        // TODO: need to decide on whether to support genotyped VCF input and make consistent with GenotypeSVs
        final GenotypesContext genotypes = SVCallRecordUtils.populateGenotypesForMissingSamplesWithAlleles(record, new HashSet<>(masterSampleList), false, ploidyTable, outputHeader);
        for (int i = 0; i < masterSampleList.size(); i++) {
            final String sample = masterSampleList.get(i);
            if (!genotypes.containsSample(sample)) {
                throw new IllegalArgumentException("Sample " + sample + " does not exist in record " + record.getId());
            }
            final GenotypeBuilder builder = new GenotypeBuilder(genotypes.get(sample));
            if (depthGenotypeResults.containsKey(record.getId())) {
                final DepthEvidenceGenotyper.DepthGenotypeResult result = depthGenotypeResults.get(record.getId());
                builder.attribute(GATKSVVCFConstants.DEPTH_GENOTYPE_COPY_NUMBER_FORMAT, result.copyStates()[i]);
                builder.attribute(GATKSVVCFConstants.DEPTH_MEDIAN_COPY_RATIO, result.sampleDepths()[i]);
                builder.attribute(GATKSVVCFConstants.DEPTH_GENOTYPE_QUALITY_ATTRIBUTE, (int) result.genotypeQuals()[i]);
            }
            if (discordantPairGenotypeResults.containsKey(record.getId())) {
                final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult result = discordantPairGenotypeResults.get(record.getId());
                builder.attribute(GATKSVVCFConstants.DISCORDANT_PAIR_GENOTYPE_ATTRIBUTE, result.genotypes()[i]);
                builder.attribute(GATKSVVCFConstants.DISCORDANT_PAIR_GENOTYPE_QUALITY_ATTRIBUTE, result.genotypeQuals()[i]);
            }
            if (splitReadGenotypeResults.containsKey(record.getId())) {
                final SplitReadEvidenceGenotyper.SplitReadGenotypeResult result = splitReadGenotypeResults.get(record.getId());
                builder.attribute(GATKSVVCFConstants.SPLIT_READ_GENOTYPE_ATTRIBUTE, result.genotypes()[i]);
                if (result.genotypeQuals() != null) {
                    builder.attribute(GATKSVVCFConstants.SPLIT_READ_GENOTYPE_QUALITY_ATTRIBUTE, result.genotypeQuals()[i]);
                }
            }
            newGenotypeList.add(builder.make());
        }
        final GenotypesContext newGenotypes = GenotypesContext.create(newGenotypeList);
        final SVCallRecord regenotypedCall = SVCallRecordUtils.copyCallWithNewGenotypes(record, newGenotypes);
        final VariantContext variant = SVCallRecordUtils.getVariantBuilder(regenotypedCall).make();
        writer.add(variant);
    }

    @Override
    protected void afterNthPass(final int n) {
        if (n == 0 && discordantPairCollectionEnabled()) {
            discordantPairGenotyper.aggregateOverlapCheckIntervals();
        } else if (n == 1 && discordantPairCollectionEnabled()) {
            discordantPairGenotyper.finalizeFirstPass();
        } else if (n == 2 && discordantPairCollectionEnabled()) {
            discordantPairParameters = discordantPairGenotyper.finalizeSecondPass();
            writeDiscordantPairParameters(discordantPairParameters);
        } else if (n == 4 && splitReadCollectionEnabled()) {
            splitReadGenotyper.finalizeFirstPass();
        } else if (n == 5 && splitReadCollectionEnabled()) {
            splitReadParameters = splitReadGenotyper.finalizeSecondPass();
        } else if (n == 6 && splitReadCollectionEnabled()) {
            splitReadFrequencyCutoffs = splitReadGenotyper.finalizeThirdPass();
            writeSplitReadParameters(new SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics(splitReadParameters, splitReadFrequencyCutoffs));
        }
    }

    private void writeDiscordantPairParameters(final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeParameters parameters) {
        try (final TableWriter<DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeParameters> tableWriter = TableUtils.writer(getTablePath(".pe_geno_params.tsv").toPath(), DiscordantPairEvidenceGenotyper.DiscordantPairTableParser.CUTOFFS_COLUMNS, new DiscordantPairEvidenceGenotyper.DiscordantPairTableParser()::composeCutoffsLine)) {
            tableWriter.writeRecord(parameters);
        } catch (IOException e) {
            throw new GATKException("Error while writing PE cutoffs table", e);
        }
    }

    private void writeSplitReadParameters(final SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics parameters) {
        try (final TableWriter<SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics> tableWriter = TableUtils.writer(getTablePath(".sr_geno_params.tsv").toPath(), SplitReadEvidenceGenotyper.SplitReadTableParser.CUTOFFS_COLUMNS, new SplitReadEvidenceGenotyper.SplitReadTableParser()::composeCutoffsLine)) {
            tableWriter.writeRecord(parameters);
        } catch (IOException e) {
            throw new GATKException("Error while writing PE cutoffs table", e);
        }
    }

    private void applyReadDepth(final SVCallRecord record) {
        // Must be a CNV
        final GATKSVVCFConstants.StructuralVariantAnnotationType svtype = record.getType();
        if (svtype != GATKSVVCFConstants.StructuralVariantAnnotationType.DEL && svtype != GATKSVVCFConstants.StructuralVariantAnnotationType.DUP && svtype != GATKSVVCFConstants.StructuralVariantAnnotationType.CNV) {
            return;
        }
        final DepthMatrix depthMatrix = loader.load(new SimpleInterval(record.getContigA(), record.getPositionA(), record.getPositionB()), sampleMedians);

        final DepthEvidenceGenotyper.DepthGenotypeResult genotypeResult = depthGenotyper.genotype(depthMatrix);
        if (depthGenotypeResults.containsKey(record.getContigA())) {
            throw new UserException.BadInput("Duplicate variant ID: " + record.getId());
        }
        if (genotypeResult != null) {
            depthGenotypeResults.put(record.getId(), genotypeResult);
        }
    }

    private void applyDiscordantPairFirstPass(final SVCallRecord record) {
        if (discordantPairCollectionEnabled()
                && depthGenotypeResults.containsKey(record.getId())
                && retainDiscordantPairTrainingRecord(record)) {
            final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = depthGenotypeResults.get(record.getId());
            if (discordantPairGenotyper.trainableRecord(record, depthResult, pesrExclusionEngine)) {
                final List<DiscordantPairEvidence> discordantPairEvidence = discordantPairCollector.collectEvidence(record);
                discordantPairGenotyper.addFirstPass(record, discordantPairEvidence, depthResult, masterSampleList);
            }
        }
    }

    private void applyDiscordantPairSecondPass(final SVCallRecord record) {
        if (discordantPairCollectionEnabled()
                && depthGenotypeResults.containsKey(record.getId())
                && retainDiscordantPairTrainingRecord(record)) {
            final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = depthGenotypeResults.get(record.getId());
            discordantPairGenotyper.addSecondPass(record, depthResult, masterSampleList);
        }
    }

    private void applyDiscordantPairThirdPass(final SVCallRecord record) {
        if (discordantPairCollectionEnabled()
                && (!emitTrainingSubsetOnly || retainDiscordantPairTrainingRecord(record))) {
            final List<DiscordantPairEvidence> discordantPairEvidence = discordantPairCollector.collectEvidence(record);
            final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult genotypeResult = discordantPairGenotyper.genotype(record, discordantPairEvidence, discordantPairParameters, masterSampleList);
            if (discordantPairGenotypeResults.containsKey(record.getContigA())) {
                throw new UserException.BadInput("Duplicate variant ID: " + record.getId());
            }
            if (genotypeResult != null) {
                discordantPairGenotypeResults.put(record.getId(), genotypeResult);
            }
        }
    }

    private boolean discordantPairCollectionEnabled() {
        return discordantPairsFile != null;
    }

    private boolean splitReadCollectionEnabled() {
        return splitReadsFile != null;
    }

    private void applySplitReadFirstPass(final SVCallRecord record) {
        if (splitReadCollectionEnabled()
                && depthGenotypeResults.containsKey(record.getId())
                && discordantPairGenotypeResults.containsKey(record.getId())
                && retainSplitReadTrainingRecord(record)) {
            final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = depthGenotypeResults.get(record.getId());
            if (splitReadGenotyper.trainableRecord(record, discordantPairGenotyper, pesrExclusionEngine)) {
                final List<SplitReadEvidence> startSplitReads = splitReadStartCollector.collectEvidence(record);
                final List<SplitReadEvidence> endSplitReads = splitReadEndCollector.collectEvidence(record);
                splitReadGenotyper.addFirstPass(record, startSplitReads, endSplitReads, depthResult, masterSampleList);
            }
        }
    }

    private void applySplitReadSecondPass(final SVCallRecord record) {
        if (splitReadCollectionEnabled()
                && depthGenotypeResults.containsKey(record.getId())
                && retainSplitReadTrainingRecord(record)) {
            final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = depthGenotypeResults.get(record.getId());
            splitReadGenotyper.addSecondPass(record, depthResult, masterSampleList);
        }
    }

    private void applySplitReadThirdPass(final SVCallRecord record) {
        if (splitReadCollectionEnabled()
                && (!emitTrainingSubsetOnly || retainSplitReadTrainingRecord(record))) {
            final List<SplitReadEvidence> startSplitReads = splitReadStartCollector.collectEvidence(record);
            final List<SplitReadEvidence> endSplitReads = splitReadEndCollector.collectEvidence(record);
            final DepthEvidenceGenotyper.DepthGenotypeResult depthGenotype = depthGenotypeResults.get(record.getId());
            final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult discordantPairGenotype = discordantPairGenotypeResults.get(record.getId());
            final SplitReadEvidenceGenotyper.SplitReadGenotypeResult genotypeResult = splitReadGenotyper.genotypeTraining(record, startSplitReads, endSplitReads, depthGenotype, discordantPairGenotype, splitReadParameters, masterSampleList);
            if (splitReadGenotypeResults.containsKey(record.getContigA())) {
                throw new UserException.BadInput("Duplicate variant ID: " + record.getId());
            }
            if (genotypeResult != null) {
                splitReadGenotypeResults.put(record.getId(), genotypeResult);
            }
        }
    }

    private GATKPath getTablePath(final String suffix) {
        return new GATKPath(Path.of(tablesDir.toString(), tableBaseName + suffix).toString());
    }

    @Override
    public Object onTraversalSuccess() {
        if (writer != null) {
            writer.close();
        }
        return null;
    }

    private VCFHeader createHeader(VCFHeader header) {
        ploidyTable = new PloidyTable(ploidyTablePath.toPath());
        final List<String> featureSamples = ((SVFeaturesHeader) depthSource.getHeader()).getSampleNames();
        for (final String sample : featureSamples) {
            Utils.validate(ploidyTable.contains(sample), "Ploidy table does not contain sample " + sample + " from the depth file");
        }
        header = new VCFHeader(header.getMetaDataInInputOrder(), featureSamples);
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.DEPTH_COPY_STATE_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Depth genotyping copy state"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.DEPTH_MEDIAN_COPY_RATIO, 1, VCFHeaderLineType.Float, "Median read depth copy ratio"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.DEPTH_GENOTYPE_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Depth genotyping quality"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 1, VCFHeaderLineType.Integer, "Expected copy number for ref genotype"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.DISCORDANT_PAIR_GENOTYPE_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Discordant pair genotype"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.DISCORDANT_PAIR_GENOTYPE_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Discordant pair genotyping quality"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.SPLIT_READ_GENOTYPE_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read genotype"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.SPLIT_READ_GENOTYPE_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read genotyping quality"));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY));
        return header;
    }
}