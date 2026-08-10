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
 *         Genotyping cutoff tables: separate RD tables for depth-only and PESR variants, plus PE and SR tables
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
    public static final String MIN_PESR_SIZE_LONG_NAME = "min-pesr-size";
    public static final String MIN_SR_QUALITY_LONG_NAME = "sr-quality";
    public static final String DISCORDANT_PAIR_QUERY_LOOKAHEAD_LONG_NAME = "pe-query-lookahead";
    public static final String SPLIT_READ_QUERY_LOOKAHEAD_LONG_NAME = "sr-query-lookahead";
    public static final String DEPTH_MIN_SEPARATION_LONG_NAME = "rd-depth-min-separation";
    public static final String PESR_MIN_SEPARATION_LONG_NAME = "rd-pesr-min-separation";
    public static final String OUTPUT_TRAINING_VCF_LONG_NAME = "output-training-vcf";

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
            doc = "Output VCF (required when --" + OUTPUT_TRAINING_VCF_LONG_NAME + " is true)",
            optional = true
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
            fullName = MIN_PESR_SIZE_LONG_NAME,
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
        public int discordantPairQueryLookahead = 0;

        @Argument(
            fullName = SPLIT_READ_QUERY_LOOKAHEAD_LONG_NAME,
            doc = "Number of bases to prefetch after SR evidence query cache misses",
            minValue = 0,
            optional = true
        )
        public int splitReadQueryLookahead = 0;

    @Argument(
            fullName = DEPTH_MIN_SEPARATION_LONG_NAME,
            doc = "Minimum RD median-ratio separation between copy states 1 and 2 for depth-only variants. " +
                    "Enforces that copy-state 1 upper bound <= 1 - sep and copy-state 2 upper bound >= 1 + sep. " +
                    "Higher values produce more conservative depth-only genotyping.",
            minValue = 0,
            optional = true
    )
    public double depthMinSeparation = 0.;

    @Argument(
            fullName = PESR_MIN_SEPARATION_LONG_NAME,
            doc = "Minimum RD median-ratio separation between copy states 1 and 2 for PESR variants. " +
                    "Enforces that copy-state 1 upper bound <= 1 - sep and copy-state 2 upper bound >= 1 + sep. " +
                    "Typically smaller than the depth-only separation since PE/SR evidence compensates.",
            minValue = 0,
            optional = true
    )
    public double pesrMinSeparation = 0.;

    @Argument(
            fullName = OUTPUT_TRAINING_VCF_LONG_NAME,
            doc = "Write genotyped training VCF output. When false (the default), only parameter " +
                    "tables are produced and the -O argument is not required. Enable for debugging " +
                    "or integration testing.",
            optional = true
    )
    public boolean outputTrainingVcf = false;

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
    private DepthEvidenceGenotyper depthGenotyper; // initial genotyper with default cutoffs, used for training
    private DepthEvidenceGenotyper depthOnlyGenotyper; // trained genotyper for depth-only variants
    private DepthEvidenceGenotyper pesrDepthGenotyper; // trained genotyper for PESR variants
    private List<String> masterSampleList;
    private FeatureDataSource<DepthEvidence> depthSource;

    // Training-phase only: depth results keyed by variant ID for the trainable subset
    private Map<String, DepthEvidenceGenotyper.DepthGenotypeResult> trainingDepthResults;

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

    // Traversal tallies, written to the SR cutoff diagnostics file. Observational only.
    private long diagPhase1DepthResults = 0;
    private long diagPhase1PeTrainable = 0;
    private long diagPhase1SrTrainable = 0;
    private long diagPhase2VariantsVisited = 0;
    private long diagPhase2VariantsProcessed = 0;

    private static final String PESR_EXCLUSION_STRATIFICATION = "pesrex";

    protected int numberOfPasses() {
        // The actual traversal is managed by the overridden traverse() below.
        // This value is retained for compatibility with the parent contract.
        return 2;
    }

    /**
     * Two-phase streaming architecture that minimizes memory consumption:
     *
     * Phase 1 (Training): Reads the VCF once, identifies trainable variants, and learns
     * PE/SR genotyping parameters. Only trainable variants (typically ~5% of total) and their
     * depth genotype results are cached in memory. After training completes, all training
     * data is freed.
     *
     * Phase 2 has two modes depending on {@code --output-training-vcf}:
     *
     * Phase 2a (Full Genotype + Write, when --output-training-vcf is true): Re-reads the VCF
     * from disk, streaming each variant through depth/PE/SR genotyping with the trained
     * parameters, and writes immediately.
     *
     * Phase 2b (SR Histogram Only, when --output-training-vcf is false): Re-reads the VCF
     * from disk but only queries SR evidence for recovery histogram accumulation, skipping
     * depth and PE evidence queries entirely.
     *
     * In both modes, SR recovery statistics are maintained as fixed-size histograms (352 bytes
     * total) rather than unbounded lists. No per-variant data accumulates across records.
     *
     * Memory usage scales with the number of trainable variants (~100-200K) rather than the
     * total VCF size (3-4M), reducing peak memory from ~50+ GB to ~2-4 GB.
     */
    @Override
    public void traverse() {
        final CountingVariantFilter countingVariantFilter = makeVariantFilter();
        final CountingReadFilter readFilter = makeReadFilter();

        // ========== PHASE 1: TRAINING ==========
        // Scan the VCF once: register PE overlaps, depth-genotype all CNVs, and cache
        // only the trainable subset for PE/SR training passes.

        logger.info("Phase 1: Training — scanning VCF for overlap registration and depth genotyping");
        trainingDepthResults = new HashMap<>();
        final Map<String, DepthEvidenceGenotyper.DepthGenotypeResult> allDepthResults = new HashMap<>();
        final List<SVCallRecord> trainablePERecords = new ArrayList<>();

        StreamSupport.stream(getSpliteratorForDrivingVariants(), false)
                .filter(countingVariantFilter)
                .forEach(variant -> {
                    progressMeter.update(new SimpleInterval(variant));
                    // Register ALL PE-eligible variants for overlap checking.
                    if (discordantPairCollectionEnabled()) {
                        discordantPairGenotyper.registerVariantForOverlapCheck(variant);
                    }
                    final SVCallRecord record = SVCallRecordUtils.create(variant, dictionary);
                    final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = computeDepthGenotype(record);
                    if (depthResult != null) {
                        allDepthResults.put(record.getId(), depthResult);
                    }
                });
        diagPhase1DepthResults = allDepthResults.size();
        logger.info("Phase 1: Depth genotyping complete (" + allDepthResults.size() + " CNV results)");

        // Build the overlap detector from registered intervals
        if (discordantPairCollectionEnabled()) {
            discordantPairGenotyper.aggregateOverlapCheckIntervals();
        }
        progressMeter.reset();

        // Identify trainable records: filter allDepthResults to only those passing trainableRecord()
        // and cache lightweight SVCallRecords for the training passes.
        // We need a second VCF scan because trainableRecord() requires the overlap detector
        // which is only available after aggregateOverlapCheckIntervals().
        logger.info("Phase 1: Identifying trainable variants");
        StreamSupport.stream(getSpliteratorForDrivingVariants(), false)
                .filter(countingVariantFilter)
                .forEach(variant -> {
                    progressMeter.update(new SimpleInterval(variant));
                    final String id = variant.getID();
                    if (!allDepthResults.containsKey(id)) {
                        return;
                    }
                    final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = allDepthResults.get(id);
                    final SVCallRecord record = SVCallRecordUtils.create(variant, dictionary);
                    if (discordantPairCollectionEnabled()
                            && discordantPairGenotyper.trainableRecord(record, depthResult, pesrExclusionEngine)) {
                        trainablePERecords.add(record);
                        trainingDepthResults.put(id, depthResult);
                    }
                });
        // Free the full depth results map — only trainingDepthResults is needed for training
        allDepthResults.clear();
        diagPhase1PeTrainable = trainablePERecords.size();
        logger.info("Phase 1: Found " + trainablePERecords.size() + " PE-trainable variants");
        progressMeter.reset();

        // PE training: first pass (evidence collection for trainable records)
        logger.info("Phase 1: PE evidence collection");
        for (final SVCallRecord record : trainablePERecords) {
            final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = trainingDepthResults.get(record.getId());
            final List<DiscordantPairEvidence> evidence = discordantPairCollector.collectEvidence(record);
            discordantPairGenotyper.addFirstPass(record, evidence, depthResult, masterSampleList);
            progressMeter.update(getProgressInterval(record));
        }
        if (discordantPairCollectionEnabled()) {
            discordantPairGenotyper.finalizeFirstPass();
        }
        progressMeter.reset();

        // PE training: second pass (parameter estimation)
        logger.info("Phase 1: PE parameter estimation");
        for (final SVCallRecord record : trainablePERecords) {
            final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = trainingDepthResults.get(record.getId());
            discordantPairGenotyper.addSecondPass(record, depthResult, masterSampleList);
            progressMeter.update(getProgressInterval(record));
        }
        if (discordantPairCollectionEnabled()) {
            discordantPairParameters = discordantPairGenotyper.finalizeSecondPass();
            writeDiscordantPairParameters(discordantPairParameters);
        }
        progressMeter.reset();

        // PE genotype trainable records + SR first pass (merged)
        // We need PE genotype results for trainable records to feed into SR first pass.
        logger.info("Phase 1: PE genotyping + SR evidence collection (merged)");
        final Map<String, DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult> trainingPEResults = new HashMap<>();
        for (final SVCallRecord record : trainablePERecords) {
            if (discordantPairCollectionEnabled()) {
                final List<DiscordantPairEvidence> peEvidence = discordantPairCollector.collectEvidence(record);
                final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult peResult =
                        discordantPairGenotyper.genotype(record, peEvidence, discordantPairParameters, masterSampleList);
                if (peResult != null) {
                    trainingPEResults.put(record.getId(), peResult);
                }
            }
            if (splitReadCollectionEnabled() && trainingPEResults.containsKey(record.getId())) {
                final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = trainingDepthResults.get(record.getId());
                final boolean peEligible = true; // already filtered to PE-trainable records
                if (splitReadGenotyper.trainableRecord(record, peEligible, pesrExclusionEngine)) {
                    diagPhase1SrTrainable++;
                    final List<SplitReadEvidence> startSplitReads = splitReadStartCollector.collectEvidence(record);
                    final List<SplitReadEvidence> endSplitReads = splitReadEndCollector.collectEvidence(record);
                    splitReadGenotyper.addFirstPass(record, startSplitReads, endSplitReads, depthResult, masterSampleList);
                }
            }
            progressMeter.update(getProgressInterval(record));
        }
        if (splitReadCollectionEnabled()) {
            splitReadGenotyper.finalizeFirstPass();
        }
        // PE training data is no longer needed
        if (discordantPairCollectionEnabled()) {
            discordantPairGenotyper.clearTrainingData();
        }
        trainingPEResults.clear();
        progressMeter.reset();

        // SR training: second pass (parameter estimation, trainable records only)
        logger.info("Phase 1: SR parameter estimation");
        for (final SVCallRecord record : trainablePERecords) {
            if (splitReadCollectionEnabled()) {
                final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = trainingDepthResults.get(record.getId());
                final boolean peEligible = true; // already filtered
                if (splitReadGenotyper.trainableRecord(record, peEligible, pesrExclusionEngine)) {
                    splitReadGenotyper.addSecondPass(record, depthResult, masterSampleList);
                }
            }
            progressMeter.update(getProgressInterval(record));
        }
        if (splitReadCollectionEnabled()) {
            splitReadParameters = splitReadGenotyper.finalizeSecondPass();
        }
        progressMeter.reset();

        // Free all training data — only trained parameters survive
        trainingDepthResults.clear();
        trainingDepthResults = null;
        trainablePERecords.clear();
        logger.info("Phase 1 complete — all training parameters learned");

        if (outputTrainingVcf) {
            // ========== PHASE 2a: FULL GENOTYPE + WRITE ==========
            // Re-read the VCF from disk and process each variant with constant memory.
            // Depth genotypes are re-computed from the tabix-indexed RD file.
            // PE and SR genotypes use the trained parameters. SR recovery stats accumulate in
            // fixed-size histograms (352 bytes). No per-variant data persists across records.

            logger.info("Phase 2: Streaming genotype + write");
            StreamSupport.stream(getSpliteratorForDrivingVariants(), false)
                    .filter(countingVariantFilter)
                    .forEach(variant -> {
                        progressMeter.update(new SimpleInterval(variant));
                        diagPhase2VariantsVisited++;
                        diagPhase2VariantsProcessed++;
                        final SVCallRecord record = SVCallRecordUtils.create(variant, dictionary);

                        // Depth genotype (re-computed from tabix)
                        final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = computeDepthGenotype(record);

                        // PE genotype
                        DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult peResult = null;
                        if (discordantPairCollectionEnabled()) {
                            final List<DiscordantPairEvidence> peEvidence = discordantPairCollector.collectEvidence(record);
                            peResult = discordantPairGenotyper.genotype(record, peEvidence, discordantPairParameters, masterSampleList);
                        }

                        // SR genotype + recovery histogram accumulation
                        SplitReadEvidenceGenotyper.SplitReadGenotypeResult srResult = null;
                        if (splitReadCollectionEnabled()) {
                            final List<SplitReadEvidence> startSplitReads = splitReadStartCollector.collectEvidence(record);
                            final List<SplitReadEvidence> endSplitReads = splitReadEndCollector.collectEvidence(record);
                            srResult = splitReadGenotyper.genotypeTraining(record, startSplitReads, endSplitReads, depthResult, peResult, splitReadParameters, masterSampleList);
                        }

                        // Write immediately — no accumulation
                        writeGenotypes(record, depthResult, peResult, srResult);
                    });
            logger.info("Phase 2: Genotype + write complete");
        } else {
            // ========== PHASE 2b: SR HISTOGRAM ACCUMULATION ONLY ==========
            // When not writing the training VCF, Phase 2 only needs SR evidence queries
            // for recovery histogram accumulation. Depth and PE evidence queries are skipped
            // because the histogram's pass flag depends only on SVType (CNV <-> depthGenotype
            // != null) and PE GQ which is always >= 1. See accumulateHistogramOnly() for details.

            logger.info("Phase 2: SR histogram accumulation");
            StreamSupport.stream(getSpliteratorForDrivingVariants(), false)
                    .filter(countingVariantFilter)
                    .forEach(variant -> {
                        progressMeter.update(new SimpleInterval(variant));
                        diagPhase2VariantsVisited++;
                        if (!splitReadCollectionEnabled()) {
                            return;
                        }
                        diagPhase2VariantsProcessed++;
                        final SVCallRecord record = SVCallRecordUtils.create(variant, dictionary);
                        final GATKSVVCFConstants.StructuralVariantAnnotationType svtype = record.getType();
                        final boolean isCNV = svtype == GATKSVVCFConstants.StructuralVariantAnnotationType.DEL
                                || svtype == GATKSVVCFConstants.StructuralVariantAnnotationType.DUP
                                || svtype == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV;
                        final List<SplitReadEvidence> startSplitReads = splitReadStartCollector.collectEvidence(record);
                        final List<SplitReadEvidence> endSplitReads = splitReadEndCollector.collectEvidence(record);
                        splitReadGenotyper.accumulateHistogramOnly(record, startSplitReads, endSplitReads, isCNV, splitReadParameters, masterSampleList);
                    });
            logger.info("Phase 2: SR histogram accumulation complete");
        }

        // Finalize SR recovery cutoffs from histograms and write params
        if (splitReadCollectionEnabled()) {
            try {
                splitReadFrequencyCutoffs = splitReadGenotyper.finalizeThirdPass();
                writeSplitReadParameters(new SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics(splitReadParameters, splitReadFrequencyCutoffs));
            } finally {
                // Always emit diagnostics: a rejected grid is exactly when they are needed, and
                // this tool must still exit zero for Cromwell to delocalize them.
                writeSplitReadCutoffDiagnostics();
            }
            logSplitReadCutoffSelectionOutcomes();
        }

        logger.info(countingVariantFilter.getSummaryLine());
        logger.info(readFilter.getSummaryLine());
    }

    private SimpleInterval getProgressInterval(final SVCallRecord record) {
        return new SimpleInterval(record.getContigA(), record.getPositionA(), record.getPositionA());
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
        if (outputTrainingVcf) {
            Utils.validate(outputVcf != null, "Output VCF path (-O) is required when --" + OUTPUT_TRAINING_VCF_LONG_NAME + " is true");
            writer = createVCFWriter(outputVcf);
        }
        outputHeader = createHeader(getHeaderForVariants());
        if (writer != null) {
            writer.writeHeader(outputHeader);
        }
        masterSampleList = outputHeader.getSampleNamesInOrder();
        depthGenotyper = new DepthEvidenceGenotyper(null, masterSampleList, maxQual, dictionary);
        trainCopyNumberSites();

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
        final List<DepthEvidenceGenotyper.CopyStateStats> baseCopyStateStats = depthGenotyper.train(genotypeResults, numTrainingStates);

        // Apply minimum separation constraints to produce separate cutoffs for depth-only and PESR variants.
        // This mirrors the old pipeline (v1.1 TrainRDGenotyping.wdl / UpdateCutoff) which applied different
        // RD_Median_Separation values from RF cutoffs to the trained boundaries.
        final List<DepthEvidenceGenotyper.CopyStateStats> depthOnlyStats = applySeparation(baseCopyStateStats, depthMinSeparation);
        final List<DepthEvidenceGenotyper.CopyStateStats> pesrStats = applySeparation(baseCopyStateStats, pesrMinSeparation);

        writeDepthCutoffs(depthOnlyStats, ".rd_depth_geno_params.tsv");
        writeDepthCutoffs(pesrStats, ".rd_pesr_geno_params.tsv");

        depthOnlyGenotyper = new DepthEvidenceGenotyper(depthOnlyStats, masterSampleList, maxQual, dictionary);
        pesrDepthGenotyper = new DepthEvidenceGenotyper(pesrStats, masterSampleList, maxQual, dictionary);
        logger.info("Training completed (depth-only separation=" + depthMinSeparation + ", PESR separation=" + pesrMinSeparation + ")");
    }

    /**
     * Applies minimum-separation constraints to trained copy-state boundaries.
     * For copy state 1 (single-copy deletion), the upper bound is capped at {@code 1 - minSeparation}.
     * For copy state 2 (diploid/reference), the upper bound is raised to at least {@code 1 + minSeparation}.
     * This enforces a minimum gap around the diploid ratio of 1.0, making genotyping more or less
     * conservative depending on the separation value.
     */
    @VisibleForTesting
    static List<DepthEvidenceGenotyper.CopyStateStats> applySeparation(
            final List<DepthEvidenceGenotyper.CopyStateStats> stats, final double minSeparation) {
        if (minSeparation == 0) {
            return stats;
        }
        final List<DepthEvidenceGenotyper.CopyStateStats> adjusted = new ArrayList<>(stats.size());
        for (final DepthEvidenceGenotyper.CopyStateStats s : stats) {
            double upperBound = s.upperBound();
            if (s.copyState() == 1 && upperBound > 1.0 - minSeparation) {
                upperBound = 1.0 - minSeparation;
            } else if (s.copyState() == 2 && upperBound < 1.0 + minSeparation) {
                upperBound = 1.0 + minSeparation;
            }
            adjusted.add(new DepthEvidenceGenotyper.CopyStateStats(s.copyState(), s.mean(), s.stdDev(), upperBound));
        }
        return adjusted;
    }

    private void writeDepthCutoffs(final List<DepthEvidenceGenotyper.CopyStateStats> stats, final String suffix) {
        try (final TableWriter<DepthEvidenceGenotyper.CopyStateStats> tableWriter = TableUtils.writer(
                getTablePath(suffix).toPath(),
                DepthEvidenceGenotyper.DepthTableParser.CUTOFFS_COLUMNS,
                new DepthEvidenceGenotyper.DepthTableParser()::composeCutoffsLine)) {
            tableWriter.writeAllRecords(stats);
        } catch (IOException e) {
            throw new GATKException("Error while writing RD cutoffs table " + suffix, e);
        }
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

    /**
     * Required by the {@link MultiplePassVariantWalker} contract but never invoked because
     * {@link #traverse()} is overridden. Retained for interface compliance.
     */
    @Override
    public void nthPassApply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext, final int n) {
        throw new GATKException("nthPassApply should not be called — traverse() is overridden");
    }

    /**
     * Computes depth genotype for a CNV record. Returns null for non-CNV types.
     */
    private DepthEvidenceGenotyper.DepthGenotypeResult computeDepthGenotype(final SVCallRecord record) {
        final GATKSVVCFConstants.StructuralVariantAnnotationType svtype = record.getType();
        if (svtype != GATKSVVCFConstants.StructuralVariantAnnotationType.DEL
                && svtype != GATKSVVCFConstants.StructuralVariantAnnotationType.DUP
                && svtype != GATKSVVCFConstants.StructuralVariantAnnotationType.CNV) {
            return null;
        }
        final DepthMatrix depthMatrix = loader.load(
                new SimpleInterval(record.getContigA(), record.getPositionA(), record.getPositionB()), sampleMedians);
        final DepthEvidenceGenotyper genotyper = record.isDepthOnly() ? depthOnlyGenotyper : pesrDepthGenotyper;
        return genotyper.genotype(depthMatrix);
    }

    /**
     * Writes genotype results for a single record. Accepts pre-computed results directly
     * rather than looking them up from maps, enabling streaming without accumulation.
     */
    public void writeGenotypes(final SVCallRecord record,
                               final DepthEvidenceGenotyper.DepthGenotypeResult depthResult,
                               final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult peResult,
                               final SplitReadEvidenceGenotyper.SplitReadGenotypeResult srResult) {
        final ArrayList<Genotype> newGenotypeList = new ArrayList<>(masterSampleList.size());
        // TODO: need to decide on whether to support genotyped VCF input and make consistent with GenotypeSVs
        final GenotypesContext genotypes = SVCallRecordUtils.populateGenotypesForMissingSamplesWithAlleles(
                record, new HashSet<>(masterSampleList), false, ploidyTable, outputHeader);
        for (int i = 0; i < masterSampleList.size(); i++) {
            final String sample = masterSampleList.get(i);
            if (!genotypes.containsSample(sample)) {
                throw new IllegalArgumentException("Sample " + sample + " does not exist in record " + record.getId());
            }
            final GenotypeBuilder builder = new GenotypeBuilder(genotypes.get(sample));
            if (depthResult != null) {
                builder.attribute(GATKSVVCFConstants.DEPTH_GENOTYPE_COPY_NUMBER_FORMAT, depthResult.copyStates()[i]);
                builder.attribute(GATKSVVCFConstants.DEPTH_MEDIAN_COPY_RATIO, depthResult.sampleDepths()[i]);
                builder.attribute(GATKSVVCFConstants.DEPTH_GENOTYPE_QUALITY_ATTRIBUTE, (int) depthResult.genotypeQuals()[i]);
            }
            if (peResult != null) {
                builder.attribute(GATKSVVCFConstants.DISCORDANT_PAIR_GENOTYPE_ATTRIBUTE, peResult.genotypes()[i]);
                builder.attribute(GATKSVVCFConstants.DISCORDANT_PAIR_GENOTYPE_QUALITY_ATTRIBUTE, peResult.genotypeQuals()[i]);
            }
            if (srResult != null) {
                builder.attribute(GATKSVVCFConstants.SPLIT_READ_GENOTYPE_ATTRIBUTE, srResult.genotypes()[i]);
                if (srResult.genotypeQuals() != null) {
                    builder.attribute(GATKSVVCFConstants.SPLIT_READ_GENOTYPE_QUALITY_ATTRIBUTE, srResult.genotypeQuals()[i]);
                }
            }
            newGenotypeList.add(builder.make());
        }
        final GenotypesContext newGenotypes = GenotypesContext.create(newGenotypeList);
        final SVCallRecord regenotypedCall = SVCallRecordUtils.copyCallWithNewGenotypes(record, newGenotypes);
        final VariantContext variant = SVCallRecordUtils.getVariantBuilder(regenotypedCall).make();
        writer.add(variant);
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

    /**
     * Write the SR frequency-cutoff diagnostic report to
     * {@code <output-dir>/<output-name>.sr_cutoff_diagnostics.txt}.
     *
     * <p>The cutoffs in {@code .sr_geno_params.tsv} record only the winning grid cell, which
     * is not enough to tell a real optimum from a degenerate grid: if the frac histograms
     * concentrate in one bin, every cell scores alike and the argmax returns the first cell,
     * (0.0, 0.0). This report carries the histograms, the fully scored grid, and the tie/NaN
     * counts needed to distinguish those cases. It goes to a file rather than the log so it
     * survives as a workflow output.</p>
     *
     * <p>Diagnostics are best-effort: a failure here must not fail a training run that has
     * already produced its parameter tables.</p>
     */
    private void writeSplitReadCutoffDiagnostics() {
        final GATKPath path = getTablePath(".sr_cutoff_diagnostics.txt");
        try (final java.io.Writer out = java.nio.file.Files.newBufferedWriter(path.toPath(), Charset.defaultCharset())) {
            out.write("## RUN_CONFIGURATION\n");
            out.write("tool\tTrainSVGenotyping\n");
            out.write("output_name\t" + tableBaseName + '\n');
            out.write("num_samples\t" + masterSampleList.size() + '\n');
            out.write("sr_quality_cutoff\t" + minSplitReadQuality + '\n');
            out.write("pe_quality_cutoff\t" + minDiscordantPairQuality + '\n');
            out.write("min_pesr_size\t" + minPesrSize + '\n');
            out.write("output_training_vcf\t" + outputTrainingVcf + '\n');
            out.write("phase2_mode\t" + (outputTrainingVcf ? "2a_full_genotype" : "2b_histogram_only") + '\n');
            out.write("max_quality\t" + maxQual + '\n');
            out.write("num_training_states\t" + numTrainingStates + '\n');
            out.write("depth_min_separation\t" + depthMinSeparation + '\n');
            out.write("pesr_min_separation\t" + pesrMinSeparation + '\n');

            out.write("## TRAVERSAL_TALLIES\n");
            out.write("phase1_depth_genotyped_cnvs\t" + diagPhase1DepthResults + '\n');
            out.write("phase1_pe_trainable_variants\t" + diagPhase1PeTrainable + '\n');
            out.write("phase1_sr_trainable_variants\t" + diagPhase1SrTrainable + '\n');
            out.write("phase2_variants_visited\t" + diagPhase2VariantsVisited + '\n');
            out.write("phase2_variants_processed\t" + diagPhase2VariantsProcessed + '\n');

            out.write(splitReadGenotyper.cutoffDiagnosticsReport());
            logger.info("Wrote SR cutoff diagnostics to " + path);
        } catch (final IOException | RuntimeException e) {
            logger.warn("Could not write SR cutoff diagnostics to " + path + ": " + e.getMessage());
        }
    }

    /**
     * Log any SR cutoff selection rejection prominently.
     *
     * <p>This tool does not fail on a rejection, because Cromwell delocalizes task outputs only
     * on success and failing here would strand the diagnostics report. Enforcement is the
     * ValidateSRCutoffs task, which reads the delocalized report.</p>
     */
    private void logSplitReadCutoffSelectionOutcomes() {
        for (final SplitReadEvidenceGenotyper.SelectionOutcome outcome : splitReadGenotyper.cutoffSelectionOutcomes()) {
            if (outcome.rejected()) {
                logger.warn("*** SR frequency cutoff selection REJECTED (" + outcome.status() + "). "
                        + "Cutoffs fell back to 0.0, which disables SR background filtering. "
                        + outcome.detail());
            }
        }
    }

    private boolean discordantPairCollectionEnabled() {
        return discordantPairsFile != null;
    }

    private boolean splitReadCollectionEnabled() {
        return splitReadsFile != null;
    }

    @Override
    protected void afterNthPass(final int n) {
        // All pass logic is handled inline in traverse(); this is a no-op.
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