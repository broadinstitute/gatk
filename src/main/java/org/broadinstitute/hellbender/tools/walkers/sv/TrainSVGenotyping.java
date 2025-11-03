package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
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
import java.util.*;
import java.util.function.Function;

/**
 * <p>This tool assesses depth evidence for structural variants (SVs),annotating records with statistical metrics that
 * can be used to assess a variant's quality. The input VCF should
 * contain multiple samples with GT fields populated. Note that this tool only considers carrier status and does not
 * differentiate heterozygous from homozygous variant genotypes. For paired-end, split-read, and B-allele frequency
 * evidence metrics, see the AggregateSVEvidence tool.</p>
 *
 * <p>Detailed methodology can be found in the supplement of <a href="https://doi.org/10.1038/s41586-020-2287-8">Collins et al. 2020</a>.</p>
 *
 * <p>Briefly, for each variant the median binned read counts are aggreagated across samples. Carrier and control
 * samples are then compared using a permuted t-test. However, if the test is underpowered, </p>
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         SV VCF
 *     </li>
 *     <li>
 *         Depth evidence file
 *     </li>
 *     <li>
 *         Median binned read counts table
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         SV VCF
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk AggregateDepthEvidence \
 *      -V input.vcf.gz \
 *      -O output.vcf.gz \
 *      --median-coverage median_coverage.tsv \
 *      --rd-file all_samples.rd.txt.gz
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        summary = "Read depth evidence assessment tool for copy number variants",
        oneLineSummary = "Read depth evidence assessment tool for copy number variants",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
public final class TrainSVGenotyping extends MultiplePassVariantWalker {

    public static final String DEPTH_EVIDENCE_FILE_PATH_LONG_NAME = "rd-file";
    public static final String MEDIAN_COUNTS_FILE_PATH_LONG_NAME = AggregateSVEvidence.MEDIAN_COVERAGE_LONG_NAME;
    public static final String NUM_TRAINING_STATES_LONG_NAME = "n-training-states";
    public static final String CUTOFFS_OUTPUT_LONG_NAME = "cutoffs-output";
    public static final String TRAINING_INTERVALS_LONG_NAME = "training-intervals";
    public static final String PESR_EXCLUSION_INTERVALS_LONG_NAME = "pesr-exclusion-intervals";
    public static final String MIN_PE_QUALITY_LONG_NAME = "pe-quality";
    public static final String MIN_PESER_SIZE_LONG_NAME = "min-pesr-size";
    public static final String MIN_SR_QUALITY_LONG_NAME = "sr-quality";

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
    public double minDiscordantPairQuality = 20;

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
    public double minSplitReadQuality = 20;

    @Argument(
            fullName = AggregateDepthEvidence.MAX_QUALITY_LONG_NAME,
            doc = "Max quality score",
            minValue = 1
    )
    public int maxQual = 999;

    @Argument(
            fullName = AggregateDepthEvidence.POWER_THRESHOLD_LONG_NAME,
            doc = "Power threshold for permuted T tests",
            minValue = 0,
            maxValue = 1
    )
    public double powerThreshold = 0.8;

    @Argument(
            fullName = CUTOFFS_OUTPUT_LONG_NAME,
            doc = "Enables genotype training and produces cutoffs table (.tsv) at the designated path"
    )
    public GATKPath cutoffsOutputPath;

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

    public static final String COPY_STATE_COLUMN = "copy_state";
    public static final String MEAN_COLUMN = "mean";
    public static final String STD_DEV_COLUMN = "sd";
    public static final String CUTOFFS_COLUMN = "cutoffs";
    private static final TableColumnCollection CUTOFFS_COLUMNS = new TableColumnCollection(Arrays.asList(COPY_STATE_COLUMN, MEAN_COLUMN, STD_DEV_COLUMN, CUTOFFS_COLUMN));

    private SAMSequenceDictionary dictionary;
    private VariantContextWriter writer;
    private Map<String, Double> sampleMedians;
    private PloidyTable ploidyTable;
    private DepthMatrixLoader loader;
    private DepthEvidenceGenotyper depthGenotyper;
    private List<String> masterSampleList;
    private FeatureDataSource<DepthEvidence> evidenceDataSource;
    private List<DepthEvidenceGenotyper.CopyStateStats> trainedCopyStateStats;

    private Map<String, DepthEvidenceGenotyper.DepthGenotypeResult> depthGenotypeResults;
    private Map<String, DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult> discordantPairGenotypeResults;

    private FeatureDataSource<DiscordantPairEvidence> discordantPairSource;
    private DiscordantPairEvidenceAggregator discordantPairCollector;
    private DiscordantPairEvidenceGenotyper discordantPairGenotyper;
    private DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeParameters discordantPairParameters;

    private FeatureDataSource<SplitReadEvidence> splitReadSource;
    private SplitReadEvidenceAggregator splitReadStartCollector;
    private SplitReadEvidenceAggregator splitReadEndCollector;
    private SplitReadEvidenceGenotyper splitReadGenotyper;

    private SVStratificationEngine pesrExclusionEngine;

    private static final int DISCORDANT_PAIR_QUERY_LOOKAHEAD = 0;
    private static final int SPLIT_READ_QUERY_LOOKAHEAD = 0;
    private static final String PESR_EXCLUSION_STRATIFICATION = "pesrex";

    protected int numberOfPasses() {
        return 6;
    }

    private void initializeDiscordantPairCollection() {
        if (discordantPairCollectionEnabled()) {
            discordantPairSource = new FeatureDataSource<>(
                    discordantPairsFile.toString(),
                    "discordantPairsFile",
                    DISCORDANT_PAIR_QUERY_LOOKAHEAD,
                    DiscordantPairEvidence.class,
                    cloudPrefetchBuffer,
                    cloudIndexPrefetchBuffer);
            discordantPairCollector = new DiscordantPairEvidenceAggregator(discordantPairSource, dictionary, innerWindow, outerWindow);
        }
    }

    private void initializeSplitReadCollection() {
        if (splitReadCollectionEnabled()) {
            splitReadSource = new FeatureDataSource<>(
                    splitReadsFile.toString(),
                    "splitReadsFile",
                    SPLIT_READ_QUERY_LOOKAHEAD,
                    SplitReadEvidence.class,
                    cloudPrefetchBuffer,
                    cloudIndexPrefetchBuffer);
            splitReadStartCollector = new SplitReadEvidenceAggregator(splitReadSource, dictionary, 0, true);
            splitReadEndCollector = new SplitReadEvidenceAggregator(splitReadSource, dictionary, 0, false);
        }
    }

    @Override
    public void onTraversalStart() {
        evidenceDataSource = new FeatureDataSource<>(depthEvidenceFile.toString());
        sampleMedians = loadMedianSampleCoverageTable();
        dictionary = getBestAvailableSequenceDictionary();
        loader = new DepthMatrixLoader(evidenceDataSource, numBins, largeVariantSize, largeVariantPoints, largeVariantWindow, dictionary);
        writer = createVCFWriter(outputVcf);
        final VCFHeader header = createHeader(getHeaderForVariants());
        writer.writeHeader(header);
        masterSampleList = header.getSampleNamesInOrder();
        depthGenotyper = new DepthEvidenceGenotyper(null, masterSampleList, dictionary);
        trainCopyNumberSites();

        depthGenotypeResults = new HashMap<>();
        discordantPairGenotypeResults = new HashMap<>();

        if (splitReadCollectionEnabled()) {
            Utils.validate(discordantPairCollectionEnabled(), "Discordant pairs file must be provided for split read training");
        }

        initializeDiscordantPairCollection();
        initializeSplitReadCollection();
        discordantPairGenotyper = new DiscordantPairEvidenceGenotyper(sampleMedians, minDiscordantPairQuality, minPesrSize, PESREvidenceTester.DEPTH_BASIS, maxQual);
        splitReadGenotyper = new SplitReadEvidenceGenotyper(sampleMedians, minSplitReadQuality, minPesrSize, PESREvidenceTester.DEPTH_BASIS, maxQual);

        if (pesrExclusionIntervalsPath != null) {
            pesrExclusionEngine = new SVStratificationEngine(dictionary);
            final GenomeLocParser genomeLocParser = new GenomeLocParser(dictionary);
            final GenomeLocSortedSet genomeLocs = IntervalUtils.loadIntervals(Collections.singletonList(pesrExclusionIntervalsPath.toString()), IntervalSetRule.UNION, IntervalMergingRule.OVERLAPPING_ONLY, 0, genomeLocParser);
            final List<Locatable> intervals = Collections.unmodifiableList(genomeLocs.toList());
            pesrExclusionEngine.addTrack(PESR_EXCLUSION_STRATIFICATION, intervals);
            pesrExclusionEngine.addStratification(PESR_EXCLUSION_STRATIFICATION, null, null, null, Collections.singleton(PESR_EXCLUSION_STRATIFICATION));
        }
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
        try (final TableWriter<DepthEvidenceGenotyper.CopyStateStats> tableWriter = TableUtils.writer(cutoffsOutputPath.toPath(), CUTOFFS_COLUMNS, this::composeCutoffsLine)) {
            tableWriter.writeAllRecords(trainedCopyStateStats);
        } catch (IOException e) {
            throw new GATKException("Error while writing cutoffs table", e);
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
            applySplitRead(record);
        } else if (n == 5) {
            writeGenotypes(record);
        } else {
            throw new GATKException("Unexpected number of passes: " + n);
        }
    }

    public void writeGenotypes(final SVCallRecord record) {
        final ArrayList<Genotype> newGenotypeList = new ArrayList<>(masterSampleList.size());
        final GenotypesContext genotypes = record.getGenotypes();
        for (int i = 0; i < masterSampleList.size(); i++) {
            final String sample = masterSampleList.get(i);
            if (!genotypes.containsSample(sample)) {
                throw new IllegalArgumentException("Sample " + sample + " does not exist in record " + record.getId());
            }
            final GenotypeBuilder builder = new GenotypeBuilder(genotypes.get(sample));
            if (depthGenotypeResults.containsKey(record.getId())) {
                final DepthEvidenceGenotyper.DepthGenotypeResult result = depthGenotypeResults.get(record.getId());
                final DepthEvidenceGenotyper.DepthGenotype depthGenotype =  result.genotypeQuals().get(i);
                builder.attribute(GATKSVVCFConstants.DEPTH_GENOTYPE_COPY_NUMBER_FORMAT, depthGenotype.copyState());
                builder.attribute(GATKSVVCFConstants.DEPTH_GENOTYPE_QUALITY_ATTRIBUTE, (int) depthGenotype.quality());
            }
            if (discordantPairGenotypeResults.containsKey(record.getId())) {
                final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult result = discordantPairGenotypeResults.get(record.getId());
                builder.attribute(GATKSVVCFConstants.DISCORDANT_PAIR_GENOTYPE_ATTRIBUTE, result.genotypes()[i]);
                builder.attribute(GATKSVVCFConstants.DISCORDANT_PAIR_GENOTYPE_QUALITY_ATTRIBUTE, result.genotypeQuals()[i]);
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
        } else if (n == 4 && splitReadCollectionEnabled()) {
            splitReadGenotyper.finalizeFirstPass();
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
                && depthGenotypeResults.containsKey(record.getId())) {
            final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = depthGenotypeResults.get(record.getId());
            if (discordantPairGenotyper.trainableRecord(record, depthResult, pesrExclusionEngine)) {
                final List<DiscordantPairEvidence> discordantPairEvidence = discordantPairCollector.collectEvidence(record);
                discordantPairGenotyper.addFirstPass(record, discordantPairEvidence, depthResult, masterSampleList);
            }
        }
    }

    private void applyDiscordantPairSecondPass(final SVCallRecord record) {
        if (discordantPairCollectionEnabled()
                && depthGenotypeResults.containsKey(record.getId())) {
            final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = depthGenotypeResults.get(record.getId());
            discordantPairGenotyper.addSecondPass(record, depthResult, masterSampleList);
        }
    }

    private void applyDiscordantPairThirdPass(final SVCallRecord record) {
        if (discordantPairCollectionEnabled()) {
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

    private void applySplitRead(final SVCallRecord record) {
        if (splitReadCollectionEnabled()
                && depthGenotypeResults.containsKey(record.getId())
                && discordantPairGenotypeResults.containsKey(record.getId())) {
            final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = depthGenotypeResults.get(record.getId());
            final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult discorantPairResult = discordantPairGenotypeResults.get(record.getId());
            if (splitReadGenotyper.trainableRecord(record, depthResult, discorantPairResult, pesrExclusionEngine)) {
                final List<SplitReadEvidence> startSplitReads = splitReadStartCollector.collectEvidence(record);
                final List<SplitReadEvidence> endSplitReads = splitReadEndCollector.collectEvidence(record);
                splitReadGenotyper.addFirstPass(record, startSplitReads, endSplitReads, depthResult, masterSampleList);
            }
        }
    }

    @Override
    public Object onTraversalSuccess() {
        if (writer != null) {
            writer.close();
        }
        return null;
    }

    private void composeCutoffsLine(DepthEvidenceGenotyper.CopyStateStats stats, DataLine dataLine) {
        dataLine.append(stats.copyState());
        dataLine.append(stats.mean());
        dataLine.append(stats.stdDev());
        dataLine.append(stats.upperBound());
    }

    private Function<DataLine, DepthEvidenceGenotyper.CopyStateStats> tableParser(TableColumnCollection columns, Function<String, RuntimeException> exceptionFactory) {
        // Check for expected columns
        for (final String column : CUTOFFS_COLUMNS.names()) {
            if (!columns.contains(column)) {
                throw exceptionFactory.apply("Missing column " + column);
            }
        }
        // Check there are no extra columns
        if (columns.columnCount() != CUTOFFS_COLUMNS.columnCount()) {
            throw exceptionFactory.apply("Expected " + columns.columnCount() + " columns but found " + columns.columnCount());
        }
        return this::parseTableLine;
    }

    private DepthEvidenceGenotyper.CopyStateStats parseTableLine(final DataLine dataLine) {
        final int copyState = Integer.parseInt(dataLine.get(COPY_STATE_COLUMN));
        final double mean = Double.parseDouble(dataLine.get(MEAN_COLUMN));
        final double stdDev = Double.parseDouble(dataLine.get(STD_DEV_COLUMN));
        final double cutoff = Double.parseDouble(dataLine.get(CUTOFFS_COLUMN));
        return new DepthEvidenceGenotyper.CopyStateStats(copyState, mean, stdDev, cutoff);
    }

    private VCFHeader createHeader(VCFHeader header) {
        final List<String> samples = header.getSampleNamesInOrder();
        if (samples.isEmpty()) {
            Utils.nonNull(ploidyTablePath, "Ploidy table required for sites-only VCFs");
            ploidyTable = new PloidyTable(ploidyTablePath.toPath());
            final List<String> featureSamples = ((SVFeaturesHeader) evidenceDataSource.getHeader()).getSampleNames();
            for (final String sample : samples) {
                Utils.validate(ploidyTable.contains(sample), "Ploidy table does not contain sample " + sample + " from the depth file");
            }
            header = new VCFHeader(header.getMetaDataInInputOrder(), featureSamples);
        }
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.DEPTH_VARIANT_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Float, "Depth genotyping variant quality"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.DEPTH_COPY_STATE_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Depth genotyping copy state"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.DEPTH_GENOTYPE_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Depth genotyping quality"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 1, VCFHeaderLineType.Integer, "Expected copy number for ref genotype"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.DISCORDANT_PAIR_GENOTYPE_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Discordant pair genotype"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.DISCORDANT_PAIR_GENOTYPE_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Discordant pair genotyping quality"));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY));
        return header;
    }
}