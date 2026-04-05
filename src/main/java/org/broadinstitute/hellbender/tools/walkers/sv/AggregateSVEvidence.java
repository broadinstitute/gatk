package org.broadinstitute.hellbender.tools.walkers.sv;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.io.IOUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.*;
import org.broadinstitute.hellbender.tools.sv.aggregation.*;
import org.broadinstitute.hellbender.tools.sv.cluster.PloidyTable;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.nio.charset.Charset;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * <p>This tool assesses split read (SR), discordant paired-end (PE), and B-allele frequency (BAF) evidence for structural variants (SVs),
 * annotating records with statistical metrics that can be used to assess a variant's quality. The input VCF should
 * contain multiple samples with GT fields populated. Note that this tool only considers carrier status and does not
 * differentiate heterozygous from homozygous variant genotypes. For read depth evidence metrics, see
 * the AggregateDepthEvidence tool.</p>
 *
 * <p>Detailed methodology can be found in the supplement of <a href="https://doi.org/10.1038/s41586-020-2287-8">Collins et al. 2020</a>.</p>
 *
 * <p>Briefly, for each variant the supporting split reads and discordant pairs are counted. Phred-scaled quality scores
 * (SRQ, PEQ, PESRQ) are then computed based on a Poisson test of the observed median carrier sample signal against
 * background. The raw fraction of median SR signal attributed to carriers (SRCS, PECS, PESRCS) is also annotated as an additional
 * metric to assess concordance between detected evidence and genotypes.</p>
 *
 * <p>During SR aggregation, breakpoint refinement is performed (SR1POS, SR2POS) by maximizing the quality score
 * over all positions within a small window around each end of the variant.</p>
 *
 * <p>Bi-allelic copy number variants are also assessed using BAF evidence. Deletions are annotated with the ratio of heterozygous
 * SNPs in carrier samples to in controls (BAF_HET_RATIO). Duplications are assessed by comparing the distribution of
 * BAFs across SNPs with a Kolmogorov-Smirnov test statistic (BAF_KS_STAT), which is used to compute a quality
 * score (BAF_KS_Q).</p>
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         SV VCF
 *     </li>
 *     <li>
 *         PE evidence file (optional)
 *     </li>
 *     <li>
 *         SR evidence file (optional)
 *     </li>
 *     <li>
 *         BAF evidence file (optional)
 *     </li>
 *     <li>
 *         Median binned read counts table
 *     </li>
 *     <li>
 *         Ploidy table
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
 *     gatk AggregateSVEvidence \
 *      -V input.vcf.gz \
 *      -O output.vcf.gz \
 *      --median-coverage median_coverage.tsv \
 *      --ploidy-table ploidy_table.tsv \
 *      --baf-file all_samples.baf.txt.gz \
 *      --split-reads-file all_samples.sr.txt.gz \
 *      --discordant-pairs-file all_samples.pe.txt.gz
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */

@CommandLineProgramProperties(
        summary = "Annotate SVs with supporting evidence metrics",
        oneLineSummary = "Annotate SVs with supporting evidence metrics",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@BetaFeature
@DocumentedFeature
public final class AggregateSVEvidence extends VariantWalker {
    public static final String SPLIT_READ_LONG_NAME = "split-reads-file";
    public static final String DISCORDANT_PAIRS_LONG_NAME = "discordant-pairs-file";
    public static final String BAF_LONG_NAME = "baf-file";
    public static final String MEDIAN_COVERAGE_LONG_NAME = "median-coverage";
    public static final String PE_INNER_WINDOW_LONG_NAME = "pe-inner-window";
    public static final String PE_OUTER_WINDOW_LONG_NAME = "pe-outer-window";
    public static final String SR_WINDOW_LONG_NAME = "sr-window";
    public static final String SR_CROSSOVER_LONG_NAME = "sr-insertion-crossover";
    public static final String BAF_MIN_SIZE_LONG_NAME = "baf-min-size";
    public static final String BAF_MAX_SIZE_LONG_NAME = "baf-max-size";
    public static final String BAF_KS_SEED_LONG_NAME = "baf-ks-seed";
    public static final String BAF_GMM_SEED_LONG_NAME = "baf-gmm-seed";
    public static final String BAF_PADDING_FRACTION_LONG_NAME = "baf-padding-fraction";
    public static final String MIN_BAF_COUNT_LONG_NAME = "min-baf-count";
    public static final String X_CHROMOSOME_LONG_NAME = "x-chromosome-name";
    public static final String Y_CHROMOSOME_LONG_NAME = "y-chromosome-name";

    @Argument(
            doc = "Split reads evidence file (indexed and ending in .sr.txt.gz)",
            fullName = SPLIT_READ_LONG_NAME,
            optional = true
    )
    private GATKPath splitReadsFile;

    @Argument(
            doc = "Discordant pairs evidence file (indexed and ending in .pe.txt.gz)",
            fullName = DISCORDANT_PAIRS_LONG_NAME,
            optional = true
    )
    private GATKPath discordantPairsFile;

    @Argument(
            doc = "B-allele frequency (BAF) evidence file (indexed and ending in .baf.txt.gz)",
            fullName = BAF_LONG_NAME,
            optional = true
    )
    private GATKPath bafFile;

    @Argument(
            doc = "Tab-delimited table with sample IDs in the first row and median binned read counts in the second row.",
            fullName = MEDIAN_COVERAGE_LONG_NAME
    )
    private GATKPath medianCoverageFile;

    @Argument(
            doc = "Output VCF",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private String outputFile;

    /**
     * Expected format is tab-delimited and contains a header with the first column SAMPLE and remaining columns
     * contig names. Each row corresponds to a sample, with the sample ID in the first column and contig ploidy
     * integers in their respective columns. This information is used to determine sample sex.
     */
    @Argument(
            doc = "Sample ploidy table (.tsv)",
            fullName = SVCluster.PLOIDY_TABLE_LONG_NAME
    )
    private GATKPath ploidyTablePath;

    /**
     * Paired-end window size downstream of the first mate and upstream of the second mate.
     */
    @Argument(
            doc = "Inner discordant pair window size (bp)",
            fullName = PE_INNER_WINDOW_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private int innerWindow = 50;

    /**
     * Paired-end window size upstream of the first mate and downstream of the second mate.
     */
    @Argument(
            doc = "Outer discordant pair window size (bp)",
            fullName = PE_OUTER_WINDOW_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private int outerWindow = 500;

    @Argument(
            doc = "Split read window size (bp)",
            fullName = SR_WINDOW_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private int splitReadWindow = 50;

    /**
     * The split read signature for INS/DEL identified by searching for right-clipped (+ stranded) reads
     * upstream of the breakpoint and left-clipped (- stranded) reads downstream. In some instances, the left- and
     * right-clipped reads "cross over" such that the right-clipped position is downstream. This parameter adjusts the
     * maximum allowed distance that left- and right-clipped positions may cross over.
     */
    @Argument(
            doc = "Max split read crossover distance (bp) for INS/DEL",
            fullName = SR_CROSSOVER_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private int splitReadCrossover = SplitReadEvidenceTester.DEFAULT_MAX_CROSS_DISTANCE;

    @Argument(
            doc = "Minimum variant size for BAF aggregation",
            fullName = BAF_MIN_SIZE_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private int bafMinSize = 5000;

    @Argument(
            doc = "Maximum variant size for BAF aggregation",
            fullName = BAF_MAX_SIZE_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private int bafMaxSize = 10000000;

    @Argument(
            doc = "Seed for BAF Kolmogorov-Smirnov test PRNG",
            fullName = BAF_KS_SEED_LONG_NAME,
            optional = true
    )
    private long bafKsSeed = 93083432L;

    @Argument(
            doc = "Seed for BAF deletion GMM initialization",
            fullName = BAF_GMM_SEED_LONG_NAME,
            optional = true
    )
    private int bafGmmSeed = 0;

    @Argument(
            doc = "BAF flanking region size as a fraction of variant size",
            fullName = BAF_PADDING_FRACTION_LONG_NAME,
            minValue = 0.01,
            optional = true
    )
    private double bafPaddingFraction = 1.0;

    @Argument(
            doc = "Minimum number of BAF values required in carrier and non-carrier groups, used for duplication BAF assessment.",
            fullName = MIN_BAF_COUNT_LONG_NAME,
            minValue = 1,
            optional = true
    )
    private int minBafCount = 2;

    @Argument(
            doc = "X chromosome name",
            fullName = X_CHROMOSOME_LONG_NAME,
            optional = true
    )
    private String xChromosomeName = "chrX";

    @Argument(
            doc = "Y chromosome name",
            fullName = Y_CHROMOSOME_LONG_NAME,
            optional = true
    )
    private String yChromosomeName = "chrY";

    private SAMSequenceDictionary dictionary;
    private VariantContextWriter writer;

    // SR
    private FeatureDataSource<SplitReadEvidence> splitReadSource;
    private SplitReadEvidenceAggregator startSplitCollector;
    private SplitReadEvidenceAggregator endSplitCollector;
    private SplitReadEvidenceTester splitReadEvidenceTester;

    // PE
    private FeatureDataSource<DiscordantPairEvidence> discordantPairSource;
    private DiscordantPairEvidenceAggregator discordantPairCollector;
    private DiscordantPairEvidenceTester discordantPairEvidenceTester;

    // PESR
    private PESREvidenceTester pesrEvidenceTester;

    // BAF
    private FeatureDataSource<BafEvidence> bafSource;
    private BafEvidenceAggregator bafCollector;
    private BafHetRatioTester bafHetRatioTester;
    private BafKolmogorovSmirnovTester bafKolmogorovSmirnovTester;

    // Metadata
    private Map<String,Double> sampleCoverageMap;
    private Set<String> samples;
    private VCFHeader header;
    private Set<String> maleSamples;
    private Set<String> femaleSamples;

    private static final int BAF_QUERY_LOOKAHEAD = 0;
    private static final int SPLIT_READ_QUERY_LOOKAHEAD = 0;
    private static final int DISCORDANT_PAIR_QUERY_LOOKAHEAD = 0;

    @Override
    public void onTraversalStart() {
        progressMeter.setRecordsBetweenTimeChecks(100);
        dictionary = getBestAvailableSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }
        samples = new LinkedHashSet<>(getHeaderForVariants().getSampleNamesInOrder());
        if (!splitReadCollectionEnabled() && !discordantPairCollectionEnabled() && !bafCollectionEnabled()) {
            throw new UserException.BadInput("At least one evidence file must be provided");
        }
        sampleCoverageMap = loadSampleCoverage(medianCoverageFile);
        if (discordantPairCollectionEnabled()) {
            initializeDiscordantPairCollection();
        }
        if (splitReadCollectionEnabled()) {
            initializeSplitReadCollection();
            if (discordantPairCollectionEnabled()) {
                pesrEvidenceTester = new PESREvidenceTester(sampleCoverageMap);
            }
        }
        if (bafCollectionEnabled()) {
            initializeBAFCollection();
        }
        writer = createVCFWriter(Paths.get(outputFile));
        header = getVCFHeader();
        writer.writeHeader(header);
        initializeSampleSexSets();
    }

    /**
     * Read ploidy table and assign sample sex
     */
    private void initializeSampleSexSets() {
        final PloidyTable table = new PloidyTable(ploidyTablePath.toPath());
        maleSamples = header.getGenotypeSamples().stream()
                .filter(s -> table.get(s, xChromosomeName) == 1 && table.get(s, yChromosomeName) == 1)
                .collect(Collectors.toSet());
        femaleSamples = header.getGenotypeSamples().stream()
                .filter(s -> table.get(s, xChromosomeName) == 2 && table.get(s, yChromosomeName) == 0)
                .collect(Collectors.toSet());
    }

    private void initializeDiscordantPairCollection() {
        initializeDiscordantPairDataSource();
        discordantPairCollector = new DiscordantPairEvidenceAggregator(discordantPairSource, dictionary, innerWindow, outerWindow);
        discordantPairEvidenceTester = new DiscordantPairEvidenceTester(sampleCoverageMap, dictionary);
    }

    private void initializeSplitReadCollection() {
        initializeSplitReadEvidenceDataSource();
        splitReadEvidenceTester = new SplitReadEvidenceTester(sampleCoverageMap, splitReadCrossover, dictionary);
        startSplitCollector = new SplitReadEvidenceAggregator(splitReadSource, dictionary, splitReadWindow, true);
        endSplitCollector = new SplitReadEvidenceAggregator(splitReadSource, dictionary, splitReadWindow, false);
    }

    private void initializeBAFCollection() {
        initializeBAFEvidenceDataSource();
        bafCollector = new BafEvidenceAggregator(bafSource, dictionary, bafPaddingFraction);
        bafHetRatioTester = new BafHetRatioTester(bafGmmSeed);
        bafKolmogorovSmirnovTester = new BafKolmogorovSmirnovTester(minBafCount, bafKsSeed);
    }

    private void initializeDiscordantPairDataSource() {
        discordantPairSource = new FeatureDataSource<>(
                discordantPairsFile.toString(),
                "discordantPairsFile",
                DISCORDANT_PAIR_QUERY_LOOKAHEAD,
                DiscordantPairEvidence.class,
                cloudPrefetchBuffer,
                cloudIndexPrefetchBuffer);
    }

    private void initializeSplitReadEvidenceDataSource() {
        splitReadSource = new FeatureDataSource<>(
                splitReadsFile.toString(),
                "splitReadsFile",
                SPLIT_READ_QUERY_LOOKAHEAD,
                SplitReadEvidence.class,
                0,  // end positions can skip around and cause major prefetching slowdowns
                cloudIndexPrefetchBuffer);
    }

    private void initializeBAFEvidenceDataSource() {
        bafSource = new FeatureDataSource<>(
                bafFile.toString(),
                "bafFile",
                BAF_QUERY_LOOKAHEAD,
                BafEvidence.class,
                cloudPrefetchBuffer,
                cloudIndexPrefetchBuffer);
    }

    /**
     * Loads median coverage table
     */
    public static Map<String, Double> loadSampleCoverage(final GATKPath path) {
        final String fileString = path.toString();
        final List<String> lines = IOUtils.readLines(BucketUtils.openFile(fileString), Charset.defaultCharset());
        Utils.validate(lines.size() >= 2, "Median coverage file must contain at least two lines");
        final String[] samples = lines.get(0).split("\t");
        final String[] values = lines.get(1).split("\t");
        Utils.validate(samples.length == values.length,
                "Median file's first two lines must have the same number of columns");
        final Map<String, Double> sampleCoverageMap = new HashMap<>();
        try {
            for (int i = 0; i < samples.length; i++) {
                sampleCoverageMap.put(samples[i], Double.valueOf(values[i]));
            }
        } catch (final NumberFormatException nfe) {
            throw new UserException.BadInput("Encountered NumberFormatException while parsing median coverage file");
        }
        return sampleCoverageMap;
    }

    @Override
    public void closeTool() {
        if (discordantPairSource != null) {
            discordantPairSource.close();
        }
        if (splitReadSource != null) {
            splitReadSource.close();
        }
        if (bafSource != null) {
            bafSource.close();
        }
        if (writer != null) {
            writer.close();
        }
        super.closeTool();
    }

    /**
     * For skipping records that are not testable types
     */
    private boolean validRecordType(final SVCallRecord call) {
        return call.getType() != GATKSVVCFConstants.StructuralVariantAnnotationType.CNV;
    }

    /**
     * Determines which variants to run PE testing on
     */
    private boolean useDiscordantPairEvidence(final SVCallRecord call) {
        return !call.isDepthOnly() && call.getType() != GATKSVVCFConstants.StructuralVariantAnnotationType.INS;
    }

    /**
     * Determines which variants to run SR testing on
     */
    private boolean useSplitReadEvidence(final SVCallRecord call) {
        return !call.isDepthOnly();
    }

    /**
     * Determines which variants to run BAF testing on
     */
    private boolean useBafEvidence(final SVCallRecord call) {
        final Integer length = call.getLength();
        return (call.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.DEL || call.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.DUP)
                && length != null && length >= bafMinSize && length <= bafMaxSize;
    }

    /**
     * Returns set of samples to exclude for evidence stat calculations by sex. Prefer female carriers on chrX when
     * available, and male carriers on chrY.
     */
    public Set<String> getSamplesToExcludeForStatsBySex(final SVCallRecord record) {
        // TODO paired BND records may cause problems here
        if (record.getContigA().equals(xChromosomeName)) {
            Utils.validate(femaleSamples != null, "Ploidy table must be provided to call on X chromosome");
            final Set<String> carrierSamples = record.getCarrierSampleSet();
            final Collection<String> femaleCarriers = Sets.intersection(carrierSamples, femaleSamples);
            if (femaleCarriers.isEmpty()) {
                // Exclude no samples for X chromosome when there are no female carriers
                return Collections.emptySet();
            } else {
                // If there are female carriers on X, exclude males
                return maleSamples;
            }
        } else if (record.getContigA().equals(yChromosomeName)) {
            Utils.validate(maleSamples != null, "Ploidy table must be provided to call on Y chromosome");
            // Always exclude females for Y chromosome
            return femaleSamples;
        } else {
            // Default autosome
            return Collections.emptySet();
        }
    }

    /**
     * Perform evidence testing and write output
     */
    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext,
                                final ReferenceContext referenceContext, final FeatureContext featureContext) {
        SVCallRecord record = SVCallRecordUtils.create(variant, dictionary);
        if (validRecordType(record)) {
            final Set<String> excludedSamples = getSamplesToExcludeForStatsBySex(record);
            final Set<String> allSamples = Sets.difference(record.getAllSamples(), excludedSamples);
            final Set<String> carrierSamples = Sets.intersection(record.getCarrierSampleSet(), allSamples);
            final Set<String> backgroundSamples = Sets.difference(allSamples, carrierSamples);
            if (bafCollectionEnabled() && useBafEvidence(record)) {
                // BAF
                final List<BafEvidence> bafEvidence = bafCollector.collectEvidence(record).stream().filter(baf -> allSamples.contains(baf.getSample())).collect(Collectors.toList());
                if (record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.DEL) {
                    final BafHetRatioTester.BafDelResult result = bafHetRatioTester.test(record, bafEvidence, allSamples, carrierSamples, record.getLength());
                    record = bafHetRatioTester.applyToRecord(record, result);
                } else if (record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.DUP) {
                    final BafKolmogorovSmirnovTester.KSTestResult result = bafKolmogorovSmirnovTester.test(record, bafEvidence, carrierSamples);
                    record = bafKolmogorovSmirnovTester.applyToRecord(record, result);
                }
            }
            DiscordantPairEvidenceTester.DiscordantPairTestResult discordantPairResult = null;
            if (discordantPairCollectionEnabled() && useDiscordantPairEvidence(record)) {
                // PE
                final List<DiscordantPairEvidence> discordantPairEvidence = discordantPairCollector.collectEvidence(record);
                discordantPairResult = discordantPairEvidenceTester.test(record, discordantPairEvidence, carrierSamples, backgroundSamples);
                record = discordantPairEvidenceTester.applyToRecord(record, discordantPairResult);
            }
            SplitReadEvidenceTester.SplitReadTestResult splitReadResult = null;
            if (splitReadCollectionEnabled() && useSplitReadEvidence(record)) {
                // SR
                final List<SplitReadEvidence> startSplitReadEvidence = startSplitCollector.collectEvidence(record);
                final List<SplitReadEvidence> endSplitReadEvidence = endSplitCollector.collectEvidence(record);
                splitReadResult = splitReadEvidenceTester.test(record, startSplitReadEvidence, endSplitReadEvidence, carrierSamples, backgroundSamples);
                record = splitReadEvidenceTester.applyToRecord(record, splitReadResult);
            }
            if (discordantPairResult != null && splitReadResult != null) {
                // Combined PE/SR
                final PESREvidenceTester.PESRTestResult result = pesrEvidenceTester.test(splitReadResult, discordantPairResult, carrierSamples, backgroundSamples);
                record = pesrEvidenceTester.applyToRecord(record, result);
            }
        }
        writer.add(SVCallRecordUtils.getVariantBuilder(record).make());
    }

    @Override
    public Object onTraversalSuccess() {
        return super.onTraversalSuccess();
    }

    private VCFHeader getVCFHeader() {
        final VCFHeader header = new VCFHeader(getDefaultToolVCFHeaderLines(), samples);
        header.setSequenceDictionary(dictionary);
        for (final VCFHeaderLine line : getHeaderForVariants().getMetaDataInInputOrder()) {
            header.addMetaDataLine(line);
        }
        if (bafCollectionEnabled()) {
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.BAF_HET_RATIO_ATTRIBUTE, 1, VCFHeaderLineType.Float, "BAF deletion het SNP ratio (mean of 10^(-log-ratio) for carriers)"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.BAF_DEL_LOGLIK_ATTRIBUTE, 1, VCFHeaderLineType.Float, "BAF deletion GMM log-likelihood (negated, higher = more anomalous)"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.BAF_KS_STAT_ATTRIBUTE, 1, VCFHeaderLineType.Float, "BAF Kolmogorov-Smirnov test statistic"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.BAF_KS_Q_ATTRIBUTE, 1, VCFHeaderLineType.Float, "BAF Kolmogorov-Smirnov test quality (-log10 p-value)"));
        }
        if (discordantPairCollectionEnabled()) {
            header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.DISCORDANT_PAIR_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Discordant pair count"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.DISCORDANT_PAIR_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Float, "Discordant pair quality"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.DISCORDANT_PAIR_CARRIER_SIGNAL_ATTRIBUTE, 1, VCFHeaderLineType.Float, "Carrier sample discordant pair signal"));
        }
        if (splitReadCollectionEnabled()) {
            header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.FIRST_SPLIT_READ_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read count at start"));
            header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.SECOND_SPLIT_READ_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read count at end"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.FIRST_SPLIT_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Float, "Split read quality at start"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.SECOND_SPLIT_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Float, "Split read quality at end"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TOTAL_SPLIT_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Float, "Split read quality for both ends"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.FIRST_SPLIT_CARRIER_SIGNAL_ATTRIBUTE, 1, VCFHeaderLineType.Float, "Carrier sample split read signal at start"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.SECOND_SPLIT_CARRIER_SIGNAL_ATTRIBUTE, 1, VCFHeaderLineType.Float, "Carrier sample split read signal at end"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TOTAL_SPLIT_CARRIER_SIGNAL_ATTRIBUTE, 1, VCFHeaderLineType.Float, "Carrier sample split read signal for both ends"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.FIRST_SPLIT_POSITION_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read position at start"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.SECOND_SPLIT_POSITION_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read position at end"));
            if (discordantPairCollectionEnabled()) {
                header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.PESR_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Float, "Combined PE/SR quality"));
                header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.PESR_CARRIER_SIGNAL_ATTRIBUTE, 1, VCFHeaderLineType.Float, "Combined PE/SR carrier signal"));
            }
        }
        return header;
    }

    private boolean splitReadCollectionEnabled() {
        return splitReadsFile != null;
    }

    private boolean discordantPairCollectionEnabled() {
        return discordantPairsFile != null;
    }

    private boolean bafCollectionEnabled() {
        return bafFile != null;
    }
}
