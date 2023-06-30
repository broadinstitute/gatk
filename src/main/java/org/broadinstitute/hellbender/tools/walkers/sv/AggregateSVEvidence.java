package org.broadinstitute.hellbender.tools.walkers.sv;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.io.IOUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.*;
import org.broadinstitute.hellbender.tools.sv.aggregation.*;
import org.broadinstitute.hellbender.tools.sv.cluster.PloidyTable;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.cloud.LocalizedFeatureSubset;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Annotates records with PE/SR/BAF evidence and performs breakpoint refinement with SR.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         SV VCF
 *     </li>
 *     <li>
 *         PE evidence file
 *     </li>
 *     <li>
 *         SR evidence file
 *     </li>
 *     <li>
 *         BAF evidence file
 *     </li>
 *     <li>
 *         Mean depth table
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
 *     gatk AggregateSVEvidence
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */

@CommandLineProgramProperties(
        summary = "Annotate SVs with supporting evidence",
        oneLineSummary = "Annotate SVs with supporting evidence",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@ExperimentalFeature
@DocumentedFeature
public final class AggregateSVEvidence extends TwoPassVariantWalker {
    public static final String SPLIT_READ_LONG_NAME = "split-reads-file";
    public static final String DISCORDANT_PAIRS_LONG_NAME = "discordant-pairs-file";
    public static final String BAF_LONG_NAME = "baf-file";
    public static final String SAMPLE_COVERAGE_LONG_NAME = "sample-coverage";
    public static final String PE_INNER_WINDOW_LONG_NAME = "pe-inner-window";
    public static final String PE_OUTER_WINDOW_LONG_NAME = "pe-outer-window";
    public static final String SR_WINDOW_LONG_NAME = "sr-window";
    public static final String SR_CROSSOVER_LONG_NAME = "sr-insertion-crossover";
    public static final String BAF_MIN_SIZE_LONG_NAME = "baf-min-size";
    public static final String BAF_MAX_SIZE_LONG_NAME = "baf-max-size";
    public static final String BAF_PADDING_FRACTION_LONG_NAME = "baf-padding-fraction";
    public static final String MIN_SNP_CARRIERS_LONG_NAME = "min-snp-carriers";
    public static final String MIN_BAF_COUNT_LONG_NAME = "min-baf-count";
    public static final String P_SNP_LONG_NAME = "p-snp";
    public static final String P_MAX_HOMOZYGOUS_LONG_NAME = "p-max-homozygous";
    public static final String X_CHROMOSOME_LONG_NAME = "x-chromosome-name";
    public static final String Y_CHROMOSOME_LONG_NAME = "y-chromosome-name";

    @Argument(
            doc = "Split reads evidence file",
            fullName = SPLIT_READ_LONG_NAME,
            optional = true
    )
    private GATKPath splitReadsFile;

    @Argument(
            doc = "Discordant pairs evidence file",
            fullName = DISCORDANT_PAIRS_LONG_NAME,
            optional = true
    )
    private GATKPath discordantPairsFile;

    @Argument(
            doc = "B-allele frequency (BAF) evidence file",
            fullName = BAF_LONG_NAME,
            optional = true
    )
    private GATKPath bafFile;

    @Argument(
            doc = "Tab-delimited table with sample IDs in the first column and expected per-base coverage per sample " +
                    "in the second column.",
            fullName = SAMPLE_COVERAGE_LONG_NAME
    )
    private GATKPath sampleCoverageFile;

    @Argument(
            doc = "Output VCF",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private String outputFile;

    @Argument(
            doc = "Inner discordant pair window size (bp)",
            fullName = PE_INNER_WINDOW_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private int innerWindow = 50;

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
    private int splitReadWindow = 200;

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
    private int bafMinSize = 0;

    @Argument(
            doc = "Maximum variant size for BAF aggregation",
            fullName = BAF_MAX_SIZE_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private int bafMaxSize = 10000000;

    @Argument(
            doc = "BAF flanking region size as a fraction of variant size",
            fullName = BAF_PADDING_FRACTION_LONG_NAME,
            minValue = 0.01,
            optional = true
    )
    private double bafPaddingFraction = 1.0;

    @Argument(
            doc = "Minimum number of SNP carriers, used for filtering loci during duplication BAF assessment.",
            fullName = MIN_SNP_CARRIERS_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private int minSnpCarriers = 5;

    @Argument(
            doc = "Minimum number of BAF values required in carrier and non-carrier groups.",
            fullName = MIN_BAF_COUNT_LONG_NAME,
            minValue = 1,
            optional = true
    )
    private int minBafCount = 2;

    @Argument(
            doc = "Baseline expected SNPs per locus, used for filtering deletions in likely regions of homozygosity " +
                    "during BAF assessment.",
            fullName = P_SNP_LONG_NAME,
            minValue = 0.,
            maxValue = 1.,
            optional = true
    )
    private double pSnp = 0.001;

    @Argument(
            doc = "Significance threshold for binomial test on region homozygosity, used for filtering deletions in " +
                    "likely regions of homozygosity during BAF assessment.",
            fullName = P_MAX_HOMOZYGOUS_LONG_NAME,
            minValue = 0.,
            maxValue = 1.,
            optional = true
    )
    private double pMaxHomozygous = 0.05;

    /**
     * Expected format is tab-delimited and contains a header with the first column SAMPLE and remaining columns
     * contig names. Each row corresponds to a sample, with the sample ID in the first column and contig ploidy
     * integers in their respective columns.
     */
    @Argument(
            doc = "Sample ploidy table (.tsv). Required only if the input VCF contains allosomal records.",
            fullName = SVCluster.PLOIDY_TABLE_LONG_NAME,
            optional = true
    )
    private GATKPath ploidyTablePath;

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

    @Argument(
            doc = "Output compression level",
            fullName = PrintSVEvidence.COMPRESSION_LEVEL_NAME,
            minValue = 0, maxValue = 9, optional = true
    )
    private int compressionLevel = 4;

    private SAMSequenceDictionary dictionary;
    private VariantContextWriter writer;

    private FeatureDataSource<SplitReadEvidence> splitReadSource;
    private SplitReadEvidenceAggregator startSplitCollector;
    private SplitReadEvidenceAggregator endSplitCollector;
    private SplitReadEvidenceTester splitReadEvidenceTester;

    private FeatureDataSource<DiscordantPairEvidence> discordantPairSource;
    private DiscordantPairEvidenceAggregator discordantPairCollector;
    private DiscordantPairEvidenceTester discordantPairEvidenceTester;

    private PESREvidenceTester pesrEvidenceTester;

    private FeatureDataSource<BafEvidence> bafSource;
    private BafEvidenceAggregator bafCollector;
    private BafHetRatioTester bafHetRatioTester;
    private BafKolmogorovSmirnovTester bafKolmogorovSmirnovTester;

    private Map<String,Double> sampleCoverageMap;
    private Set<String> samples;
    private VCFHeader header;
    private Set<String> maleSamples;
    private Set<String> femaleSamples;

    private Collection<SimpleInterval> discordantPairIntervals;
    private Collection<SimpleInterval> startSplitReadIntervals;
    private Collection<SimpleInterval> endSplitReadIntervals;
    private Collection<SimpleInterval> bafIntervals;

    private LocalizedFeatureSubset<SplitReadEvidence> localizedEndSplitReads;

    private Collection<VariantContext> outputBuffer;

    private final int BAF_QUERY_LOOKAHEAD = 0;
    private final int SPLIT_READ_QUERY_LOOKAHEAD = 0;
    private final int DISCORDANT_PAIR_QUERY_LOOKAHEAD = 0;

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
        loadSampleCoverage();
        if (discordantPairCollectionEnabled()) {
            initializeDiscordantPairCollection();
        }
        if (splitReadCollectionEnabled()) {
            initializeSplitReadCollection();
            if (discordantPairCollectionEnabled()) {
                pesrEvidenceTester = new PESREvidenceTester(sampleCoverageMap, dictionary);
            }
        }
        if (bafCollectionEnabled()) {
            initializeBAFCollection();
        }
        discordantPairIntervals = new ArrayList<>();
        startSplitReadIntervals = new ArrayList<>();
        endSplitReadIntervals = new ArrayList<>();
        bafIntervals = new ArrayList<>();
        outputBuffer = new ArrayList<>();
        writer = createVCFWriter(Paths.get(outputFile));
        header = getVCFHeader();
        writer.writeHeader(header);
        if (ploidyTablePath != null) {
            initializeSampleSexSets();
        }
    }

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
    }

    private void initializeBAFCollection() {
        initializeBAFEvidenceDataSource();
        bafCollector = new BafEvidenceAggregator(bafSource, dictionary, bafPaddingFraction);
        bafHetRatioTester = new BafHetRatioTester(pSnp, pMaxHomozygous);
        bafKolmogorovSmirnovTester = new BafKolmogorovSmirnovTester(minSnpCarriers, minBafCount);
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

    private void loadSampleCoverage() {
        final String fileString = sampleCoverageFile.toString();
        try {
            sampleCoverageMap = IOUtils.readLines(BucketUtils.openFile(fileString), Charset.defaultCharset()).stream()
                    .map(line -> line.split("\t"))
                    .collect(Collectors.toMap(tokens -> tokens[0], tokens -> Double.valueOf(tokens[1])));
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(fileString, e);
        }
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

    private boolean validRecordType(final SVCallRecord call) {
        return call.getType() != GATKSVVCFConstants.StructuralVariantAnnotationType.CNV;
    }

    private boolean useDiscordantPairEvidence(final SVCallRecord call) {
        return !call.isDepthOnly() && call.getType() != GATKSVVCFConstants.StructuralVariantAnnotationType.INS;
    }

    private boolean useSplitReadEvidence(final SVCallRecord call) {
        return !call.isDepthOnly();
    }

    private boolean useBafEvidence(final SVCallRecord call) {
        final Integer length = call.getLength();
        return (call.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.DEL || call.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.DUP)
                && length != null && length >= bafMinSize && length <= bafMaxSize;
    }

    // Collect intervals for the evidence
    @Override
    public void firstPassApply(final VariantContext variant, final ReadsContext readsContext,
                      final ReferenceContext referenceContext, final FeatureContext featureContext) {
        final SVCallRecord call = SVCallRecordUtils.create(variant);
        if (!validRecordType(call)) {
            return;
        }
        if (discordantPairCollectionEnabled() && useDiscordantPairEvidence(call)) {
            discordantPairIntervals.add(discordantPairCollector.getEvidenceQueryInterval(call));
        }
        if (splitReadCollectionEnabled() && useSplitReadEvidence(call)) {
            startSplitReadIntervals.add(SplitReadEvidenceAggregator.getStartEvidenceQueryInterval(call, splitReadWindow, dictionary));
            endSplitReadIntervals.add(SplitReadEvidenceAggregator.getEndEvidenceQueryInterval(call, splitReadWindow, dictionary));
        }
        if (bafCollectionEnabled() && useBafEvidence(call)) {
            bafIntervals.add(bafCollector.getEvidenceQueryInterval(call));
        }
    }

    // Localize evidence over the intervals we found
    @Override
    public void afterFirstPass() {
        if (discordantPairCollectionEnabled()) {
            discordantPairCollector.setCacheIntervals(discordantPairIntervals);
        }
        if (splitReadCollectionEnabled()) {
            startSplitCollector = new SplitReadEvidenceAggregator(splitReadSource, dictionary, splitReadWindow, true);
            startSplitCollector.setCacheIntervals(startSplitReadIntervals);
            logger.info("Localizing split read evidence at END loci...");
            localizedEndSplitReads = new LocalizedFeatureSubset<>(splitReadsFile, endSplitReadIntervals, "endSplitReads", 0,
                    cloudPrefetchBuffer, cloudIndexPrefetchBuffer, new ArrayList<>(samples), compressionLevel, dictionary);
            logger.info("Finished localization of " + localizedEndSplitReads.getRecordCount() + " feature records");
            // Does not use interval caching to avoid thrashing due distant ENDs
            endSplitCollector = new SplitReadEvidenceAggregator(localizedEndSplitReads.getLocalizedDataSource(), dictionary, splitReadWindow, false);
        }
        if (bafCollectionEnabled()) {
            bafCollector.setCacheIntervals(bafIntervals);
        }
    }

    /**
     * Returns set of samples to exclude for evidence stat calculations by sex.
     * @param record
     * @return set of sample ids or null
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

    // Perform evidence testing
    @Override
    public void secondPassApply(final VariantContext variant, final ReadsContext readsContext,
                                final ReferenceContext referenceContext, final FeatureContext featureContext) {
        SVCallRecord record = SVCallRecordUtils.create(variant);
        if (validRecordType(record)) {
            final Set<String> excludedSamples = getSamplesToExcludeForStatsBySex(record);
            flushOutputBuffer(record.getPositionAInterval());
            final Set<String> allSamples = Sets.difference(record.getAllSamples(), excludedSamples);
            final Set<String> carrierSamples = Sets.intersection(record.getCarrierSampleSet(), allSamples);
            final Set<String> backgroundSamples = Sets.difference(allSamples, carrierSamples);
            if (bafCollectionEnabled() && useBafEvidence(record)) {
                // BAF
                final List<BafEvidence> bafEvidence = bafCollector.collectEvidence(record).stream().filter(baf -> allSamples.contains(baf.getSample())).collect(Collectors.toList());
                if (record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.DEL) {
                    final Double result = bafHetRatioTester.test(record, bafEvidence, allSamples, carrierSamples, (int) bafPaddingFraction * record.getLength());
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
        outputBuffer.add(SVCallRecordUtils.getVariantBuilder(record).make());
    }

    @Override
    public Object onTraversalSuccess() {
        outputBuffer.stream().sorted(IntervalUtils.getDictionaryOrderComparator(dictionary)).forEach(writer::add);
        outputBuffer.clear();
        return super.onTraversalSuccess();
    }

    private void flushOutputBuffer(final SimpleInterval currentLocus) {
        outputBuffer.stream()
                .filter(v -> !variantIsActive(v, currentLocus))
                .sorted(IntervalUtils.getDictionaryOrderComparator(dictionary))
                .forEach(writer::add);
        outputBuffer = outputBuffer.stream()
                .filter(v -> variantIsActive(v, currentLocus))
                .collect(Collectors.toList());
    }

    private boolean variantIsActive(final VariantContext variant, final SimpleInterval currentLocus) {
        return variant.getContig().equals(currentLocus.getContig()) && variant.getStart() >= currentLocus.getStart() - splitReadWindow;
    }

    private VCFHeader getVCFHeader() {
        final VCFHeader header = new VCFHeader(getDefaultToolVCFHeaderLines(), samples);
        header.setSequenceDictionary(dictionary);
        for (final VCFHeaderLine line : getHeaderForVariants().getMetaDataInInputOrder()) {
            header.addMetaDataLine(line);
        }
        if (bafCollectionEnabled()) {
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.BAF_HET_RATIO_ATTRIBUTE, 1, VCFHeaderLineType.Float, "Log ratio of non-carrier to carrier het count"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.BAF_KS_STAT_ATTRIBUTE, 1, VCFHeaderLineType.Float, "BAF Kolmogorov-Smirnov test statistic"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.BAF_KS_Q_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "BAF Kolmogorov-Smirnov test phred-scaled p-value"));
        }
        if (discordantPairCollectionEnabled()) {
            header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.DISCORDANT_PAIR_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Discordant pair count"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.DISCORDANT_PAIR_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Discordant pair quality"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.DISCORDANT_PAIR_CARRIER_SIGNAL_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Carrier sample discordant pair signal"));
        }
        if (splitReadCollectionEnabled()) {
            header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.FIRST_SPLIT_READ_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read count at start"));
            header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.SECOND_SPLIT_READ_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read count at end"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.FIRST_SPLIT_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read quality at start"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.SECOND_SPLIT_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read quality at end"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TOTAL_SPLIT_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read quality for both ends"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.FIRST_SPLIT_CARRIER_SIGNAL_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Carrier sample split read signal at start"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.SECOND_SPLIT_CARRIER_SIGNAL_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Carrier sample split read signal at end"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TOTAL_SPLIT_CARRIER_SIGNAL_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Carrier sample split read signal for both ends"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.FIRST_SPLIT_POSITION_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read position at start"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.SECOND_SPLIT_POSITION_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read position at end"));
            if (discordantPairCollectionEnabled()) {
                header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.PESR_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Combined PE/SR quality"));
                header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.PESR_CARRIER_SIGNAL_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Combined PE/SR carrier signal"));
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
