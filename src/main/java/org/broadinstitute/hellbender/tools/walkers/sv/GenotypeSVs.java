package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.*;
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
import org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVCollapser;
import org.broadinstitute.hellbender.tools.sv.cluster.PloidyTable;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterWalker;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.codecs.DepthEvidenceCodec;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.io.IOException;
import java.nio.charset.Charset;
import java.util.*;

/**
 * <p>Genotypes SVs.</p>
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
 *     <li>
 *         Depth and PE/SR genotyping tables generated with TrainSVGenotyping
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Genotyped VCF
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk GenotypeSVs \
 *      -V variants.vcf.gz \
 *      -O genotyped.vcf.gz \
*       --median-coverage median_coverage.tsv \
 *      --rd-file samples.rd.txt.gz \
 *      --discordant-pairs-file samples.pe.txt.gz \
 *      --split-reads-file samples.sr.txt.gz \
 *      --sequence-dictionary Homo_sapiens_assembly38.dict \
 *      --ploidy-table ploidy_table.tsv \
 *      --pesr-exclusion-intervals pesr_list.bed \
 *      --depth-exclusion-intervals depth_list.bed \
 *      --rd-depth-table samples.rd_depth_geno_params.tsv \
 *      --rd-pesr-table samples.rd_pesr_geno_params.tsv \
 *      --pe-table samples.pe_geno_params.tsv \
 *      --sr-table samples.sr_geno_params.tsv
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        summary = "Genotypes SVs",
        oneLineSummary = "Genotypes SVs",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
public final class GenotypeSVs extends MultiplePassVariantWalker {

    public static final String DEPTH_CUTOFFS_TABLE_LONG_NAME = "rd-depth-table";
    public static final String PESR_DEPTH_CUTOFFS_TABLE_LONG_NAME = "rd-pesr-table";
    public static final String DISCORDANT_PAIR_CUTOFFS_TABLE_LONG_NAME = "pe-table";
    public static final String SPLIT_READ_CUTOFFS_TABLE_LONG_NAME = "sr-table";
    public static final String SPLIT_READ_INSERTION_MEDIAN_HOM_LONG_NAME = "sr-median-hom-ins";
    public static final String SPLIT_READ_HOM_CUTOFF_MULTIPLIER_LONG_NAME = "sr-hom-cutoff-multiplier";

    @Argument(
            fullName = TrainSVGenotyping.DEPTH_EVIDENCE_FILE_PATH_LONG_NAME,
            doc = "Indexed read counts file ending with " + DepthEvidenceCodec.FORMAT_SUFFIX + ".gz"
    )
    public GATKPath depthEvidenceFile;

    @Argument(
            fullName = TrainSVGenotyping.MEDIAN_COUNTS_FILE_PATH_LONG_NAME,
            doc = "Median counts file"
    )
    public GATKPath medianFile;

    @Argument(
            fullName = DEPTH_CUTOFFS_TABLE_LONG_NAME,
            doc = "Depth genotyping cutoffs table for depth-only variants"
    )
    public GATKPath depthCutoffsTablePath;

    @Argument(
            fullName = PESR_DEPTH_CUTOFFS_TABLE_LONG_NAME,
            doc = "Depth genotyping cutoffs table for PESR variants"
    )
    public GATKPath pesrDepthCutoffsTablePath;

    @Argument(
            doc = "Intervals to exclude for depth evidence genotyping",
            fullName = TrainSVGenotyping.DEPTH_EXCLUSION_INTERVALS_LONG_NAME,
            optional = true
    )
    protected GATKPath depthExclusionIntervalsPath;

    @Argument(
            doc = "Intervals to exclude for PE/SR genotyping",
            fullName = TrainSVGenotyping.PESR_EXCLUSION_INTERVALS_LONG_NAME,
            optional = true
    )
    protected GATKPath pesrExclusionIntervalsPath;

    @Argument(
            fullName = DISCORDANT_PAIR_CUTOFFS_TABLE_LONG_NAME,
            doc = "Discordant pair genotyping cutoffs table",
            optional = true
    )
    public GATKPath discordantPairCutoffsTablePath;

    @Argument(
            fullName = SPLIT_READ_CUTOFFS_TABLE_LONG_NAME,
            doc = "Split read genotyping cutoffs table",
            optional = true
    )
    public GATKPath splitReadCutoffsTablePath;

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
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output VCF"
    )
    public GATKPath outputVcf;

    @Argument(
            fullName = SPLIT_READ_INSERTION_MEDIAN_HOM_LONG_NAME,
            doc = "Median homozygous split read count for insertions",
            minValue = 1
    )
    public int medianHomIns = 105;

    @Argument(
            fullName = SPLIT_READ_HOM_CUTOFF_MULTIPLIER_LONG_NAME,
            doc = "Median homozygous split read cutoff multiplier",
            minValue = 1
    )
    public double medianHomCutoffMultiplier = 1.6;

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
    public int numBins = 100000;

    @Argument(
            fullName = TrainSVGenotyping.MIN_PE_QUALITY_LONG_NAME,
            doc = "Discordant pair quality threshold, used for setting evidence count threshold",
            minValue = 1
    )
    public double minDiscordantPairQuality = 30;

    @Argument(
            fullName = TrainSVGenotyping.MIN_SR_QUALITY_LONG_NAME,
            doc = "Split read quality threshold, used for setting evidence count threshold",
            minValue = 1
    )
    public double minSplitReadQuality = 30;

    @Argument(
            fullName = AggregateDepthEvidence.MAX_QUALITY_LONG_NAME,
            doc = "Max quality score",
            minValue = 1
    )
    public int maxQual = 999;

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
    private DepthEvidenceGenotyper depthOnlyGenotyper;
    private DepthEvidenceGenotyper pesrDepthGenotyper;
    private List<String> masterSampleList;
    private FeatureDataSource<DepthEvidence> depthSource;

    private FeatureDataSource<DiscordantPairEvidence> discordantPairSource;
    private DiscordantPairEvidenceAggregator discordantPairCollector;
    private DiscordantPairEvidenceGenotyper discordantPairGenotyper;
    private DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeParameters discordantPairParameters;

    private FeatureDataSource<SplitReadEvidence> splitReadSource;
    private SplitReadEvidenceAggregator splitReadStartCollector;
    private SplitReadEvidenceAggregator splitReadEndCollector;
    private SplitReadEvidenceGenotyper splitReadGenotyper;
    private SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics splitReadParameters;

    private OverlapDetector<SimpleInterval> depthExclusionIntervals;

    private static final int DISCORDANT_PAIR_QUERY_LOOKAHEAD = 0;
    private static final int SPLIT_READ_QUERY_LOOKAHEAD = 0;

    protected int numberOfPasses() {
        return 1;
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
        depthSource = new FeatureDataSource<>(depthEvidenceFile.toString());
        sampleMedians = loadMedianSampleCoverageTable();
        dictionary = getBestAvailableSequenceDictionary();

        if (depthExclusionIntervalsPath != null) {
            final GenomeLocParser genomeLocParser = new GenomeLocParser(dictionary);
            final GenomeLocSortedSet genomeLocs = IntervalUtils.loadIntervals(Collections.singletonList(depthExclusionIntervalsPath.toString()), IntervalSetRule.UNION, IntervalMergingRule.OVERLAPPING_ONLY, 0, genomeLocParser);
            depthExclusionIntervals = OverlapDetector.create(genomeLocs.stream().map(SimpleInterval::new).toList());
        }

        loader = new DepthMatrixLoader(depthSource, numBins, largeVariantSize, largeVariantPoints, largeVariantWindow, depthExclusionIntervals, dictionary);
        writer = createVCFWriter(outputVcf);
        outputHeader = createHeader(getHeaderForVariants());
        writer.writeHeader(outputHeader);
        masterSampleList = outputHeader.getSampleNamesInOrder();
        try (final TableReader<DepthEvidenceGenotyper.CopyStateStats> reader = TableUtils.reader(depthCutoffsTablePath.toPath(), new DepthEvidenceGenotyper.DepthTableParser()::tableParser)) {
            depthOnlyGenotyper = new DepthEvidenceGenotyper(reader.toList(), masterSampleList, maxQual, dictionary);
        } catch (final IOException e) {
            throw new GATKException("Error while reading depth-only cutoffs table", e);
        }
        try (final TableReader<DepthEvidenceGenotyper.CopyStateStats> reader = TableUtils.reader(pesrDepthCutoffsTablePath.toPath(), new DepthEvidenceGenotyper.DepthTableParser()::tableParser)) {
            pesrDepthGenotyper = new DepthEvidenceGenotyper(reader.toList(), masterSampleList, maxQual, dictionary);
        } catch (final IOException e) {
            throw new GATKException("Error while reading PESR depth cutoffs table", e);
        }
        try (final TableReader<DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeParameters> reader = TableUtils.reader(discordantPairCutoffsTablePath.toPath(), new DiscordantPairEvidenceGenotyper.DiscordantPairTableParser()::tableParser)) {
            discordantPairParameters = reader.readRecord();
        } catch (final IOException e) {
            throw new GATKException("Error while reading discordant pair cutoffs table", e);
        }
        try (final TableReader<SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics> reader = TableUtils.reader(splitReadCutoffsTablePath.toPath(), new SplitReadEvidenceGenotyper.SplitReadTableParser()::tableParser)) {
            splitReadParameters = reader.readRecord();
        } catch (final IOException e) {
            throw new GATKException("Error while reading split read cutoffs table", e);
        }
        if (splitReadCollectionEnabled()) {
            Utils.validate(discordantPairCollectionEnabled(), "Discordant pairs file must be provided for split read training");
            Utils.validateArg(splitReadCutoffsTablePath != null, "Must provide split read cutoffs table with split read evidence file");
        }
        if (discordantPairCollectionEnabled()) {
            Utils.validateArg(discordantPairCutoffsTablePath != null, "Must provide discordant pair cutoffs table with discordant pair evidence file");
        }

        initializeDiscordantPairCollection();
        initializeSplitReadCollection();
        discordantPairGenotyper = new DiscordantPairEvidenceGenotyper(sampleMedians, minDiscordantPairQuality, null, PESREvidenceTester.DEPTH_BASIS, maxQual);
        splitReadGenotyper = new SplitReadEvidenceGenotyper(sampleMedians, masterSampleList.size(), minSplitReadQuality, null, PESREvidenceTester.DEPTH_BASIS, maxQual);

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
            final DepthEvidenceGenotyper.DepthGenotypeResult depth = genotypeDepth(record);
            final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult discordantPair = genotypeDiscordantPairs(record);
            final SplitReadEvidenceGenotyper.SplitReadGenotypeResult splitRead = genotypeSplitReads(record);
            final Double combinedVariantQual = Math.max(combineVariantQuals(record, depth, discordantPair, splitRead), 1); // VARQ > 0
            final List<GenotypeAndQuality> combinedGenotypes = combineGenotypes(record, depth, discordantPair, splitRead);
            writeGenotypes(record, depth, discordantPair, splitRead, combinedVariantQual, combinedGenotypes);
        }
    }

    private record GenotypeAndQuality(int genotype, int quality, List<GATKSVVCFConstants.EvidenceTypes> evidence) {}

    private List<GenotypeAndQuality> combineGenotypes(final SVCallRecord record, final DepthEvidenceGenotyper.DepthGenotypeResult depthGenotypeResult,
                                       final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult discordantPairGenotypeResult,
                                       final SplitReadEvidenceGenotyper.SplitReadGenotypeResult splitReadGenotypeResult) {
        final int numSamples = masterSampleList.size();
        final List<GenotypeAndQuality> result = new ArrayList<>(numSamples);
        for (int i = 0; i < numSamples; i++) {
            if (record.isDepthOnly()) {
                final int genotype = genotypeStateFromCopyState(depthGenotypeResult.copyStates()[i], record.getType());
                final double depthQuality = depthGenotypeResult.genotypeQuals()[i];
                result.add(new GenotypeAndQuality(genotype, (int) depthQuality, Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.RD)));
            } else {
                final int peGenotype = discordantPairGenotypeResult.genotypes()[i];
                final int peQuality = discordantPairGenotypeResult.genotypeQuals()[i];
                final int srGenotype = splitReadGenotypeResult.genotypes()[i];
                final int srQuality = splitReadGenotypeResult.genotypeQuals()[i];
                final GenotypeAndQuality pesrResult = getPESRGenotype(peGenotype, peQuality, srGenotype, srQuality);
                if (!record.isSimpleCNV()) {
                    result.add(pesrResult);
                } else {
                    final int depthCopyState = depthGenotypeResult.copyStates()[i];
                    final double depthQuality = depthGenotypeResult.genotypeQuals()[i];
                    final int peCopyState;
                    final int srCopyState;
                    final int pesrCopyState;
                    if (record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.DEL) {
                        peCopyState = Math.max(0, 2 - peGenotype);
                        srCopyState = Math.max(0, 2 - srGenotype);
                        pesrCopyState = Math.max(0, 2 - pesrResult.genotype);
                    } else {
                        peCopyState = 2 + peGenotype;
                        srCopyState = 2 + srGenotype;
                        pesrCopyState = 2 + pesrResult.genotype;
                    }
                    if (record.getLength() != null && record.getLength() >= 5000) {
                        final List<GATKSVVCFConstants.EvidenceTypes> evidence;
                        if (depthCopyState == peCopyState && depthCopyState == srCopyState) {
                            evidence = List.of(GATKSVVCFConstants.EvidenceTypes.RD, GATKSVVCFConstants.EvidenceTypes.PE, GATKSVVCFConstants.EvidenceTypes.SR);
                        } else if (depthCopyState == 2) {
                            if (depthCopyState == peCopyState) {
                                    evidence = List.of(GATKSVVCFConstants.EvidenceTypes.RD, GATKSVVCFConstants.EvidenceTypes.PE);
                            } else if (depthCopyState == srCopyState) {
                                evidence = List.of(GATKSVVCFConstants.EvidenceTypes.RD, GATKSVVCFConstants.EvidenceTypes.SR);
                            } else {
                                evidence = List.of(GATKSVVCFConstants.EvidenceTypes.RD);
                            }
                        } else {
                            // depthCopyState != 2
                            if (peCopyState != 2) {
                                if (srCopyState != 2) {
                                    evidence = List.of(GATKSVVCFConstants.EvidenceTypes.RD, GATKSVVCFConstants.EvidenceTypes.PE, GATKSVVCFConstants.EvidenceTypes.SR);
                                } else {
                                    evidence = List.of(GATKSVVCFConstants.EvidenceTypes.RD, GATKSVVCFConstants.EvidenceTypes.PE);
                                }
                            } else if (srCopyState != 2) {
                                evidence = List.of(GATKSVVCFConstants.EvidenceTypes.RD, GATKSVVCFConstants.EvidenceTypes.SR);
                            } else {
                                evidence = Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.RD);
                            }
                        }
                        final int quality = adjustQuality((int) depthQuality, depthCopyState, pesrResult.quality, pesrCopyState);
                        final int genotype = genotypeStateFromCopyState(depthCopyState, record.getType());
                        result.add(new GenotypeAndQuality(genotype, quality, evidence));
                    } else if (record.getLength() != null && record.getLength() >= 1000) {
                        final List<GATKSVVCFConstants.EvidenceTypes> evidence = new ArrayList<>(3);
                        if (depthCopyState == peCopyState && depthCopyState == srCopyState) {
                            evidence.add(GATKSVVCFConstants.EvidenceTypes.RD);
                            evidence.add(GATKSVVCFConstants.EvidenceTypes.PE);
                            evidence.add(GATKSVVCFConstants.EvidenceTypes.SR);
                        } else if (pesrCopyState == 2) {
                            if (depthCopyState == 2) {
                                evidence.add(GATKSVVCFConstants.EvidenceTypes.RD);
                            }
                            if (peCopyState == 2) {
                                evidence.add(GATKSVVCFConstants.EvidenceTypes.PE);
                            }
                            if (srCopyState == 2) {
                                evidence.add(GATKSVVCFConstants.EvidenceTypes.SR);
                            }
                        } else {
                            // pesrCopyState != 2
                            if (depthCopyState != 2) {
                                evidence.add(GATKSVVCFConstants.EvidenceTypes.RD);
                            }
                            if (peCopyState != 2) {
                                evidence.add(GATKSVVCFConstants.EvidenceTypes.PE);
                            }
                            if (srCopyState != 2) {
                                evidence.add(GATKSVVCFConstants.EvidenceTypes.SR);
                            }
                        }
                        final int quality = adjustQuality(pesrResult.quality, pesrCopyState, (int) depthQuality, depthCopyState);
                        result.add(new GenotypeAndQuality(pesrResult.genotype, quality, evidence));
                    } else {
                        final List<GATKSVVCFConstants.EvidenceTypes> evidence = new ArrayList<>(3);
                        if (depthCopyState == peCopyState && depthCopyState == srCopyState) {
                            evidence.add(GATKSVVCFConstants.EvidenceTypes.RD);
                            evidence.add(GATKSVVCFConstants.EvidenceTypes.PE);
                            evidence.add(GATKSVVCFConstants.EvidenceTypes.SR);
                        } else if (pesrCopyState == 2) {
                            if (depthCopyState == 2) {
                                evidence.add(GATKSVVCFConstants.EvidenceTypes.RD);
                            }
                            if (peCopyState == 2) {
                                evidence.add(GATKSVVCFConstants.EvidenceTypes.PE);
                            }
                            if (srCopyState == 2) {
                                evidence.add(GATKSVVCFConstants.EvidenceTypes.SR);
                            }
                        } else {
                            // pesrCopyState != 2
                            if (depthCopyState != 2 && peCopyState != 2 && srCopyState != 2) {
                                evidence.add(GATKSVVCFConstants.EvidenceTypes.RD);
                                evidence.add(GATKSVVCFConstants.EvidenceTypes.PE);
                                evidence.add(GATKSVVCFConstants.EvidenceTypes.SR);
                            } else {
                                if (pesrCopyState == depthCopyState) {
                                    if (peCopyState != 2) {
                                        evidence.add(GATKSVVCFConstants.EvidenceTypes.RD);
                                        evidence.add(GATKSVVCFConstants.EvidenceTypes.PE);
                                    } else if (srCopyState != 2) {
                                        evidence.add(GATKSVVCFConstants.EvidenceTypes.RD);
                                        evidence.add(GATKSVVCFConstants.EvidenceTypes.SR);
                                    }
                                } else {
                                    if (peCopyState == 2 && srCopyState != 2) {
                                        evidence.add(GATKSVVCFConstants.EvidenceTypes.SR);
                                    } else if (peCopyState != 2 && srCopyState == 2) {
                                        evidence.add(GATKSVVCFConstants.EvidenceTypes.PE);
                                    } else {
                                        evidence.add(GATKSVVCFConstants.EvidenceTypes.PE);
                                        evidence.add(GATKSVVCFConstants.EvidenceTypes.SR);
                                    }
                                }
                            }
                        }
                        final int quality;
                        if (pesrResult.quality == maxQual) {
                            quality = pesrResult.quality;
                        } else if (peCopyState == srCopyState) {
                            if (peQuality >= srQuality) {
                                // TODO should this use peQuality?
                                quality = (int) (pesrResult.quality + (0.5 * (maxQual - pesrResult.quality) * (srQuality / (double) maxQual)));
                            } else {
                                quality = (int) (pesrResult.quality + (0.5 * (maxQual - pesrResult.quality) * (peQuality / (double) maxQual)));
                            }
                        } else {
                            quality = pesrResult.quality;
                        }
                        result.add(new GenotypeAndQuality(pesrResult.genotype, quality, evidence));
                    }
                }
            }
        }
        return result;
    }

    private int genotypeStateFromCopyState(final int copyState, final GATKSVVCFConstants.StructuralVariantAnnotationType svType) {
        if (svType == GATKSVVCFConstants.StructuralVariantAnnotationType.DEL) {
            return Math.max(0, 2 - copyState);
        } else {
            return Math.max(0, copyState - 2);
        }
    }

    private int adjustQuality(final int qualityA, final int stateA, final int qualityB, final int stateB) {
        if (qualityA == maxQual) {
            return qualityA;
        } else if (stateA == stateB) {
            return (int) (qualityA + (0.5 * (maxQual - qualityA) * (qualityB / (double) maxQual)));
        } else {
            return (int) (qualityA - (0.5 * (qualityA - 1) * (qualityB / (double) maxQual)));
        }
    }

    private static GenotypeAndQuality getPESRGenotype(final int peState, final int peQuality, final int srState, final int srQuality) {
        if (peState == srState) {
            return new GenotypeAndQuality(peState, Math.max(peQuality, srQuality), List.of(GATKSVVCFConstants.EvidenceTypes.PE, GATKSVVCFConstants.EvidenceTypes.SR));
        } else if (peState > 0 && srState == 0) {
            return new GenotypeAndQuality(peState, peQuality, Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.PE));
        } else if (peState == 0 && srState > 0) {
            return new GenotypeAndQuality(srState, srQuality, Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.SR));
        } else if (peQuality >= srQuality) {
            return new GenotypeAndQuality(peState, peQuality, List.of(GATKSVVCFConstants.EvidenceTypes.PE, GATKSVVCFConstants.EvidenceTypes.SR));
        } else {
            return new GenotypeAndQuality(srState, srQuality, List.of(GATKSVVCFConstants.EvidenceTypes.PE, GATKSVVCFConstants.EvidenceTypes.SR));
        }
    }

    private double combineVariantQuals(final SVCallRecord record, final DepthEvidenceGenotyper.DepthGenotypeResult depthGenotypeResult,
                                       final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult discordantPairGenotypeResult,
                                       final SplitReadEvidenceGenotyper.SplitReadGenotypeResult splitReadGenotypeResult) {
        final Double depthQual = depthGenotypeResult == null ? null : depthGenotypeResult.variantQual();
        final Double peQual = discordantPairGenotypeResult == null ? null : discordantPairGenotypeResult.variantQual();
        final Double srQual = splitReadGenotypeResult == null ? null : splitReadGenotypeResult.variantQual();
        if (record.isDepthOnly()) {
            return depthQual;
        } else {
            Utils.nonNull(discordantPairGenotypeResult);
            Utils.nonNull(splitReadGenotypeResult);
            if (!record.isSimpleCNV()) {
                return Math.max(peQual, srQual);
            } else if (record.getLength() >= 5000) {
                Utils.nonNull(depthGenotypeResult);
                if (depthQual == maxQual) {
                    return depthQual;
                } else if (peQual >= srQual) {
                    return depthQual + ((maxQual - depthQual) * (peQual / maxQual) * 0.5);
                } else {
                    return depthQual + ((maxQual - depthQual) * (srQual / maxQual) * 0.5);
                }
            } else {
                // Under 5kbp
                if (peQual == maxQual || srQual == maxQual) {
                    return maxQual;
                } else if (peQual >= srQual) {
                    return peQual + ((maxQual - peQual) * (depthQual / maxQual) * 0.5);
                } else {
                    return srQual + ((maxQual - srQual) * (depthQual / maxQual) * 0.5);
                }
            }
        }
    }

    public void writeGenotypes(final SVCallRecord record, final DepthEvidenceGenotyper.DepthGenotypeResult depthGenotypeResult,
                               final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult discordantPairGenotypeResult,
                               final SplitReadEvidenceGenotyper.SplitReadGenotypeResult splitReadGenotypeResult,
                               final Double combinedVariantQual,
                               final List<GenotypeAndQuality> combinedGenotypes) {
        final ArrayList<Genotype> newGenotypeList = new ArrayList<>(masterSampleList.size());
        int maxGenotype = 0;
        for (int i = 0; i < masterSampleList.size(); i++) {
            final String sample = masterSampleList.get(i);
            final GenotypeBuilder builder = new GenotypeBuilder(sample);
            builder.attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, ploidyTable.get(sample, record.getContigA()));
            if (depthGenotypeResult != null) {
                builder.attribute(GATKSVVCFConstants.DEPTH_GENOTYPE_COPY_NUMBER_FORMAT, depthGenotypeResult.copyStates()[i]);
                builder.attribute(GATKSVVCFConstants.DEPTH_MEDIAN_COPY_RATIO, depthGenotypeResult.sampleDepths()[i]);
                builder.attribute(GATKSVVCFConstants.DEPTH_GENOTYPE_QUALITY_ATTRIBUTE, rescaleGq((int) depthGenotypeResult.genotypeQuals()[i]));
            }
            if (discordantPairGenotypeResult != null) {
                builder.attribute(GATKSVVCFConstants.DISCORDANT_PAIR_GENOTYPE_ATTRIBUTE, discordantPairGenotypeResult.genotypes()[i]);
                builder.attribute(GATKSVVCFConstants.DISCORDANT_PAIR_GENOTYPE_QUALITY_ATTRIBUTE, rescaleGq(discordantPairGenotypeResult.genotypeQuals()[i]));
            }
            if (splitReadGenotypeResult != null) {
                builder.attribute(GATKSVVCFConstants.SPLIT_READ_GENOTYPE_ATTRIBUTE, splitReadGenotypeResult.genotypes()[i]);
                builder.attribute(GATKSVVCFConstants.SPLIT_READ_GENOTYPE_QUALITY_ATTRIBUTE, rescaleGq(splitReadGenotypeResult.genotypeQuals()[i]));
            }
            builder.GQ(rescaleGq(combinedGenotypes.get(i).quality));
            final int genotypeInt = combinedGenotypes.get(i).genotype;
            if (genotypeInt > maxGenotype) {
                maxGenotype = genotypeInt;
            }
            final List<Allele> altAlleles = record.getAltAlleles();
            if (altAlleles.size() != 1) {
                throw new UserException.BadInput("Tool only supports biallelic variants but record " + record.getId() + " does not have exactly one alt allele");
            }
            builder.alleles(CanonicalSVCollapser.makeBiallelicList(altAlleles.get(0), record.getRefAllele(), Math.min(genotypeInt, 2), 2));
            final String evidence = String.join(",", combinedGenotypes.get(i).evidence.stream().sorted().map(GATKSVVCFConstants.EvidenceTypes::name).toList());
            builder.attribute(GATKSVVCFConstants.GENOTYPE_EVIDENCE, evidence);
            newGenotypeList.add(builder.make());
        }
        final GenotypesContext newGenotypes = GenotypesContext.create(newGenotypeList);
        final SVCallRecord regenotypedCall = SVCallRecordUtils.copyCallWithNewGenotypes(record, newGenotypes);
        final Map<String, Object> attributes = new HashMap<>(regenotypedCall.getAttributes());
        if (splitReadGenotypeResult != null) {
            attributes.put(GATKSVVCFConstants.BOTHSIDES_SUPPORT_ATTRIBUTE, splitReadGenotypeResult.bothsidePass());
            attributes.put(GATKSVVCFConstants.HIGH_SR_BACKGROUND_ATTRIBUTE, splitReadGenotypeResult.backgroundFail());
        }
        if (combinedVariantQual != null) {
            attributes.put(GATKSVVCFConstants.VARIANT_GENOTYPE_QUALITY_ATTRIBUTE, combinedVariantQual);
        }
        final SVCallRecord newAttributeCall = SVCallRecordUtils.copyCallWithNewAttributes(regenotypedCall, attributes);
        final VariantContextBuilder builder = SVCallRecordUtils.getVariantBuilder(newAttributeCall);
        if (record.isDepthOnly() && record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.DUP && maxGenotype > 2) {
            final List<Allele> newAlleles = new ArrayList<>(2);
            newAlleles.add(record.getRefAllele());
            newAlleles.add(Allele.SV_SIMPLE_CNV);
            builder.alleles(newAlleles);
            builder.attribute(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.StructuralVariantAnnotationType.CNV);
            final List<Genotype> noCallGenotypes = new ArrayList<>(builder.getGenotypes().size());
            for (final Genotype g : builder.getGenotypes()) {
                final GenotypeBuilder gBuilder = new GenotypeBuilder(g);
                gBuilder.alleles(Collections.nCopies(g.getPloidy(), Allele.NO_CALL));
                noCallGenotypes.add(gBuilder.make());
            }
            builder.genotypes(noCallGenotypes);
        }
        writer.add(builder.make());
    }

    private static int getDuplicationState(final Genotype g) {
        return (int) g.getExtendedAttribute(GATKSVVCFConstants.DEPTH_COPY_STATE_ATTRIBUTE) - g.getPloidy();
    }

    private int rescaleGq(final int gq) {
        // Note GATK automatically truncates GQ at 99
        double gqScaleFactor = 99. / (double) maxQual;
        return (int) Math.round(gqScaleFactor * gq);
    }

    @Override
    protected void afterNthPass(final int n) {
    }

    private DepthEvidenceGenotyper.DepthGenotypeResult genotypeDepth(final SVCallRecord record) {
        // Must be a CNV
        final GATKSVVCFConstants.StructuralVariantAnnotationType svtype = record.getType();
        if (svtype != GATKSVVCFConstants.StructuralVariantAnnotationType.DEL && svtype != GATKSVVCFConstants.StructuralVariantAnnotationType.DUP && svtype != GATKSVVCFConstants.StructuralVariantAnnotationType.CNV) {
            return null;
        }
        final DepthMatrix depthMatrix = loader.load(new SimpleInterval(record.getContigA(), record.getPositionA(), record.getPositionB()), sampleMedians);
        // Use depth cutoffs for depth-only variants or any CNV >= 5kbp, matching v1.1 size-based routing.
        // This is consistent with combineGenotypes(), which uses depth as the primary genotype for >= 5kbp.
        final boolean useLargeVariantDepthCutoffs = record.getLength() != null && record.getLength() >= 5000;
        final DepthEvidenceGenotyper genotyper = (record.isDepthOnly() || useLargeVariantDepthCutoffs) ? depthOnlyGenotyper : pesrDepthGenotyper;
        return genotyper.genotype(depthMatrix);
    }

    private DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult genotypeDiscordantPairs(final SVCallRecord record) {
        if (record.isDepthOnly()) {
            return null;
        }
        final List<DiscordantPairEvidence> discordantPairEvidence = discordantPairCollector.collectEvidence(record);
        return discordantPairGenotyper.genotype(record, discordantPairEvidence, discordantPairParameters, masterSampleList);
    }

    private SplitReadEvidenceGenotyper.SplitReadGenotypeResult genotypeSplitReads(final SVCallRecord record) {
        if (record.isDepthOnly()) {
            return null;
        }
        final List<SplitReadEvidence> startSplitReads = splitReadStartCollector.collectEvidence(record);
        final List<SplitReadEvidence> endSplitReads = splitReadEndCollector.collectEvidence(record);
        return splitReadGenotyper.genotype(record, startSplitReads, endSplitReads, splitReadParameters, medianHomIns, medianHomCutoffMultiplier, masterSampleList);
    }

    private boolean discordantPairCollectionEnabled() {
        return discordantPairsFile != null;
    }

    private boolean splitReadCollectionEnabled() {
        return splitReadsFile != null;
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
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.BOTHSIDES_SUPPORT_ATTRIBUTE, 1, VCFHeaderLineType.Flag, "Variant has read-level support for both sides of breakpoint"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.HIGH_SR_BACKGROUND_ATTRIBUTE, 1, VCFHeaderLineType.Flag, "High number of SR splits in background samples indicating messy region"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.VARIANT_GENOTYPE_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Overall variant genotype quality"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.DEPTH_COPY_STATE_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Depth genotype copy state"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.DEPTH_MEDIAN_COPY_RATIO, 1, VCFHeaderLineType.Float, "Median read depth copy ratio"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.DEPTH_GENOTYPE_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Depth genotype quality"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 1, VCFHeaderLineType.Integer, "Expected copy number for ref genotype"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.DISCORDANT_PAIR_GENOTYPE_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Discordant pair genotype"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.DISCORDANT_PAIR_GENOTYPE_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Discordant pair genotype quality"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.SPLIT_READ_GENOTYPE_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read genotype"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.SPLIT_READ_GENOTYPE_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read genotype quality"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.GENOTYPE_EVIDENCE, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Classes of evidence supporting final genotype"));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_QUALITY_KEY));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY));
        return header;
    }
}