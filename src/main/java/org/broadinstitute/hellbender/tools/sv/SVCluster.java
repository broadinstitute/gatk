package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.io.IOUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.SVCallRecordCodec;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * Clusters SVs with similar breakpoints based on coordinates and supporting evidence.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         Unclustered structural variants from
 *     </li>
 *     <li>
 *         PE evidence file
 *     </li>
 *     <li>
 *         SR evidence file
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Structural variant VCF
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk SVCluster
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */

@CommandLineProgramProperties(
        summary = "Clusters structural variants",
        oneLineSummary = "Clusters structural variants",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public final class SVCluster extends GATKTool {
    public static final String SPLIT_READ_LONG_NAME = "split-reads-file";
    public static final String DISCORDANT_PAIRS_LONG_NAME = "discordant-pairs-file";
    public static final String SAMPLE_COVERAGE_LONG_NAME = "sample-coverage";
    public static final String MIN_SIZE_LONG_NAME = "min-size";

    @Argument(
            doc = "Split reads file",
            fullName = SPLIT_READ_LONG_NAME
    )
    private String splitReadsFile;

    @Argument(
            doc = "Discordant pairs file",
            fullName = DISCORDANT_PAIRS_LONG_NAME
    )
    private String discordantPairsFile;

    @Argument(
            doc = "Sample coverage tsv",
            fullName = SAMPLE_COVERAGE_LONG_NAME
    )
    private String sampleCoverageFile;

    @Argument(
            doc = "Input file",
            fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME,
            shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME
    )
    private String inputFile;

    @Argument(
            doc = "Output file",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private String outputFile;

    @Argument(
            doc = "Min event size",
            fullName = MIN_SIZE_LONG_NAME,
            minValue = 0,
            maxValue = Integer.MAX_VALUE,
            optional = true
    )
    private int minEventSize = 50;

    private SAMSequenceDictionary dictionary;

    private final Map<String,IntervalTree> whitelistedIntervalTreeMap = new HashMap<>();
    private FeatureDataSource<SVCallRecord> reader;
    private VariantContextWriter writer;
    private int variantsWritten;

    private FeatureDataSource<SplitReadEvidence> splitReadSource;
    private FeatureDataSource<DiscordantPairEvidence> discordantPairSource;

    private SVDepthOnlyCallDefragmenter defragmenter;
    private List<SVCallRecord> nonDepthRawCallsBuffer;
    private SVClusterEngine clusterEngine;
    private BreakpointRefiner breakpointRefiner;
    private SVEvidenceCollector evidenceCollector;
    private Map<String,Double> sampleCoverageMap;
    private List<String> samplesList;
    private String currentContig;

    public static String END_CONTIG_ATTRIBUTE = "CHR2";
    public static String END_POS_ATTRIBUTE = VCFConstants.END_KEY;
    public static String SVLEN_ATTRIBUTE = GATKSVVCFConstants.SVLEN;
    public static String SVTYPE_ATTRIBUTE = VCFConstants.SVTYPE;
    public static String STRANDS_ATTRIBUTE = "STRANDS";
    public static String START_STRAND_ATTRIBUTE = "SS";
    public static String END_STRAND_ATTRIBUTE = "ES";
    public static String ALGORITHMS_ATTRIBUTE = "ALGORITHMS";
    public static String START_SPLIT_READ_COUNT_ATTRIBUTE = "SSR";
    public static String END_SPLIT_READ_COUNT_ATTRIBUTE = "ESR";
    public static String DISCORDANT_PAIR_COUNT_ATTRIBUTE = "PE";

    public static String DEPTH_ALGORITHM = "depth";

    private final int SPLIT_READ_QUERY_LOOKAHEAD = 0;
    private final int DISCORDANT_PAIR_QUERY_LOOKAHEAD = 0;
    private final int INPUT_QUERY_LOOKAHEAD = 10000;

    @Override
    public void onTraversalStart() {
        dictionary = getBestAvailableSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }
        loadSampleCoverage();
        initializeSplitReadEvidenceDataSource();
        initializeDiscordantPairDataSource();
        loadIntervalTree();

        defragmenter = new SVDepthOnlyCallDefragmenter(dictionary);
        nonDepthRawCallsBuffer = new ArrayList<>();
        clusterEngine = new SVClusterEngine(dictionary);
        breakpointRefiner = new BreakpointRefiner(sampleCoverageMap);
        evidenceCollector = new SVEvidenceCollector(splitReadSource, discordantPairSource, dictionary, progressMeter);

        initializeReader();
        progressMeter.setRecordsBetweenTimeChecks(100);
        writer = createVCFWriter(Paths.get(outputFile));
        writeVCFHeader();
        currentContig = null;
        variantsWritten = 0;
    }

    @Override
    public Object onTraversalSuccess() {
        reader.close();
        writer.close();
        return null;
    }

    private void initializeReader() {
        reader = new FeatureDataSource<>(
                inputFile,
                "inputFile",
                INPUT_QUERY_LOOKAHEAD,
                SVCallRecord.class,
                getDefaultCloudPrefetchBufferSize(),
                getDefaultCloudIndexPrefetchBufferSize());
        reader.setIntervalsForTraversal(getTraversalIntervals());
    }

    private void initializeSplitReadEvidenceDataSource() {
        splitReadSource = new FeatureDataSource<>(
                splitReadsFile,
                "splitReadsFile",
                SPLIT_READ_QUERY_LOOKAHEAD,
                SplitReadEvidence.class,
                cloudPrefetchBuffer,
                cloudIndexPrefetchBuffer);
    }

    private void initializeDiscordantPairDataSource() {
        discordantPairSource = new FeatureDataSource<>(
                discordantPairsFile,
                "discordantPairsFile",
                DISCORDANT_PAIR_QUERY_LOOKAHEAD,
                DiscordantPairEvidence.class,
                cloudPrefetchBuffer,
                cloudIndexPrefetchBuffer);
    }

    private void loadSampleCoverage() {
        try {
            sampleCoverageMap = IOUtils.readLines(BucketUtils.openFile(sampleCoverageFile), Charset.defaultCharset()).stream()
                    .map(line -> line.split("\t"))
                    .collect(Collectors.toMap(tokens -> tokens[0], tokens -> Double.valueOf(tokens[1])));
            samplesList = IOUtils.readLines(BucketUtils.openFile(sampleCoverageFile), Charset.defaultCharset()).stream()
                    .map(line -> line.split("\t")[0])
                    .collect(Collectors.toList());
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(sampleCoverageFile, e);
        }
    }

    private void loadIntervalTree() {
        for (final SimpleInterval interval : getTraversalIntervals()) {
            whitelistedIntervalTreeMap.putIfAbsent(interval.getContig(), new IntervalTree());
            whitelistedIntervalTreeMap.get(interval.getContig()).put(interval.getStart(), interval.getEnd(), null);
        }
    }

    @Override
    public void traverse() {
        StreamSupport.stream(Spliterators.spliteratorUnknownSize(reader.iterator(), Spliterator.ORDERED), false)
                .filter(this::isValidSize)
                .filter(this::isWhitelisted)
                .forEachOrdered(this::processRecord);
        if (!defragmenter.isEmpty()) {
            processClusters();
        }
    }

    private void processRecord(final SVCallRecord record) {
        if (!record.getContig().equals(currentContig)) {
            if (currentContig != null) {
                processClusters();
            }
            currentContig = record.getContig();
        }
        if (SVDepthOnlyCallDefragmenter.isDepthOnlyCall(record)) {
            defragmenter.add(new SVCallRecordWithEvidence(record));
        } else {
            nonDepthRawCallsBuffer.add(record);
        }
    }

    private void processClusters() {
        final List<SVCallRecordWithEvidence> defragmentedCalls = defragmenter.getOutput();
        defragmentedCalls.stream().forEachOrdered(clusterEngine::add);
        nonDepthRawCallsBuffer.stream().map(SVCallRecordWithEvidence::new).forEachOrdered(clusterEngine::add);
        nonDepthRawCallsBuffer.clear();
        final List<SVCallRecordWithEvidence> clusteredCalls = clusterEngine.getOutput();
        final List<SVCallRecordWithEvidence> callsWithEvidence = evidenceCollector.collectEvidence(clusteredCalls);
        final List<SVCallRecordWithEvidence> refinedCalls = callsWithEvidence.stream().map(breakpointRefiner::refineCall).collect(Collectors.toList());
        final List<SVCallRecordWithEvidence> finalCalls = clusterEngine.deduplicateItems(refinedCalls);
        write(finalCalls);
    }

    private void write(final List<SVCallRecordWithEvidence> calls) {
        calls.stream()
                .sorted(Comparator.comparing(c -> c.getStartAsInterval(), IntervalUtils.getDictionaryOrderComparator(dictionary)))
                .map(this::buildVariantContext)
                .forEachOrdered(writer::add);
    }

    private void writeVCFHeader() {
        final VCFHeader header = new VCFHeader(getDefaultToolVCFHeaderLines(), samplesList);
        header.setSequenceDictionary(dictionary);
        header.addMetaDataLine(new VCFInfoHeaderLine(END_CONTIG_ATTRIBUTE, 1, VCFHeaderLineType.String, "End contig"));
        header.addMetaDataLine(new VCFInfoHeaderLine(END_POS_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "End position"));
        header.addMetaDataLine(new VCFInfoHeaderLine(SVLEN_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Variant length"));
        header.addMetaDataLine(new VCFInfoHeaderLine(SVTYPE_ATTRIBUTE, 1, VCFHeaderLineType.String, "Variant type"));
        header.addMetaDataLine(new VCFInfoHeaderLine(START_STRAND_ATTRIBUTE, 1, VCFHeaderLineType.String, "Start strand"));
        header.addMetaDataLine(new VCFInfoHeaderLine(END_STRAND_ATTRIBUTE, 1, VCFHeaderLineType.String, "End strand"));
        header.addMetaDataLine(new VCFInfoHeaderLine(ALGORITHMS_ATTRIBUTE, 1, VCFHeaderLineType.String, "List of calling algorithms"));
        header.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1, VCFHeaderLineType.String, "Genotype"));
        header.addMetaDataLine(new VCFFormatHeaderLine(START_SPLIT_READ_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read count at start of variant"));
        header.addMetaDataLine(new VCFFormatHeaderLine(END_SPLIT_READ_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read count at end of variant"));
        header.addMetaDataLine(new VCFFormatHeaderLine(DISCORDANT_PAIR_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Discordant pair count"));
        writer.writeHeader(header);
    }

    public VariantContext buildVariantContext(final SVCallRecordWithEvidence call) {
        Utils.nonNull(call);
        final Allele altAllele = Allele.create("<" + call.getType().name() + ">", false);
        final Allele refAllele = Allele.REF_N;
        final VariantContextBuilder builder = new VariantContextBuilder("", call.getContig(), call.getStart(), call.getEnd(),
                Lists.newArrayList(refAllele, altAllele));
        builder.attribute(END_CONTIG_ATTRIBUTE, call.getEndContig());
        builder.attribute(END_POS_ATTRIBUTE, call.getEnd());
        builder.attribute(SVLEN_ATTRIBUTE, call.getLength());
        builder.attribute(SVTYPE_ATTRIBUTE, call.getType());
        builder.attribute(START_STRAND_ATTRIBUTE, getStrandString(call.getStartStrand()));
        builder.attribute(END_STRAND_ATTRIBUTE, getStrandString(call.getEndStrand()));
        builder.attribute(ALGORITHMS_ATTRIBUTE, call.getAlgorithms());
        final List<Genotype> genotypes = new ArrayList<>();
        final Map<String,Integer> startSplitReadCounts = getSplitReadCountsAtPosition(call.getStartSplitReadSites(), call.getStart());
        final Map<String,Integer> endSplitReadCounts = getSplitReadCountsAtPosition(call.getEndSplitReadSites(), call.getEnd());
        final Map<String,Integer> discordantPairCounts = getDiscordantPairCountsMap(call.getDiscordantPairs());
        for (final String sample : sampleCoverageMap.keySet()) {
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sample);
            genotypeBuilder.attribute(START_SPLIT_READ_COUNT_ATTRIBUTE, startSplitReadCounts.getOrDefault(sample, 0));
            genotypeBuilder.attribute(END_SPLIT_READ_COUNT_ATTRIBUTE, endSplitReadCounts.getOrDefault(sample, 0));
            genotypeBuilder.attribute(DISCORDANT_PAIR_COUNT_ATTRIBUTE, discordantPairCounts.getOrDefault(sample, 0));
            if (call.getSamples().contains(sample)) {
                genotypeBuilder.alleles(Lists.newArrayList(refAllele, altAllele));
            } else {
                genotypeBuilder.alleles(Lists.newArrayList(refAllele, refAllele));
            }
            genotypes.add(genotypeBuilder.make());
        }
        builder.genotypes(genotypes);
        builder.id(String.format("SVx%08X", variantsWritten));
        variantsWritten++;
        return builder.make();
    }

    private static String getStrandString(final boolean strand) {
        return strand ? SVCallRecordCodec.STRAND_PLUS : SVCallRecordCodec.STRAND_MINUS;
    }

    private static Map<String,Integer> getSplitReadCountsAtPosition(final List<SplitReadSite> sites, final int pos) {
        Utils.nonNull(sites);
        Utils.validateArg(pos > 0, "Non-positive position");
        if (sites.stream().map(SplitReadSite::getPosition).distinct().count() != sites.size()) {
            throw new IllegalArgumentException("Sites did not have unique positions");
        }
        return sites.stream()
                .filter(s -> s.getPosition() == pos)
                .map(SplitReadSite::getSampleCountsMap)
                .findAny()
                .orElse(Collections.emptyMap());
    }

    private Map<String,Integer> getDiscordantPairCountsMap(final Collection<DiscordantPairEvidence> discordantPairs) {
        Utils.nonNull(discordantPairs);
        return discordantPairs.stream()
                .collect(Collectors.groupingBy(DiscordantPairEvidence::getSample,
                        Collectors.reducing(0, e -> 1, Integer::sum)));
    }

    private boolean isValidSize(final SVCallRecord call) {
        return !((call.getType().equals(StructuralVariantType.DEL)
                || call.getType().equals(StructuralVariantType.DUP)
                || call.getType().equals(StructuralVariantType.INV))
                && call.getLength() < minEventSize);
    }

    private boolean isWhitelisted(final SVCallRecord call) {
        final IntervalTree startTree = whitelistedIntervalTreeMap.get(call.getContig());
        if (startTree == null || !startTree.overlappers(call.getStart(), call.getStart() + 1).hasNext()) {
            return false;
        }
        final IntervalTree endTree = whitelistedIntervalTreeMap.get(call.getEndContig());
        if (endTree == null || !endTree.overlappers(call.getEnd(), call.getEnd() + 1).hasNext()) {
            return false;
        }
        return true;
    }
}
