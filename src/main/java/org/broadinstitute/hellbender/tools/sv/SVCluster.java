package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
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
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Paths;
import java.util.*;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

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
@ExperimentalFeature
@DocumentedFeature
public final class SVCluster extends VariantWalker {
    public static final String SPLIT_READ_LONG_NAME = "split-reads-file";
    public static final String DISCORDANT_PAIRS_LONG_NAME = "discordant-pairs-file";
    public static final String SAMPLE_COVERAGE_LONG_NAME = "sample-coverage";
    public static final String MIN_SIZE_LONG_NAME = "min-size";
    public static final String DEPTH_ONLY_INCLUDE_INTERVAL_OVERLAP_LONG_NAME = "depth-include-overlap";
    public static final String MIN_DEPTH_SIZE_LONG_NAME = "min-depth-size";
    public static final String VARIANT_PREFIX_LONG_NAME = "variant-prefix";

    @Argument(
            doc = "Split reads file",
            fullName = SPLIT_READ_LONG_NAME
    )
    private GATKPath splitReadsFile;

    @Argument(
            doc = "Discordant pairs file",
            fullName = DISCORDANT_PAIRS_LONG_NAME
    )
    private GATKPath discordantPairsFile;

    @Argument(
            doc = "Sample coverage tsv",
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
            doc = "Min event size",
            fullName = MIN_SIZE_LONG_NAME,
            minValue = 0,
            maxValue = Integer.MAX_VALUE,
            optional = true
    )
    private int minEventSize = 50;

    @Argument(
            doc = "Depth-only call min included intervals overlap",
            fullName = DEPTH_ONLY_INCLUDE_INTERVAL_OVERLAP_LONG_NAME,
            minValue = 0,
            maxValue = 1,
            optional = true
    )
    private double minDepthOnlyIncludeOverlap = 0.5;

    @Argument(
            doc = "Minimum depth-only call size to emit",
            fullName = MIN_DEPTH_SIZE_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private int minDepthOnlySize = 5000;

    @Argument(
            doc = "Prefix for variant IDs",
            fullName = VARIANT_PREFIX_LONG_NAME,
            optional = true
    )
    private String variantPrefix = "SV_x";

    private SAMSequenceDictionary dictionary;
    private final Map<String,IntervalTree<Object>> includedIntervalsTreeMap = new HashMap<>();
    private VariantContextWriter writer;
    private FeatureDataSource<SplitReadEvidence> splitReadSource;
    private FeatureDataSource<DiscordantPairEvidence> discordantPairSource;
    private SVDepthOnlyCallDefragmenter defragmenter;
    private List<SVCallRecord> nonDepthRawCallsBuffer;
    private SVClusterEngine clusterEngine;
    private BreakpointRefiner breakpointRefiner;
    private PairedEndAndSplitReadEvidenceAggregator evidenceCollector;
    private SVCallRecordDeduplicator<SVCallRecordWithEvidence> deduplicator;
    private Map<String,Double> sampleCoverageMap;
    private Set<String> samples;
    private String currentContig;
    private int numVariantsWritten = 0;

    private final int SPLIT_READ_QUERY_LOOKAHEAD = 0;
    private final int DISCORDANT_PAIR_QUERY_LOOKAHEAD = 0;

    @Override
    public boolean ignoresIntervalsForTraversal() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        dictionary = getBestAvailableSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }
        samples = new LinkedHashSet<>(getHeaderForVariants().getSampleNamesInOrder());
        loadSampleCoverage();
        initializeSplitReadEvidenceDataSource();
        initializeDiscordantPairDataSource();
        loadIntervalTree();

        defragmenter = new SVDepthOnlyCallDefragmenter(dictionary);
        nonDepthRawCallsBuffer = new ArrayList<>();
        clusterEngine = new SVClusterEngineNoCNV(dictionary, false, SVClusterEngine.BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END);
        breakpointRefiner = new BreakpointRefiner(sampleCoverageMap, dictionary);
        evidenceCollector = new PairedEndAndSplitReadEvidenceAggregator(splitReadSource, discordantPairSource, dictionary, null);

        final Function<Collection<SVCallRecordWithEvidence>,SVCallRecordWithEvidence> collapser = SVCallRecordUtils::deduplicateWithRawCallAttributeWithEvidence;
        deduplicator = new SVCallRecordDeduplicator<>(collapser, dictionary);

        writer = createVCFWriter(Paths.get(outputFile));
        writeVCFHeader();
        currentContig = null;
    }

    @Override
    public Object onTraversalSuccess() {
        if (!defragmenter.isEmpty() || !nonDepthRawCallsBuffer.isEmpty()) {
            processClusters();
        }
        writer.close();
        return null;
    }

    private void initializeSplitReadEvidenceDataSource() {
        splitReadSource = new FeatureDataSource<>(
                splitReadsFile.toString(),
                "splitReadsFile",
                SPLIT_READ_QUERY_LOOKAHEAD,
                SplitReadEvidence.class,
                cloudPrefetchBuffer,
                cloudIndexPrefetchBuffer);
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

    private void loadIntervalTree() {
        final List<SimpleInterval> intervals = getRequestedIntervals();
        if (intervals == null) {
            throw new UserException.MissingReference("Reference dictionary is required");
        }
        for (final SimpleInterval interval : intervals) {
            includedIntervalsTreeMap.putIfAbsent(interval.getContig(), new IntervalTree<>());
            includedIntervalsTreeMap.get(interval.getContig()).put(interval.getStart(), interval.getEnd(), null);
        }
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext,
                      final ReferenceContext referenceContext, final FeatureContext featureContext) {
        final Predicate<Genotype> genotypeFilter = g -> SVCallRecord.isRawCall(g);
        final Function<Genotype, Map<String,Object>> attributeGenerator = g -> Collections.emptyMap();
        final SVCallRecord originalCall = SVCallRecordUtils.create(variant);
        final GenotypesContext filteredGenotypes = SVCallRecordUtils.filterAndAddGenotypeAttributes(originalCall.getGenotypes(), genotypeFilter, attributeGenerator, false);
        final List<Allele> refAlleles = Arrays.asList(Allele.REF_N, Allele.REF_N);
        final List<Allele> nonRefAlleles = Arrays.asList(Allele.REF_N, Allele.NON_REF_ALLELE);
        final Predicate<Genotype> nonRefPredicate = g -> SVCallRecord.isRawCall(g);
        final GenotypesContext genotypes = SVCallRecordUtils.predicateGenotypeAlleles(filteredGenotypes, nonRefPredicate, nonRefAlleles, refAlleles);
        final SVCallRecord call = SVCallRecordUtils.copyCallWithNewGenotypes(originalCall, genotypes);

        // Filter
        if (!SVCallRecordUtils.isValidSize(call, minEventSize) || !SVCallRecordUtils.intervalIsIncluded(call, includedIntervalsTreeMap, minDepthOnlyIncludeOverlap)) {
            return;
        }

        // Flush clusters if we hit the next contig
        if (!call.getContigA().equals(currentContig)) {
            if (currentContig != null) {
                processClusters();
            }
            currentContig = call.getContigA();
        }

        // Add to clustering buffers
        if (SVDepthOnlyCallDefragmenter.isDepthOnlyCall(call)) {
            defragmenter.add(call);
        } else {
            nonDepthRawCallsBuffer.add(call);
        }
    }

    private void processClusters() {
        logger.info("Processing contig " + currentContig);
        final Stream<SVCallRecord> defragmentedStream = defragmenter.getOutput().stream();
        final Stream<SVCallRecord> nonDepthStream = nonDepthRawCallsBuffer.stream()
                .flatMap(SVCallRecordUtils::convertInversionsToBreakends);
        //Combine and sort depth and non-depth calls because they must be added in dictionary order
        Stream.concat(defragmentedStream, nonDepthStream)
                .sorted(SVCallRecordUtils.getCallComparator(dictionary))
                .forEachOrdered(clusterEngine::add);
        nonDepthRawCallsBuffer.clear();
        logger.info("Clustering...");
        final List<SVCallRecord> clusteredCalls = clusterEngine.getOutput();
        logger.info("Aggregating evidence...");
        final List<SVCallRecordWithEvidence> callsWithEvidence = evidenceCollector.collectEvidence(clusteredCalls);
        logger.info("Filtering and refining breakpoints...");
        final List<SVCallRecordWithEvidence> refinedCalls = callsWithEvidence.stream()
                .filter(call -> !SVDepthOnlyCallDefragmenter.isDepthOnlyCall(call) || call.getLength() >= minDepthOnlySize)
                .map(breakpointRefiner::refineCall)
                .sorted(SVCallRecordUtils.getCallComparator(dictionary))
                .collect(Collectors.toList());
        logger.info("Deduplicating variants...");
        final List<SVCallRecordWithEvidence> finalCalls = deduplicator.deduplicateItems(refinedCalls);
        logger.info("Writing to file...");
        write(finalCalls);
        logger.info("Contig " + currentContig + " successfully processed");
    }

    private void write(final List<SVCallRecordWithEvidence> calls) {
        calls.stream()
                .sorted(SVCallRecordUtils.getCallComparator(dictionary))
                .map(this::buildVariantContext)
                .forEachOrdered(writer::add);
    }

    private void writeVCFHeader() {
        final VCFHeader header = new VCFHeader(getDefaultToolVCFHeaderLines(), samples);
        header.setSequenceDictionary(dictionary);
        for (final VCFHeaderLine line : getHeaderForVariants().getMetaDataInInputOrder()) {
            header.addMetaDataLine(line);
        }
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.START_SPLIT_READ_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read count at start of variant"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.END_SPLIT_READ_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read count at end of variant"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.DISCORDANT_PAIR_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Discordant pair count"));
        writer.writeHeader(header);
    }

    public VariantContext buildVariantContext(final SVCallRecordWithEvidence call) {
        final GenotypesContext filledGenotypes = SVCallRecordUtils.fillMissingSamplesWithEmptyGenotypes(call.getGenotypes(), samples);

        final Function<Genotype, Map<String,Object>> attributeGenerator =
                g -> Collections.singletonMap(GATKSVVCFConstants.RAW_CALL_ATTRIBUTE, g.hasAnyAttribute(GATKSVVCFConstants.RAW_CALL_ATTRIBUTE) ? GATKSVVCFConstants.RAW_CALL_ATTRIBUTE_TRUE : GATKSVVCFConstants.RAW_CALL_ATTRIBUTE_FALSE);
        final GenotypesContext rawCallSetGenotypes = SVCallRecordUtils.filterAndAddGenotypeAttributes(filledGenotypes, g -> true, attributeGenerator, false);

        final GenotypesContext nonCallGenotypes = SVCallRecordUtils.predicateGenotypeAlleles(rawCallSetGenotypes, g -> true, Collections.emptyList(), Collections.emptyList());
        final String newId = String.format("%s%08x", variantPrefix, numVariantsWritten++);
        final SVCallRecordWithEvidence finalCall = new SVCallRecordWithEvidence(newId, call.getContigA(), call.getPositionA(), call.getStrandA(), call.getContigB(),
                call.getPositionB(), call.getStrandB(), call.getType(), call.getLength(), call.getAlgorithms(),
                nonCallGenotypes, call.getStartSplitReadSites(), call.getEndSplitReadSites(), call.getDiscordantPairs(),
                call.getCopyNumberDistribution());
        return SVCallRecordUtils.createBuilderWithEvidence(SVCallRecordUtils.copyCallWithNewGenotypes(finalCall, nonCallGenotypes)).make();
    }

}
