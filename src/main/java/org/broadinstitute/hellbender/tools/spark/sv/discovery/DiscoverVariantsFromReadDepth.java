package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledContigPloidyCollection;
import org.broadinstitute.hellbender.tools.spark.sv.DiscoverVariantsFromReadDepthArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth.*;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import scala.Tuple2;

import java.io.*;
import java.nio.file.Paths;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * Call structural variants using SV evidence and germline CNV (gCNV) copy number calls
 *
 * <p>See StructuralVariationDiscoveryPipelineSpark and GermlineCNVCaller for information on the prerequisite pipelines.</p>
 *
 * <p>This tool is currently designed for detecting large (>5kbp) intrachromosomal events.</p>
 *
 * <h3>Methods</h3>
 *
 * <p>The tool looks for events supported by a combination of discordant read pairs (RP), split reads (SR), SV pipeline calls (SVPCs),
 * assembled breakpoint pairs (BNDs), and copy number interval calls (a posterior likelihoods).</p>
 *
 * <p>Events are resolved using a graph-based approach that integrates all of the evidence types. A sequence graph is assembled
 * from RP, SR, SVPC, and BND evidence. All possible paths on the graph are systematically enumerated, and the probability of each
 * path is calculated based on copy number posteriors. Calls are made for each event type (DEL, DUP, INV, DUP_INV) by summing
 * the probability over all paths at each reference position. Events with probability above a minimum threshold are then called.
 * This approach is able to enumerate all possible haplotype combinations and resolve simple events comprising complex ones.</p>
 *
 * <p>Note that enumerating all paths across the entire genome would be computationally intractable. Therefore, the method employs a divide-and-conquer
 * strategy, in which the full graph is partitioned into subgraphs that likely represent independent events. This permits path
 * enumeration over most of the genome. In some partitions, however, local breakpoint complexity still results in combinatorial
 * explosion of path branching. These cases are reported in the output as unresolved events (type U).</p>
 *
 * <h3>Input</h3>
 *
 * <ul>
     * <li>Structural variant call VCF</li>
     * <li>Breakpoint (BND) VCF</li>
     * <li>Germline CNV interval calls VCF</li>
     * <li>Evidence target links .bedpe file from StructuralVariationDiscoveryPipelineSpark</li>
     * <li>High-coverage intervals from StructuralVariationDiscoveryPipelineSpark (optional)</li>
     * <li>Blacklisted intervals file (optional)</li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
     * <li>Event calls file</li>
     * <li>Unresolved event intervals file/li>
     * <li>Event haplotypes file</li>
 * </ul>
 *
 * <p>The tool may also optionally produce:</p>
 *
 * <ul>
 * <li>Breakpoint graph file</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 * gatk DiscoverVariantsFromReadDepth  \
 *   --sv-calls-vcf calls.vcf \
 *   --bnd-vcf breakpoints.vcf \
 *   --evidence-target-links-file target_links.bedpe \
 *   --cnv-intervals-file genotyped-intervals.vcf \
 *   --output output/
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */

@ExperimentalFeature
@CommandLineProgramProperties(
        oneLineSummary = "Discovers structural variants using SV evidence and read depth",
        summary = "This tool uses the output from StructuralVariationDiscoveryPipelineSpark and the gCNV pipeline " +
                "to call structural variation events.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public final class DiscoverVariantsFromReadDepth extends CommandLineProgram {

    public static final String SV_CALLS_LONG_NAME = "sv-calls-vcf";
    public static final String BREAKPOINTS_LONG_NAME = "bnd-vcf";
    public static final String EVIDENCE_TARGET_LINKS_LONG_NAME = "evidence-target-links-file";
    public static final String CNV_INTERVALS_LONG_NAME = "cnv-intervals-file";
    public static final String HIGH_COVERAGE_INTERVALS_LONG_NAME = "high-coverage-intervals";
    public static final String BLACKLIST_LONG_NAME = "blacklist";
    public static final String GRAPH_OUTPUT_LONG_NAME = "write-graph";
    public static final String CONTIG_PLOIDY_CALLS_LONG_NAME = "ploidy-calls";
    public static final String STRUCTURAL_VARIANT_TABLE_LONG_NAME = "table";
    public static final String MAPPABILITY_INTERVALS_LONG_NAME = "mappable-intervals";

    public static final String RESOLVED_CALLS_FILE_NAME = "calls";
    public static final String UNRESOLVED_CALLS_FILE_NAME = "unresolved";
    public static final String HAPLOTYPES_FILE_NAME = "haplotypes";
    public static final String GRAPH_FILE_NAME = "graph";

    @ArgumentCollection
    private final DiscoverVariantsFromReadDepthArgumentCollection arguments = new DiscoverVariantsFromReadDepthArgumentCollection();

    @Argument(doc = "Sequence dictionary", fullName = "sequence-dictionary")
    private String sequenceDictionaryPath;
    @Argument(doc = "BND list (.vcf)", fullName = BREAKPOINTS_LONG_NAME, optional = true)
    private String breakpointVCFPath;
    @Argument(doc = "Structural variant calls list (.vcf)", fullName = SV_CALLS_LONG_NAME, optional = true)
    private String svCallVCFPath;
    @Argument(doc = "Evidence target links file (.bedpe)", fullName = EVIDENCE_TARGET_LINKS_LONG_NAME, optional = true)
    private String evidenceTargetLinksFilePath;
    @Argument(doc = "Germline CNV (gCNV) interval call vcfs (.vcf), can specify more than one", fullName = CNV_INTERVALS_LONG_NAME)
    private List<String> intervalCallsFilePaths;
    @Argument(doc = "High coverage intervals file path", fullName = HIGH_COVERAGE_INTERVALS_LONG_NAME, optional = true)
    private String highCoverageIntervalsPath;
    @Argument(doc = "One or more blacklisted intervals files", fullName = BLACKLIST_LONG_NAME, optional = true)
    private List<String> blacklistedPaths;
    @Argument(doc = "Output directory path (must exist)", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputDirectoryPath;
    @Argument(doc = "Sample name", shortName = StandardArgumentDefinitions.SAMPLE_NAME_SHORT_NAME, fullName = StandardArgumentDefinitions.SAMPLE_NAME_LONG_NAME)
    private String sampleName;
    @Argument(doc = "Write graph output", fullName = GRAPH_OUTPUT_LONG_NAME, optional = true)
    private boolean writeGraph = false;

    @Argument(
            doc = "Structural variants table",
            fullName = STRUCTURAL_VARIANT_TABLE_LONG_NAME,
            optional = true
    )
    private File structuralVariantTableFile;

    @Argument(
            doc = "Contig ploidy file",
            fullName = CONTIG_PLOIDY_CALLS_LONG_NAME
    )
    private File contigPloidyCallsFile;

    @Argument(
            doc = "Mappable intervals",
            fullName = MAPPABILITY_INTERVALS_LONG_NAME
    )
    private File mappabilityIntervalsFile;

    @Override
    public Object doWork() {

        //Read input data
        final SAMSequenceDictionary dictionary = ReferenceUtils.loadFastaDictionary(BucketUtils.openFile(sequenceDictionaryPath));
        final List<SVInterval> highCoverageIntervals = getIntervals(highCoverageIntervalsPath, dictionary);
        final List<SVInterval> blacklist = new ArrayList<>();
        for (final String path : blacklistedPaths) {
            blacklist.addAll(getIntervals(path, dictionary));
        }
        final Collection<VariantContext> breakpointCalls = breakpointVCFPath == null? Collections.emptyList() : readVCF(breakpointVCFPath, dictionary);
        final Collection<VariantContext> svCalls = svCallVCFPath == null ? Collections.emptyList() : readVCF(svCallVCFPath, dictionary);
        final List<SVIntervalTree<SVCopyNumberInterval>> copyNumberTrees = new ArrayList<>(intervalCallsFilePaths.size());
        for (final String path : intervalCallsFilePaths) {
            final List<SVCopyNumberInterval> copyNumberIntervals = getCalledCopyNumberIntervals(path, dictionary);
            copyNumberTrees.add(SVIntervalUtils.buildCopyNumberIntervalTree(copyNumberIntervals));
        }
        final SVMultiscaleIntervalTree<SVCopyNumberInterval> multiscaleIntervalTree = new SVMultiscaleIntervalTree<>(copyNumberTrees);
        final Collection<EvidenceTargetLink> evidenceTargetLinks = evidenceTargetLinksFilePath == null ? Collections.emptyList() : getEvidenceTargetLinks(evidenceTargetLinksFilePath, dictionary);
        final Collection<EventRecord> eventRecordCollection = structuralVariantTableFile == null ? Collections.emptyList() : EventRecordReader.readEventRecords(structuralVariantTableFile).collect(Collectors.toList());
        final SVIntervalTree<Object> mappabilityIntervals = readMappabilityIntervals(mappabilityIntervalsFile, dictionary);

        //Create graph
        final SVEvidenceIntegrator evidenceIntegrator = new SVEvidenceIntegrator(breakpointCalls, svCalls, evidenceTargetLinks, eventRecordCollection, multiscaleIntervalTree, highCoverageIntervals, blacklist, dictionary, arguments);
        final SVGraph graph = evidenceIntegrator.buildGraph();

        //Call events
        final ReadDepthSVCaller caller = new ReadDepthSVCaller(graph, multiscaleIntervalTree, arguments);
        final Tuple2<Collection<CalledSVGraphGenotype>, Collection<CalledSVGraphEvent>> callResult = caller.callEvents();
        final Collection<CalledSVGraphGenotype> haplotypes = callResult._1;
        final Collection<CalledSVGraphEvent> events = callResult._2;

        //Write output
        final Collection<CalledSVGraphEvent> resolvedEvents = events.stream().filter(CalledSVGraphEvent::isResolved).collect(Collectors.toList());
        resolvedEvents.stream().forEach(event -> applyMappability(event, mappabilityIntervals));
        markDuplicateEvents(resolvedEvents);
        final Collection<CalledSVGraphEvent> unresolvedEvents = events.stream().filter(event -> !event.isResolved()).collect(Collectors.toList());
        writeEvents(resolvedEvents, RESOLVED_CALLS_FILE_NAME, dictionary);
        writeEvents(unresolvedEvents, UNRESOLVED_CALLS_FILE_NAME, dictionary);
        writeHaplotypes(haplotypes);
        if (writeGraph) {
            writeGraph(graph);
        }
        return null;
    }

    private static void applyMappability(final CalledSVGraphEvent event, final SVIntervalTree<Object> mappabilityIntervals) {
        final SVInterval eventInterval = event.getInterval();
        final double MIN_MAPPABILITY_OVERLAP = 0.50;
        final int overlap = iteratorToStream(mappabilityIntervals.overlappers(event.getInterval()))
                .mapToInt(entry -> entry.getInterval().overlapLen(eventInterval)).sum();
        if (overlap < MIN_MAPPABILITY_OVERLAP * eventInterval.getLength()) {
            event.setLowMappability(true);
        } else {
            event.setLowMappability(false);
        }
    }

    // TODO sensitive to input order in the case of ties
    private static void markDuplicateEvents(final Collection<CalledSVGraphEvent> events) {
        final double DUPLICATE_RECIPROCAL_OVERLAP = 0.9;
        final SVIntervalTree<CalledSVGraphEvent> tree = new SVIntervalTree<>();
        for (final CalledSVGraphEvent event : events) {
            tree.put(event.getInterval(), event);
        }
        final Set<CalledSVGraphEvent> testedEvents = new HashSet<>(SVUtils.hashMapCapacity(events.size()));
        for (final CalledSVGraphEvent event : events) {
            if (testedEvents.contains(event)) continue;
            final Collection<CalledSVGraphEvent> overlappers = iteratorToStream(tree.overlappers(event.getInterval()))
                    .map(SVIntervalTree.Entry::getValue)
                    .filter(e -> SVIntervalUtils.hasReciprocalOverlap(event.getInterval(), e.getInterval(), DUPLICATE_RECIPROCAL_OVERLAP))
                    .collect(Collectors.toList());
            final double maxQual = overlappers.stream().mapToDouble(CalledSVGraphEvent::getRefQuality).max().getAsDouble();
            boolean foundMax = false;
            for (final CalledSVGraphEvent overlapper : overlappers) {
                if (overlapper.getRefQuality() < maxQual || foundMax) {
                    overlapper.setDuplicate(true);
                }
            }
            testedEvents.addAll(overlappers);
        }
    }

    private static <T> Stream<T> iteratorToStream(final Iterator<T> iter) {
        return StreamSupport.stream(Spliterators.spliteratorUnknownSize(iter, Spliterator.ORDERED), false);
    }

    private static SVIntervalTree<Object> readMappabilityIntervals(final File file, final SAMSequenceDictionary dictionary) {
        final List<SVInterval> intervals = getIntervals(file.getAbsolutePath(), dictionary);
        final SVIntervalTree<Object> tree = new SVIntervalTree<>();
        for (final SVInterval interval : intervals) {
            tree.put(interval, null);
        }
        return tree;
    }

    /**
     * Reads VCF to collection of VariantContext
     */
    private static Collection<VariantContext> readVCF(final String vcfPath, final SAMSequenceDictionary dictionary) {
        try (final FeatureDataSource<VariantContext> dataSource =
                     new FeatureDataSource<>(vcfPath, null, 0, VariantContext.class)) {
            Collection<VariantContext> variants = Utils.stream(dataSource.iterator()).collect(Collectors.toList());
            if (variants.stream().anyMatch(variant -> dictionary.getSequenceIndex(variant.getContig()) < 0)) {
                throw new UserException("Found a VCF entry contig that does not appear in the dictionary.");
            }
            return variants;
        }
    }

    /**
     * Reads intervals file to SVInterval list
     */
    private static List<SVInterval> getIntervals(final String path, final SAMSequenceDictionary dictionary) {
        if (path == null) {
            return Collections.emptyList();
        }
        final GenomeLocParser genomeLocParser = new GenomeLocParser(dictionary);
        final List<GenomeLoc> intervals = IntervalUtils.parseIntervalArguments(genomeLocParser, path);
        final GenomeLocSortedSet mergedIntervals = IntervalUtils.sortAndMergeIntervals(genomeLocParser, intervals, IntervalMergingRule.ALL);
        return mergedIntervals.toList().stream()
                .map(genomeLoc -> new SVInterval(genomeLoc.getContigIndex(), genomeLoc.getStart(), genomeLoc.getEnd()))
                .collect(Collectors.toList());
    }

    /**
     * Reads CNV interval calls to list of SVCopyNumberInterval
     */
    private static List<SVCopyNumberInterval> getCalledCopyNumberIntervals(final String path, final SAMSequenceDictionary dictionary) {
        final File file = new File(path);
        if (path.toLowerCase().endsWith(FileExtensions.VCF.toLowerCase()) || path.toLowerCase().endsWith(FileExtensions.COMPRESSED_VCF.toLowerCase())) {
            final VCFFileReader reader = new VCFFileReader(file, false);
            return Utils.stream(reader.iterator()).map(variantContext -> new SVCopyNumberInterval(variantContext, dictionary)).collect(Collectors.toList());
        }
        throw new UserException.BadInput("CNV intervals file must be " + FileExtensions.VCF + " or " + FileExtensions.COMPRESSED_VCF);
    }

    /**
     * Reads evidence target links file
     */
    private static Collection<EvidenceTargetLink> getEvidenceTargetLinks(final String path, final SAMSequenceDictionary dictionary) {
        final File file = new File(path);
        final Collection<EvidenceTargetLink> links = new ArrayList<>();
        final String[] evidenceTargetFileLines = new String(IOUtils.readFileIntoByteArray(file)).split("\n");
        for (final String line : evidenceTargetFileLines) {
            links.add(EvidenceTargetLink.fromBedpeString(line.trim(), dictionary));
        }
        return links;
    }

    /**
     * Writes events collection to BED file
     */
    private void writeEvents(final Collection<CalledSVGraphEvent> events, final String name, final SAMSequenceDictionary dictionary) {
        final String fileName = sampleName + "-" + name + ".bed";
        final String filePath = Paths.get(outputDirectoryPath, fileName).toAbsolutePath().toString();
        try (final OutputStream outputStream = BucketUtils.createFile(filePath)) {
            outputStream.write(("#" + CalledSVGraphEvent.bedHeader() + "\n").getBytes());
            for (final CalledSVGraphEvent event : events) {
                outputStream.write((event.bedString(dictionary) + "\n").getBytes());
            }
        } catch (final IOException e) {
            throw new GATKException("Error writing calls BED file", e);
        }
    }

    /**
     * Writes haplotypes collection to BED file
     */
    private void writeHaplotypes(final Collection<CalledSVGraphGenotype> haplotypes) {
        final String fileName = sampleName + "-" + HAPLOTYPES_FILE_NAME + ".bed";
        final String filePath = Paths.get(outputDirectoryPath, fileName).toAbsolutePath().toString();
        try (final OutputStream outputStream = BucketUtils.createFile(filePath)) {
            outputStream.write(("#" + CalledSVGraphGenotype.bedHeader() + "\n").getBytes());
            for (final CalledSVGraphGenotype haplotype : haplotypes) {
                final List<String> lines = haplotype.bedStrings();
                for (final String line : lines) {
                    outputStream.write((line + "\n").getBytes());
                }
            }
        } catch (final IOException e) {
            throw new GATKException("Error writing haplotypes BED file", e);
        }
    }

    /**
     * Writes graph edges to tsv file
     */
    private void writeGraph(final SVGraph graph) {
        final String fileName = sampleName + "-" + GRAPH_FILE_NAME + ".tsv";
        final String filePath = Paths.get(outputDirectoryPath, fileName).toAbsolutePath().toString();
        try (final OutputStream outputStream = BucketUtils.createFile(filePath)) {
            outputStream.write(("sourceContig\ttargetContig\tsource\ttarget\tsourceStrand\ttargetStrand\tisReference\torientation\n").getBytes());
            for (final IndexedSVGraphEdge indexedEdge : graph.getEdges()) {
                final CoordinateSVGraphEdge edge = new CoordinateSVGraphEdge(indexedEdge, graph.generateNodes());
                final String orientation = (edge.isStrandA() ? "R" : "L") + (edge.isStrandB() ? "R" : "L") + (edge.isReference() ? "_r" : "");
                outputStream.write((edge.getContigA() + "\t" + edge.getContigB() + "\t" + (edge.getContigA() + ":" + edge.getNodeAPosition()) + "\t" + (edge.getContigB() + ":" + edge.getNodeBPosition()) + "\t" + edge.isStrandA() + "\t" + edge.isStrandB() + "\t" + edge.isReference() + "\t" + orientation + "\n").getBytes());
            }
        } catch (final IOException e) {
            throw new GATKException("Error writing graph to " + filePath);
        }
    }

    private CalledContigPloidyCollection readPloidyCalls() {
        return new CalledContigPloidyCollection(contigPloidyCallsFile);
    }

    private final static class EventRecordReader {

        private enum TABLE_COLUMNS {
            CHROM,
            START,
            END,
            ID,
            SAMPLES,
            TYPE
        }

        protected static Stream<EventRecord> readEventRecords(final File file) {
            try {
                final Reader reader = new FileReader(file);
                final TableReader<EventRecord> tableReader = TableUtils.reader(reader, EventRecordReader::eventRecordFactory);
                return tableReader.stream();
            } catch (final FileNotFoundException e) {
                throw new UserException(e.getMessage());
            } catch (final IOException e) {
                throw new UserException(e.getMessage());
            }
        }

        private static Function<DataLine, EventRecord> eventRecordFactory(final TableColumnCollection columns,
                                                                          final Function<String, RuntimeException> errorFunction) {
            return EventRecordReader::eventRecordParser;
        }

        private static EventRecord eventRecordParser(final DataLine data) {
            final String chrom = data.get(TABLE_COLUMNS.CHROM);
            final int start = Integer.valueOf(data.get(TABLE_COLUMNS.START));
            final int end = Integer.valueOf(data.get(TABLE_COLUMNS.END));
            final SimpleInterval interval = new SimpleInterval(chrom, start, end);
            final String[] samples = data.get(TABLE_COLUMNS.SAMPLES).split(",");
            return new EventRecord(data.get(TABLE_COLUMNS.ID), data.get(TABLE_COLUMNS.TYPE), interval, samples);
        }
    }

}
