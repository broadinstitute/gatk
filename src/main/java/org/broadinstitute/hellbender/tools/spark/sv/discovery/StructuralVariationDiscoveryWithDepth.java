package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.antlr.v4.runtime.misc.Utils;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledCopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.IntervalUtils;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SimpleSV;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.SVReadFilter;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import scala.Tuple2;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

@ExperimentalFeature
@CommandLineProgramProperties(
        oneLineSummary = "Finds large tandem duplications",
        summary = "This tool takes a breakpoint VCF, evidence target links file, and normalized copy number ratio and " +
                        "produces a VCF containing large tandem duplications.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public class StructuralVariationDiscoveryWithDepth extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private static final double MIN_TRUTH_SET_RECIPROCAL_OVERLAP = 0.5;
    private static final int TRUTH_INTERVAL_PADDING = 300;
    private static final int HIGH_DEPTH_COVERAGE_PEAK_FACTOR = 7;
    final StructuralVariationDiscoveryArgumentCollection.SimpleVariantDiscoveryArgumentCollection arguments = new StructuralVariationDiscoveryArgumentCollection.SimpleVariantDiscoveryArgumentCollection(HIGH_DEPTH_COVERAGE_PEAK_FACTOR);
    @Argument(doc = "Breakpoint VCF", fullName = "breakpoint-vcf")
    private String breakpointVCFPath;
    @Argument(doc = "Evidence target links file", fullName = "evidence-target-links-file")
    private String evidenceTargetLinksFilePath;
    @Argument(doc = "Segment copy number ratio calls", fullName = "seg-calls-file")
    private String segmentCallsFilePath;
    @Argument(doc = "Copy number ratio tsv file", fullName = "copy-ratio-file")
    private String copyRatioFilePath;
    @Argument(doc = "Reference sequence dictionary", fullName = "dict")
    private String contigDictPath;
    @Argument(doc = "Output path", fullName = "output")
    private String outputPath;
    @Argument(doc = "Assembly BAM", fullName = "assembly-bam", optional = true)
    private String assemblyBamPath;
    @Argument(doc = "Truth set VCF", fullName = "truth-set-vcf", optional = true)
    private String truthSetVCFPath;
    @Argument(doc = "True positives BED output", fullName = "true-positives-bed", optional = true)
    private String truePositivesBedPath;
    @Argument(doc = "False positives BED output", fullName = "false-positives-bed", optional = true)
    private String falsePositivesBedPath;
    @Argument(doc = "False negatives BED output", fullName = "false-negatives-bed", optional = true)
    private String falseNegativesBedPath;
    @Argument(doc = "Performance metrics file", fullName = "performance-metrics-file", optional = true)
    private String performanceMetricsFilePath;
    @Argument(doc = "Coverage output bed file path", fullName = "coverage-output-bed", optional = true)
    private String coverageOutputBedPath;
    private ReadMetadata readMetadata;
    private SAMSequenceDictionary contigDict;

    private static String getVariantContextSVType(final VariantContext vc) {
        return vc.getAttributeAsString("SVTYPE", "");
    }

    private static Stream<SimpleSV> getFalsePositives(final Collection<SimpleSV> events, final SVIntervalTree<VariantContext> truthTree, final double minFractionOverlap) {
        return events.stream().filter(event -> {
            final Stream<VariantContext> matchingTrueEvents = IntervalUtils.getTreeNodesWithReciprocalOverlap(event.getInterval(), truthTree, minFractionOverlap);
            return matchingTrueEvents.noneMatch(vc -> event.getType().toString().equals(getVariantContextSVType(vc)));
        });
    }

    public Tuple2<List<VariantContext>, List<VariantContext>> readTruthVCF(final String truthVCF) {
        List<VariantContext> insertions = new ArrayList<>();
        List<VariantContext> deletions = new ArrayList<>();
        try (final FeatureDataSource<VariantContext> dataSource =
                     new FeatureDataSource<>(truthVCF, null, 0, VariantContext.class)) {
            for (final VariantContext vc : dataSource) {
                final int contigID = readMetadata.getContigID(vc.getContig());
                if (contigID < 0) {
                    throw new UserException("VCF contig " + vc.getContig() + " does not appear in dictionary.");
                }
                final String svType = vc.getAttributeAsString("SVTYPE", "");
                if (svType.equals("insertion") || svType.equals("INS") || svType.equals("DUP") || svType.equals("CNV")) {
                    insertions.add(vc);
                } else if (svType.equals("deletion") || svType.equals("DEL")) {
                    deletions.add(vc);
                }
            }
        }
        return new Tuple2<>(insertions, deletions);
    }

    public List<VariantContext> readVCF(final String vcfPath) {
        List<VariantContext> variants = new ArrayList<>();
        try (final FeatureDataSource<VariantContext> dataSource =
                     new FeatureDataSource<>(vcfPath, null, 0, VariantContext.class)) {
            for (final VariantContext vc : dataSource) {
                final int contigID = readMetadata.getContigID(vc.getContig());
                if (contigID < 0) {
                    throw new UserException("VCF contig " + vc.getContig() + " does not appear in dictionary.");
                }
                variants.add(vc);
            }
        }
        return variants;
    }

    public void runTool(final JavaSparkContext ctx) {

        if (truthSetVCFPath != null || truePositivesBedPath != null || falsePositivesBedPath != null || falseNegativesBedPath != null || performanceMetricsFilePath != null) {
            if (truthSetVCFPath == null || truePositivesBedPath == null || falsePositivesBedPath == null || falseNegativesBedPath == null || performanceMetricsFilePath == null) {
                throw new IllegalArgumentException("Either all or none of the truth set arguments must be specified");
            }
        }

        contigDict = ReferenceUtils.loadFastaDictionary(new File(contigDictPath));
        final SAMFileHeader header = new SAMFileHeader(contigDict);
        final SVReadFilter svReadFilter = new SVReadFilter(new StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection());
        final JavaRDD<GATKRead> singletonReadRdd = ctx.parallelize(Collections.singletonList(ArtificialReadUtils.createRandomRead(151))); //To avoid divide by zero error in ReadMetadata
        readMetadata = new ReadMetadata(Collections.emptySet(), header, 2000, singletonReadRdd, svReadFilter, logger);

        final List<GATKRead> assembly;
        if (assemblyBamPath != null) {
            final ReadsSparkSource readsSource = new ReadsSparkSource(ctx);
            assembly = readsSource.getParallelReads(assemblyBamPath, null).collect();
        } else {
            assembly = Collections.emptyList();
        }

        final List<VariantContext> breakpoints = readVCF(breakpointVCFPath);
        final File readDepthFile = new File(copyRatioFilePath);
        final File gCNVCallsFile = new File(segmentCallsFilePath);

        final CalledCopyRatioSegmentCollection copyRatioSegments = new CalledCopyRatioSegmentCollection(gCNVCallsFile);

        final List<EvidenceTargetLink> evidenceTargetLinks = new ArrayList<>();
        try {
            final String[] evidenceTargetFileLines = String.valueOf(Utils.readFile(evidenceTargetLinksFilePath)).split("\n");
            for (final String line : evidenceTargetFileLines) {
                evidenceTargetLinks.add(EvidenceTargetLink.fromBedpeString(line.trim(), readMetadata));
            }
        } catch (final IOException e) {
            throw new GATKException("Could not read evidence links file", e);
        }

        final CopyRatioCollection readDepthData = new CopyRatioCollection(readDepthFile);
        final LargeSimpleEventCaller largeSimpleEventCaller = new LargeSimpleEventCaller(breakpoints, assembly, evidenceTargetLinks, readDepthData, copyRatioSegments, readMetadata, contigDict, arguments);
        final List<SimpleSV> events = largeSimpleEventCaller.getEvents();
        LargeSimpleEventCaller.writeTandemDuplicationEventsAsBedFile(outputPath, events, readMetadata);
        compareCallsToTruthSet(events, truthSetVCFPath, truePositivesBedPath, falsePositivesBedPath, falseNegativesBedPath, performanceMetricsFilePath, largeSimpleEventCaller);

    }

    public void compareCallsToTruthSet(final Collection<SimpleSV> events,
                                       final String truthSetVCFPath,
                                       final String truePositivesPath,
                                       final String falsePositivesPath,
                                       final String falseNegativesPath,
                                       final String performanceMetricsPath,
                                       final LargeSimpleEventCaller eventCaller) {
        final List<VariantContext> truthVariants = readVCF(truthSetVCFPath);
        final SVIntervalTree<VariantContext> truthInsertionIntervalTree = createIntervalTreeFromVariantList(truthVariants);

        final SVIntervalTree<SimpleSV> calledEventTree = new SVIntervalTree<>();
        for (final SimpleSV event : events) {
            final SVInterval eventInterval = new SVInterval(event.getContigId(), event.getStart(), event.getEnd());
            calledEventTree.put(eventInterval, event);
        }

        final Map<String,Integer> numFalsePositives;
        Map<String,Integer> numTruePositives = new HashMap<>();
        Map<String,Integer> numFalseNegatives = new HashMap<>();
        try (final OutputStream truePositivesStream = BucketUtils.createFile(truePositivesPath);
             final OutputStream falsePositivesStream = BucketUtils.createFile(falsePositivesPath);
             final OutputStream falseNegativesStream = BucketUtils.createFile(falseNegativesPath)) {
            final List<SimpleSV> falsePositives = getFalsePositives(events, truthInsertionIntervalTree, MIN_TRUTH_SET_RECIPROCAL_OVERLAP).collect(Collectors.toList());
            numFalsePositives = countTypes(falsePositives);
            for (final SimpleSV event : falsePositives) {
                falsePositivesStream.write(getCallString(event).getBytes());
            }
            for (final VariantContext vc : truthVariants) {
                final int contig = readMetadata.getContigID(vc.getContig());
                final SVInterval truthInterval = new SVInterval(contig, vc.getStart(), vc.getEnd());
                final String vcType = getVariantContextSVType(vc);
                final List<SimpleSV> overlappingEvents = IntervalUtils.getTreeNodesWithReciprocalOverlap(truthInterval, calledEventTree, MIN_TRUTH_SET_RECIPROCAL_OVERLAP).collect(Collectors.toList());
                final List<SimpleSV> typedOverlappingEvents = overlappingEvents.stream().filter(event -> event.getType().toString().equals(getVariantContextSVType(vc))).collect(Collectors.toList());
                final Optional<SimpleSV> bestEvent = LargeSimpleEventCaller.getHighestScoringEventFromStream(typedOverlappingEvents.stream(), arguments.COUNTEREVIDENCE_PSEUDOCOUNT);
                if (!bestEvent.isPresent()) {
                    if (truthInterval.getLength() >= arguments.MIN_SV_SIZE + 100 && truthInterval.getLength() <= arguments.MAX_SV_SIZE) {
                        numFalseNegatives.putIfAbsent(vcType, 0);
                        numFalseNegatives.put(vcType, numFalseNegatives.get(vcType) + 1);
                        final SVInterval leftBreakpointInterval = new SVInterval(truthInterval.getContig(), truthInterval.getStart(), truthInterval.getStart());
                        final SVInterval rightBreakpointInterval = new SVInterval(truthInterval.getContig(), truthInterval.getEnd(), truthInterval.getEnd());
                        final SVInterval leftInterval = IntervalUtils.getPaddedInterval(leftBreakpointInterval, TRUTH_INTERVAL_PADDING, contigDict);
                        final SVInterval rightInterval = IntervalUtils.getPaddedInterval(rightBreakpointInterval, TRUTH_INTERVAL_PADDING, contigDict);
                        final Optional<SimpleSV> matchingEvent = eventCaller.getHighestScoringEventOnInterval(leftInterval, rightInterval, truthInterval, new SVIntervalTree<>(), null, false);
                        if (matchingEvent.isPresent()) {
                            falseNegativesStream.write(getCallString(matchingEvent.get()).getBytes());
                        } else {
                            final SimpleSVType.TYPES type = SimpleSVType.TYPES.valueOf(vcType);
                            falseNegativesStream.write(getCallString(new SimpleSV(type, vc.getStart(), vc.getEnd(), contig, vc.getContig(), 0, 0, 0, 0, Collections.emptyList(), Collections.emptyList(), null)).getBytes());
                        }
                    }
                } else {
                    final byte[] bedString = getCallString(bestEvent.get()).getBytes();
                    truePositivesStream.write(bedString);
                    numTruePositives.putIfAbsent(vcType, 0);
                    numTruePositives.put(vcType, numTruePositives.get(vcType) + 1);
                }
            }
        } catch (final IOException e) {
            throw new GATKException("Error writing truth set BED files", e);
        }

        try (final OutputStream metricsStream = BucketUtils.createFile(performanceMetricsPath)) {
            List<String> types = Arrays.asList(SimpleSVType.TYPES.values()).stream().map(type -> type.toString()).collect(Collectors.toList());
            metricsStream.write(("TYPE\tTP\tFP\tFN\tRECALL\tFDR\n").getBytes());
            int totalTP = 0;
            int totalFP = 0;
            int totalFN = 0;
            for (final String type : types) {
                final int tp = numTruePositives.getOrDefault(type, 0);
                final int fp = numFalsePositives.getOrDefault(type, 0);
                final int fn = numFalseNegatives.getOrDefault(type, 0);
                totalTP += tp;
                totalFP += fp;
                totalFN += fn;
                final double recall = tp / (double) (tp + fn);
                final double fdr = fp / (double) (fp + tp);
                metricsStream.write((type + "\t" + tp + "\t" + fp + "\t" + fn + "\t" + recall + "\t" + fdr +"\n").getBytes());
            }
            final double totalRecall = totalTP / (double) (totalTP + totalFN);
            final double totalFDR = totalFP / (double) (totalFP + totalTP);
            metricsStream.write(("TOTAL\t" + totalTP + "\t" + totalFP + "\t" + totalFN + "\t" + totalRecall + "\t" + totalFDR + "\n").getBytes());
        } catch (final IOException e) {
            throw new GATKException("Error writing performance metrics file", e);
        }
    }

    private Map<String,Integer> countTypes(final List<SimpleSV> events) {
        return events.stream().map(event -> new HashMap.SimpleImmutableEntry<>(event.getType().toString(),1)).collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue, (a, b) -> a + b));
    }

    private String getCallString(final SimpleSV event) {
        return event.getString(readMetadata) + "\n";
    }

    private SVIntervalTree<VariantContext> createIntervalTreeFromVariantList(final List<VariantContext> list) {
        final SVIntervalTree<VariantContext> tree = new SVIntervalTree<>();
        for (final VariantContext vc : list) {
            final SVInterval interval = new SVInterval(readMetadata.getContigID(vc.getContig()), vc.getStart(), vc.getEnd());
            tree.put(interval, vc);
        }
        return tree;
    }
}
