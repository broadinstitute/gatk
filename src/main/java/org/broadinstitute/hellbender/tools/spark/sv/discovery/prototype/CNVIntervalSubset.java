package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.VariantContext;
import org.antlr.v4.runtime.misc.Utils;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.SVReadFilter;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import scala.Tuple2;

import java.io.File;
import java.io.IOException;
import java.util.*;

@CommandLineProgramProperties(
        oneLineSummary = "Takes subset of CNV read count intervals",
        summary = "",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public class CNVIntervalSubset extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private static final int INTERVAL_PADDING = 10000;
    private static final int MAX_INTERVAL_SIZE = 100000;
    private static final int MIN_INTERVAL_SIZE = 500;
    @Argument(doc = "breakpoint VCF", fullName = "breakpoint-vcf")
    private String breakpointVCFPath;
    @Argument(doc = "evidence target links file", fullName = "evidence-target-links-file")
    private String evidenceTargetLinksFilePath;
    @Argument(doc = "Read depth file (hdf5 or tsv)", fullName = "read-depth-file")
    private String readDepthFilePath;
    @Argument(doc = "Output file path", fullName = "output")
    private String outputPath;
    @Argument(doc = "Contig dictionary", fullName = "contig-dict")
    private String contigDictPath;

    public void runTool(final JavaSparkContext ctx) {

        final SAMSequenceDictionary contigDict = ReferenceUtils.loadFastaDictionary(new File(contigDictPath));
        final List<VariantContext> breakpoints = readVCF(breakpointVCFPath, contigDict);

        final SAMFileHeader header = new SAMFileHeader(contigDict);
        final List<EvidenceTargetLink> evidenceTargetLinks = new ArrayList<>();
        final SVReadFilter svReadFilter = new SVReadFilter(new StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection());
        final JavaRDD<GATKRead> singletonReadRdd = ctx.parallelize(Collections.singletonList(ArtificialReadUtils.createRandomRead(151))); //To avoid divide by zero error in ReadMetadata
        final ReadMetadata readMetadata = new ReadMetadata(Collections.emptySet(), header, 2000, singletonReadRdd, svReadFilter, logger);
        try {
            final String[] evidenceTargetFileLines = String.valueOf(Utils.readFile(evidenceTargetLinksFilePath)).split("\n");
            for (final String line : evidenceTargetFileLines) {
                evidenceTargetLinks.add(EvidenceTargetLink.fromBedpeString(line.trim(), readMetadata));
            }
        } catch (final IOException e) {
            throw new GATKException("Could not read evidence links file", e);
        }

        final Map<String, VariantContext> unpairedVariants = new HashMap<>();
        final List<Tuple2<VariantContext, Integer>> pairedBreakpoints = new ArrayList<>();
        final Iterator<VariantContext> breakpointIter = breakpoints.iterator();
        while (breakpointIter.hasNext()) {
            final VariantContext vc1 = breakpointIter.next();
            final String mate = vc1.getAttributeAsString("MATEID", "");
            if (unpairedVariants.containsKey(mate)) {
                final VariantContext vc2 = unpairedVariants.remove(mate);
                if (isBreakpointPair(vc1, vc2)) {
                    if (!vc1.getContig().equals(vc2.getContig())) {
                        continue;
                    }
                    if (vc1.getStart() < vc2.getStart()) {
                        pairedBreakpoints.add(new Tuple2<>(vc1, vc2.getStart()));
                    } else {
                        pairedBreakpoints.add(new Tuple2<>(vc2, vc1.getStart()));
                    }
                } else {
                    throw new IllegalStateException("Variant mate attributes did not match: " + vc1 + "\t" + vc2);
                }
            } else {
                unpairedVariants.put(vc1.getID(), vc1);
            }
        }
        if (!unpairedVariants.isEmpty()) {
            logger.warn("There were " + unpairedVariants.size() + " unpaired breakpoint variants");
        }

        final List<GenomeLoc> intervals = new ArrayList<>(evidenceTargetLinks.size() + breakpoints.size());
        for (final Tuple2<VariantContext,Integer> breakpoint : pairedBreakpoints) {
            final VariantContext vc = breakpoint._1;
            final int contig = contigDict.getSequence(vc.getContig()).getSequenceIndex();
            final int start = vc.getStart();
            final int end = breakpoint._2;
            final int length = end - start;
            if (length <= MAX_INTERVAL_SIZE && length >= MIN_INTERVAL_SIZE) {
                final int paddedStart = Math.max(0, start - INTERVAL_PADDING);
                final int paddedEnd = Math.min(contigDict.getSequence(contig).getSequenceLength(), end + INTERVAL_PADDING);
                intervals.add(new GenomeLoc(vc.getContig(), contig, paddedStart, paddedEnd));
            }
        }
        for (final EvidenceTargetLink link : evidenceTargetLinks) {
            final int leftContig = link.getPairedStrandedIntervals().getLeft().getInterval().getContig();
            final int rightContig = link.getPairedStrandedIntervals().getRight().getInterval().getContig();
            final boolean leftStrand = link.getPairedStrandedIntervals().getLeft().getStrand();
            final boolean rightStrand = link.getPairedStrandedIntervals().getRight().getStrand();
            if (leftContig == rightContig && !leftStrand && rightStrand) {
                final int start = link.getPairedStrandedIntervals().getLeft().getInterval().getStart();
                final int end = link.getPairedStrandedIntervals().getRight().getInterval().getEnd();
                final int length = end - start;
                if (length <= MAX_INTERVAL_SIZE && length >= MIN_INTERVAL_SIZE) {
                    final int paddedStart = Math.max(0, start - INTERVAL_PADDING);
                    final int paddedEnd = Math.min(contigDict.getSequence(leftContig).getSequenceLength(), end + INTERVAL_PADDING);
                    intervals.add(new GenomeLoc(contigDict.getSequence(leftContig).getSequenceName(), leftContig, paddedStart, paddedEnd));
                }
            }
        }
        Collections.sort(intervals, IntervalUtils.getDictionaryOrderComparator(contigDict));
        final List<GenomeLoc> mergedIntervals = IntervalUtils.mergeIntervalLocations(intervals, IntervalMergingRule.ALL);

        final boolean bSimpleCount = readDepthFilePath.endsWith(".hdf5");
        final File readDepthFile = new File(readDepthFilePath);
        if (bSimpleCount) {
            final SimpleCountCollection readDepthData = SimpleCountCollection.read(readDepthFile);
            final OverlapDetector<SimpleCount> overlapDetector = readDepthData.getOverlapDetector();
            final List<SimpleCount> countsList = new ArrayList<>(mergedIntervals.size());
            for (final GenomeLoc loc : mergedIntervals) {
                final Set<SimpleCount> countSet = overlapDetector.getOverlaps(loc);
                countsList.addAll(countSet);
            }
            final SampleLocatableMetadata sampleLocatableMetadata = new SimpleSampleLocatableMetadata(readDepthFile.getName(), contigDict);
            Collections.sort(countsList, IntervalUtils.getDictionaryOrderComparator(contigDict));

            final SimpleCountCollection resultCollection = new SimpleCountCollection(sampleLocatableMetadata, countsList);
            resultCollection.writeHDF5(new File(outputPath));
        } else {
            final CopyRatioCollection readDepthData = new CopyRatioCollection(readDepthFile);
            final OverlapDetector<CopyRatio> overlapDetector = readDepthData.getOverlapDetector();
            final List<CopyRatio> countsList = new ArrayList<>(mergedIntervals.size());
            for (final GenomeLoc loc : mergedIntervals) {
                final Set<CopyRatio> countSet = overlapDetector.getOverlaps(loc);
                countsList.addAll(countSet);
            }
            final SampleLocatableMetadata sampleLocatableMetadata = new SimpleSampleLocatableMetadata(readDepthFile.getName(), contigDict);
            Collections.sort(countsList, IntervalUtils.getDictionaryOrderComparator(contigDict));

            final CopyRatioCollection resultCollection = new CopyRatioCollection(sampleLocatableMetadata, countsList);
            resultCollection.write(new File(outputPath));
        }
    }

    public static List<VariantContext> readVCF(final String vcfPath,
                                               final SAMSequenceDictionary dictionary) {
        List<VariantContext> variants = new ArrayList<>();
        try (final FeatureDataSource<VariantContext> dataSource =
                     new FeatureDataSource<>(vcfPath, null, 0, VariantContext.class)) {
            for (final VariantContext vc : dataSource) {
                final int contigID = dictionary.getSequenceIndex(vc.getContig());
                if (contigID < 0) {
                    throw new UserException("VCF contig " + vc.getContig() + " does not appear in dictionary.");
                }
                variants.add(vc);
            }
        }
        return variants;
    }

    private static boolean isBreakpointPair(final VariantContext vc1, final VariantContext vc2) {
        return vc1.getAttributeAsString("MATEID", "").equals(vc2.getID()) &&
                vc2.getAttributeAsString("MATEID", "").equals(vc1.getID());
    }
}