package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.ProgressMeter;
import org.broadinstitute.hellbender.engine.TraversalParameters;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledCopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatioSegment;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Searches for large structural variants using assembled breakpoint pairs, clustered read pair evidence, binned copy
 * ratios, and copy ratio segments. The search is performed by iterating through all the breakpoint pairs and clustered
 * read pairs and searching for an event in their vicinity.
 */
public class LargeSimpleSVCaller {

    private static final Logger logger = LogManager.getLogger(LargeSimpleSVCaller.class);

    private final JavaRDD<GATKRead> reads;
    private final SVIntervalTree<LargeSimpleSV> largeSimpleSVSVIntervalTree;
    private final SVIntervalTree<VariantContext> truthSetTree;
    private final SVIntervalTree<GATKRead> contigTree;
    private final List<Collection<SVCopyRatio>> copyRatios;
    private final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector;
    private final SAMSequenceDictionary dictionary;
    private final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromReadDepthArgumentCollection arguments;

    public LargeSimpleSVCaller(final JavaRDD<GATKRead> reads,
                               final Collection<LargeSimpleSV> largeSimpleSVCollection,
                               final Collection<VariantContext> truthSet,
                               final Collection<GATKRead> assembledContigs,
                               final CopyRatioCollection copyRatios,
                               final CalledCopyRatioSegmentCollection copyRatioSegments,
                               final SAMSequenceDictionary dictionary,
                               final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromReadDepthArgumentCollection arguments) {
        Utils.nonNull(reads, "Reads RDD cannot be null");
        Utils.nonNull(largeSimpleSVCollection, "SV collection cannot be null");
        Utils.nonNull(assembledContigs, "Contig collection cannot be null");
        Utils.nonNull(copyRatios, "Copy ratio collection cannot be null");
        Utils.nonNull(copyRatioSegments, "Copy ratio segments collection cannot be null");
        Utils.nonNull(dictionary, "Dictionary cannot be null");
        Utils.nonNull(arguments, "Parameter arguments collection cannot be null");

        this.reads = reads;
        this.dictionary = dictionary;
        this.arguments = arguments;

        logger.info("Building interval trees...");

        if (truthSet != null) {
            truthSetTree = buildVariantIntervalTree(truthSet);
        } else {
            truthSetTree = null;
        }
        contigTree = buildReadIntervalTree(assembledContigs);
        /* copyRatioOverlapDetector =  getMinimalCopyRatioCollection(copyRatios, copyRatios.getMetadata(),
                filteredEvidenceTargetLinks, pairedBreakpoints, dictionary, arguments.smallEventSize, MAX_COPY_RATIO_EVENT_SIZE,
                arguments.breakpointPadding + arguments.hmmPadding,
                arguments.evidenceTargetLinkPadding + arguments.hmmPadding).getOverlapDetector();*/

        //Fill empty intervals with neutral calls
        final List<CalledCopyRatioSegment> segments = new ArrayList<>(copyRatioSegments.getRecords());
        Collections.sort(segments, IntervalUtils.getDictionaryOrderComparator(dictionary));
        final List<CalledCopyRatioSegment> emptySegments = new ArrayList<>(segments.size());
        for (int i = 0; i < segments.size() - 1; i++) {
            final CalledCopyRatioSegment currentSegment = segments.get(i);
            final CalledCopyRatioSegment nextSegment = segments.get(i + 1);
            if (currentSegment.getContig().equals(nextSegment.getContig())
                    && currentSegment.getEnd() + 1 < nextSegment.getStart()) {
                final SimpleInterval newInterval = new SimpleInterval(currentSegment.getContig(), currentSegment.getEnd() + 1, nextSegment.getStart() - 1);
                emptySegments.add(new CalledCopyRatioSegment(new CopyRatioSegment(newInterval, 0, 0), CalledCopyRatioSegment.Call.NEUTRAL));
            }
        }
        final List<CalledCopyRatioSegment> filledSegments = new ArrayList<>(segments.size() + emptySegments.size());
        filledSegments.addAll(segments);
        filledSegments.addAll(emptySegments);
        final CalledCopyRatioSegmentCollection filledSegmentCollection = new CalledCopyRatioSegmentCollection(copyRatioSegments.getMetadata(), filledSegments);
        copyRatioSegmentOverlapDetector = filledSegmentCollection.getOverlapDetector();

        //Filter calls without overlapping segment support
        final Collection<LargeSimpleSV> filteredCalls = largeSimpleSVCollection.stream().filter(event -> {
            final Set<CalledCopyRatioSegment> overlappers = copyRatioSegmentOverlapDetector.getOverlaps(SVIntervalUtils.convertToSimpleInterval(event.getInterval(), dictionary));
            return overlappers.stream().anyMatch(segment -> segment.getCall() == (event.getEventType() == SimpleSVType.TYPES.DEL ? CalledCopyRatioSegment.Call.DELETION : CalledCopyRatioSegment.Call.AMPLIFICATION));
        }).collect(Collectors.toList());
        largeSimpleSVSVIntervalTree = buildLargeSimpleSVTree(filteredCalls);

        logger.info("Partitioning copy ratios by contig...");
        this.copyRatios = new ArrayList<>(dictionary.size());
        for (int i = 0; i < dictionary.size(); i++) {
            final double estimatedBinFraction = dictionary.getSequence(i).getSequenceLength() / (double) dictionary.getReferenceLength();
            final int estimatedNumBins = (int) (estimatedBinFraction * copyRatios.size());
            this.copyRatios.add(new ArrayList<>(estimatedNumBins));
        }
        if (!copyRatios.getMetadata().getSequenceDictionary().isSameDictionary(dictionary)) {
            throw new UserException.IncompatibleSequenceDictionaries("Copy ratio dictionary does not match sequence dictionary", "copy ratio", copyRatios.getMetadata().getSequenceDictionary(), "master", dictionary);
        }
        for (final CopyRatio copyRatio : copyRatios.getRecords()) {
            final int copyRatioContigIndex = dictionary.getSequenceIndex(copyRatio.getContig());
            if (copyRatioContigIndex == -1) {
                throw new UserException.BadInput("Copy ratio and master sequence dictionaries matched, but encountered a copy ratio with contig " + copyRatio.getContig() + " with no record");
            }
            this.copyRatios.get(copyRatioContigIndex).add(new SVCopyRatio(copyRatioContigIndex, copyRatio.getStart(), copyRatio.getEnd(), (float) copyRatio.getLog2CopyRatioValue()));
        }
    }

    /**
     * Gets minimal list of copy ratio bins that overlap relevant evidence.
     */
    public static CopyRatioCollection getMinimalCopyRatioCollection(final CopyRatioCollection copyRatioCollection,
                                                                    final SampleLocatableMetadata copyRatioMetadata,
                                                                    final Collection<EvidenceTargetLink> evidenceTargetLinks,
                                                                    final Collection<IntrachromosomalBreakpointPair> pairedBreakpoints,
                                                                    final SAMSequenceDictionary dictionary,
                                                                    final int minIntervalSize,
                                                                    final int maxIntervalSize,
                                                                    final int breakpointIntervalPadding,
                                                                    final int linkIntervalPadding) {

        //Map evidence to padded intervals
        final Stream<SVInterval> breakpointIntervalStream = pairedBreakpoints.stream().map(IntrachromosomalBreakpointPair::getInterval)
                .filter(interval -> interval.getLength() <= maxIntervalSize && interval.getLength() >= minIntervalSize)
                .map(interval -> SVIntervalUtils.getPaddedInterval(interval, breakpointIntervalPadding, dictionary));
        final Stream<SVInterval> linkIntervalStream = evidenceTargetLinks.stream()
                .filter(link -> link.getPairedStrandedIntervals().getLeft().getInterval().getContig() == link.getPairedStrandedIntervals().getRight().getInterval().getContig())
                .map(SVIntervalUtils::getOuterIntrachromosomalLinkInterval)
                .filter(interval -> interval.getLength() <= maxIntervalSize && interval.getLength() >= minIntervalSize)
                .map(interval -> SVIntervalUtils.getPaddedInterval(interval, linkIntervalPadding, dictionary));

        //Merge streams and convert intervals to GenomeLoc
        final List<GenomeLoc> intervalList = Stream.concat(breakpointIntervalStream, linkIntervalStream)
                .map(interval -> SVIntervalUtils.convertToGenomeLoc(interval, dictionary))
                .collect(Collectors.toList());

        //Merge intervals and build into a tree
        Collections.sort(intervalList, IntervalUtils.getDictionaryOrderComparator(dictionary));
        final List<GenomeLoc> mergedIntervals = IntervalUtils.mergeIntervalLocations(intervalList, IntervalMergingRule.ALL);
        final SVIntervalTree<Object> mergedIntervalTree = new SVIntervalTree<>();
        for (final GenomeLoc loc : mergedIntervals) {
            mergedIntervalTree.put(new SVInterval(loc.getContigIndex(), loc.getStart(), loc.getEnd()), null);
        }

        //Return copy ratios overlapping evidence intervals
        final List<CopyRatio> countsList = copyRatioCollection.getRecords().stream()
                .filter(copyRatio -> mergedIntervalTree.hasOverlapper(SVIntervalUtils.convertToSVInterval(copyRatio.getInterval(), dictionary)))
                .collect(Collectors.toList());
        Collections.sort(countsList, IntervalUtils.getDictionaryOrderComparator(dictionary));
        return new CopyRatioCollection(copyRatioMetadata, countsList);
    }

    /**
     * Returns true if the two BNDs in vc1 and vc2 are a valid breakpoint pair, as indicated by their MATEID attributes
     */
    private static boolean isBreakpointPair(final VariantContext vc1, final VariantContext vc2) {
        return vc1.hasAttribute(GATKSVVCFConstants.BND_MATEID_STR) &&
                vc1.getAttributeAsString(GATKSVVCFConstants.BND_MATEID_STR, "").equals(vc2.getID()) &&
                vc2.getAttributeAsString(GATKSVVCFConstants.BND_MATEID_STR, "").equals(vc1.getID());
    }

    /**
     * Returns list of ordered CopyRatio objects on the given interval
     *
     * @param interval      Interval over which to retrieve bins
     * @param copyRatioTree Tree containing all overlapping copy ratios
     * @param binsToTrim    Number of bins to trim from either side
     * @param dictionary    Sequence dictionary
     * @return List of copy ratios
     */
    @VisibleForTesting
    static List<SVCopyRatio> getCopyRatiosOnInterval(final SVInterval interval, final SVIntervalTree<Float> copyRatioTree,
                                                     final int binsToTrim, final SAMSequenceDictionary dictionary) {
        final SimpleInterval simpleInterval = SVIntervalUtils.convertToSimpleInterval(interval, dictionary);
        if (simpleInterval.size() == 0) {
            return Collections.emptyList();
        }
        final List<SVCopyRatio> copyRatios = Utils.stream(copyRatioTree.overlappers(interval))
                .map(entry -> new SVCopyRatio(entry.getInterval(), entry.getValue())).collect(Collectors.toList());

        if (copyRatios.size() <= 2 * binsToTrim) return Collections.emptyList();
        Collections.sort(copyRatios, Comparator.comparing(SVCopyRatio::getStart));
        return copyRatios.subList(binsToTrim, copyRatios.size() - binsToTrim);
    }

    private static double fractionEmpty(final SVInterval interval, final List<SVCopyRatio> copyRatios) {
        return 1.0 - (copyRatios.stream().mapToInt(copyRatio -> copyRatio.getInterval().getLength()).sum() / (double) interval.getLength());
    }

    /**
     * Builds an interval tree containing the only the aligned reads
     */
    private SVIntervalTree<GATKRead> buildReadIntervalTree(final Collection<GATKRead> reads) {
        final SVIntervalTree<GATKRead> tree = new SVIntervalTree<>();
        for (final GATKRead read : reads) {
            if (read.isUnmapped()) continue;
            final int start = read.getStart();
            final int end = read.getEnd();
            final int contig = dictionary.getSequenceIndex(read.getContig());
            tree.put(new SVInterval(contig, start, end), read);
        }
        return tree;
    }

    /**
     * Builds an interval tree of SVIntervals with null entry values
     */
    private SVIntervalTree<Object> buildSVIntervalTree(final Collection<SVInterval> intervalCollection) {
        final SVIntervalTree<Object> tree = new SVIntervalTree<>();
        for (final SVInterval interval : intervalCollection) {
            tree.put(new SVInterval(interval.getContig(), interval.getStart(), interval.getEnd()), null);
        }
        return tree;
    }

    /**
     * Builds an interval tree containing the variant calls
     */
    private SVIntervalTree<LargeSimpleSV> buildLargeSimpleSVTree(final Collection<LargeSimpleSV> largeSimpleSVCollection) {
        final SVIntervalTree<LargeSimpleSV> tree = new SVIntervalTree<>();
        for (final LargeSimpleSV sv : largeSimpleSVCollection) {
            tree.put(sv.getInterval(), sv);
        }
        return tree;
    }

    /**
     * Builds an interval tree containing the variant calls
     */
    private SVIntervalTree<VariantContext> buildVariantIntervalTree(final Collection<VariantContext> variantContexts) {
        final SVIntervalTree<VariantContext> tree = new SVIntervalTree<>();
        for (final VariantContext vc : variantContexts) {
            final int start = vc.getStart();
            final int end = vc.getEnd();
            final int contig = dictionary.getSequenceIndex(vc.getContig());
            tree.put(new SVInterval(contig, start, end), vc);
        }
        return tree;
    }

    /**
     * Builds an interval tree of target evidence links, whose intervals may be padded. Tree may be built using
     * single intervals for each link (start of left to end of right) or two intervals (one for left and one for right)
     *
     * @param links                      Target evidence links used to build the tree
     * @param padding                    Padding applied to intervals
     * @param separateLeftRightIntervals If true, inserts intervals for left and right intervals
     * @return The interval tree
     */
    private SVIntervalTree<EvidenceTargetLink> buildEvidenceIntervalTree(final Collection<EvidenceTargetLink> links, final int padding, final boolean separateLeftRightIntervals) {
        final SVIntervalTree<EvidenceTargetLink> linkTree = new SVIntervalTree<>();
        for (final EvidenceTargetLink link : links) {
            final SVInterval linkLeftInterval = link.getPairedStrandedIntervals().getLeft().getInterval();
            final SVInterval linkRightInterval = link.getPairedStrandedIntervals().getRight().getInterval();
            if (separateLeftRightIntervals) {
                final SVInterval eventIntervalLeft = new SVInterval(linkLeftInterval.getContig(), linkLeftInterval.getStart(), linkLeftInterval.getEnd());
                final SVInterval paddedEventIntervalLeft = SVIntervalUtils.getPaddedInterval(eventIntervalLeft, padding, dictionary);
                linkTree.put(paddedEventIntervalLeft, link);
                final SVInterval eventIntervalRight = new SVInterval(linkRightInterval.getContig(), linkRightInterval.getStart(), linkRightInterval.getEnd());
                final SVInterval paddedEventIntervalRight = SVIntervalUtils.getPaddedInterval(eventIntervalRight, padding, dictionary);
                linkTree.put(paddedEventIntervalRight, link);
            } else if (linkLeftInterval.getContig() == linkRightInterval.getContig()) {
                final SVInterval eventInterval = new SVInterval(linkLeftInterval.getContig(), linkLeftInterval.getStart(), linkRightInterval.getEnd());
                final SVInterval paddedEventInterval = SVIntervalUtils.getPaddedInterval(eventInterval, padding, dictionary);
                linkTree.put(paddedEventInterval, link);
            }
        }
        return linkTree;
    }

    /**
     * Returns all events. Searches by iterating over the breakpoint pairs and then the evidence target links.
     */
    public Tuple2<List<ReadDepthEvent>, ReadDepthModel.ReadDepthModelParameters> callEvents(final String inputPath, final JavaSparkContext ctx, final ProgressMeter progressMeter) {


        final SVIntervalTree<GATKRead> readsTree = new SVIntervalTree<>();
        /*final List<SimpleInterval> eventIntervals = Utils.stream(largeSimpleSVSVIntervalTree.iterator())
                .map(SVIntervalTree.Entry::getInterval)
                .map(interval -> SVIntervalUtils.convertToSimpleInterval(interval, dictionary))
                .collect(Collectors.toList());
        final List<Tuple2<SVInterval,GATKRead>> intervalsAndReads = getReads(inputPath, ctx, eventIntervals, dictionary);
        logger.info("Collecting reads...");
        for (final Tuple2<SVInterval,GATKRead> pair : intervalsAndReads) {
            readsTree.put(pair._1, pair._2);
        }
        logger.info("Retrieved " + intervalsAndReads.stream().mapToInt(pair -> pair._2.getLength()).sum() + " reads");
        */

        logger.info("Running read depth model on " + largeSimpleSVSVIntervalTree.size() + " events");
        final ReadDepthModel readDepthModel = new ReadDepthModel(readsTree, largeSimpleSVSVIntervalTree, copyRatioSegmentOverlapDetector, dictionary);
        if (truthSetTree == null) {
            final List<ReadDepthEvent> result = readDepthModel.solve(ctx);
            return new Tuple2<>(result, readDepthModel.getParameters());
        }
        return readDepthModel.train(ctx, truthSetTree);
    }

    public static List<Tuple2<SVInterval,GATKRead>> getReads(final String inputPath, final JavaSparkContext ctx, final List<SimpleInterval> intervals, final SAMSequenceDictionary dictionary) {
        final TraversalParameters traversalParameters = new TraversalParameters(intervals, false);
        final ReadsSparkSource readsSource = new ReadsSparkSource(ctx, ValidationStringency.DEFAULT_STRINGENCY);
        final JavaRDD<GATKRead> reads = readsSource.getParallelReads(inputPath, null, traversalParameters, 8000000);
        return reads.filter(read -> !read.isUnmapped()).map(read -> new Tuple2<>(SVIntervalUtils.convertToSVInterval(new SimpleInterval(read.getContig(), read.getUnclippedStart(), read.getUnclippedEnd()), dictionary), read)).collect();
    }

}
