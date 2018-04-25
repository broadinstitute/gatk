package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.math3.linear.RealVector;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.ProgressMeter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledCopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatioSegment;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.LargeSimpleSV;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.hmm.ViterbiAlgorithm;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Searches for large structural variants using assembled breakpoint pairs, clustered read pair evidence, binned copy
 * ratios, and copy ratio segments. The search is performed by iterating through all the breakpoint pairs and clustered
 * read pairs and searching for an event in their vicinity.
 */
public class LargeSimpleSVCaller {

    private static final Logger logger = LogManager.getLogger(LargeSimpleSVCaller.class);

    private final SVIntervalTree<Object> highCoverageIntervalTree;
    private final SVIntervalTree<VariantContext> structuralVariantCallTree;
    private final SVIntervalTree<EvidenceTargetLink> intrachromosomalLinkTree;
    private final SVIntervalTree<EvidenceTargetLink> interchromosomalLinkTree;
    private final SVIntervalTree<GATKRead> contigTree;
    private final List<Collection<SVCopyRatio>> copyRatios;
    private final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector;
    private final LargeSimpleSVFactory tandemDuplicationFactory;
    private final LargeDeletionFactory largeDeletionFactory;
    private final DispersedDuplicationFactory dispersedDuplicationFactory;
    private final SAMSequenceDictionary dictionary;
    private final Collection<IntrachromosomalBreakpointPair> pairedBreakpoints;
    private final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromReadDepthArgumentCollection arguments;

    public LargeSimpleSVCaller(final Collection<VariantContext> breakpoints,
                               final Collection<VariantContext> structuralVariantCalls,
                               final Collection<GATKRead> assembledContigs,
                               final Collection<EvidenceTargetLink> evidenceTargetLinks,
                               final CopyRatioCollection copyRatios,
                               final CalledCopyRatioSegmentCollection copyRatioSegments,
                               final List<SVInterval> highCoverageIntervals,
                               final SAMSequenceDictionary dictionary,
                               final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromReadDepthArgumentCollection arguments) {
        Utils.nonNull(breakpoints, "Breakpoint collection cannot be null");
        Utils.nonNull(assembledContigs, "Contig collection cannot be null");
        Utils.nonNull(evidenceTargetLinks, "Evidence target link collection cannot be null");
        Utils.nonNull(copyRatios, "Copy ratio collection cannot be null");
        Utils.nonNull(copyRatioSegments, "Copy ratio segments collection cannot be null");
        Utils.nonNull(dictionary, "Dictionary cannot be null");
        Utils.nonNull(arguments, "Parameter arguments collection cannot be null");
        Utils.nonNull(highCoverageIntervals, "High coverage intervals list cannot be null");

        this.dictionary = dictionary;
        this.arguments = arguments;

        logger.info("Building interval trees...");

        pairedBreakpoints = getIntrachromosomalBreakpointPairs(breakpoints);
        structuralVariantCallTree = buildVariantIntervalTree(structuralVariantCalls);
        highCoverageIntervalTree = buildSVIntervalTree(highCoverageIntervals);

        final Collection<EvidenceTargetLink> filteredEvidenceTargetLinks = new ArrayList<>(evidenceTargetLinks);

        final Collection<EvidenceTargetLink> intrachromosomalEvidenceTargetLinks = getIntrachromosomalLinks(filteredEvidenceTargetLinks);
        final Collection<EvidenceTargetLink> interchromosomalEvidenceTargetLinks = getInterchromosomalLinks(filteredEvidenceTargetLinks);
        intrachromosomalLinkTree = buildEvidenceIntervalTree(intrachromosomalEvidenceTargetLinks, 0, false);
        interchromosomalLinkTree = buildEvidenceIntervalTree(interchromosomalEvidenceTargetLinks, 0, true);

        contigTree = buildReadIntervalTree(assembledContigs);
        /* copyRatioOverlapDetector =  getMinimalCopyRatioCollection(copyRatios, copyRatios.getMetadata(),
                filteredEvidenceTargetLinks, pairedBreakpoints, dictionary, arguments.minEventSize, MAX_COPY_RATIO_EVENT_SIZE,
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

        logger.info("Initializing event factories...");
        tandemDuplicationFactory = new LargeTandemDuplicationFactory(intrachromosomalLinkTree, interchromosomalLinkTree,
                structuralVariantCallTree, contigTree, arguments, dictionary);
        largeDeletionFactory = new LargeDeletionFactory(intrachromosomalLinkTree, interchromosomalLinkTree,
                structuralVariantCallTree, contigTree, arguments, dictionary);
        dispersedDuplicationFactory = new DispersedDuplicationFactory(intrachromosomalLinkTree, interchromosomalLinkTree,
                structuralVariantCallTree, contigTree, arguments, dictionary);
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
     * Returns only links that are on the same chromosome
     */
    private static Collection<EvidenceTargetLink> getIntrachromosomalLinks(final Collection<EvidenceTargetLink> links) {
        return links.stream().filter(link -> link.getPairedStrandedIntervals().getLeft().getInterval().getContig() ==
                link.getPairedStrandedIntervals().getRight().getInterval().getContig()).collect(Collectors.toList());
    }

    /**
     * Returns only links that are on different chromosomes
     */
    private static Collection<EvidenceTargetLink> getInterchromosomalLinks(final Collection<EvidenceTargetLink> links) {
        return links.stream().filter(link -> link.getPairedStrandedIntervals().getLeft().getInterval().getContig() !=
                link.getPairedStrandedIntervals().getRight().getInterval().getContig()).collect(Collectors.toList());
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
     * Returns paired breakpoints on the same chromosome
     */
    private Collection<IntrachromosomalBreakpointPair> getIntrachromosomalBreakpointPairs(final Collection<VariantContext> breakpoints) {
        final Map<String, VariantContext> unpairedVariants = new HashMap<>();
        final Collection<IntrachromosomalBreakpointPair> pairedBreakpoints = new ArrayList<>(breakpoints.size() / 2);
        final Iterator<VariantContext> breakpointIter = breakpoints.iterator();
        while (breakpointIter.hasNext()) {
            final VariantContext vc1 = breakpointIter.next();
            if (!vc1.hasAttribute(GATKSVVCFConstants.BND_MATEID_STR)) continue;
            final String mate = vc1.getAttributeAsString(GATKSVVCFConstants.BND_MATEID_STR, "");
            if (unpairedVariants.containsKey(mate)) {
                final VariantContext vc2 = unpairedVariants.remove(mate);
                if (isBreakpointPair(vc1, vc2)) {
                    if (vc1.getContig().equals(vc2.getContig())) {
                        final int contig = dictionary.getSequenceIndex(vc1.getContig());
                        final VariantContext first;
                        final VariantContext second;
                        if (vc1.getStart() < vc2.getStart()) {
                            first = vc1;
                            second = vc2;
                        } else {
                            first = vc2;
                            second = vc1;
                        }
                        final int start = first.getStart();
                        final int end = second.getStart();
                        final Collection<String> firstContigs = first.getAttributeAsStringList(GATKSVVCFConstants.CONTIG_NAMES, "");
                        final Collection<String> secondContigs = second.getAttributeAsStringList(GATKSVVCFConstants.CONTIG_NAMES, "");
                        pairedBreakpoints.add(new IntrachromosomalBreakpointPair(contig, start, end, firstContigs, secondContigs));
                    }
                } else {
                    throw new IllegalStateException("Variant mate attributes did not match: " + vc1 + "\t" + vc2);
                }
            } else {
                unpairedVariants.put(vc1.getID(), vc1);
            }
        }
        if (!unpairedVariants.isEmpty()) {
            logger.warn("There were " + unpairedVariants.size() + " unpaired breakpoint variants with a " + GATKSVVCFConstants.BND_MATEID_STR + " attribute.");
        }
        return pairedBreakpoints;
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
    public Tuple2<Collection<ReadDepthEvent>, List<VariantContext>> callEvents(final JavaSparkContext ctx, final ProgressMeter progressMeter) {

        if (progressMeter != null) {
            progressMeter.setRecordLabel("intervals");
        }

        //In order to prevent duplicate calls, keep calls in a tree and filter subsequent intervals with sufficient overlap
        final SVIntervalTree<LargeSimpleSV> calledEventTree = new SVIntervalTree<>();

        //Filter unsupported existing calls
        for (final SVIntervalTree.Entry<VariantContext> entry : structuralVariantCallTree) {
            final VariantContext variantContext = entry.getValue();
            final SVInterval interval = entry.getInterval();
            if (variantContext.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "").equals(GATKSVVCFConstants.SYMB_ALT_ALLELE_DEL)
                    && interval.getLength() >= arguments.minEventSize) {
                final SVInterval leftInterval = new SVInterval(interval.getContig(), interval.getStart(), interval.getStart());
                final SVInterval rightInterval = new SVInterval(interval.getContig(), interval.getEnd(), interval.getEnd());
                final Collection<LargeSimpleSV> events = getEventsOnInterval(leftInterval, rightInterval, interval, null, arguments.breakpointPadding);
                for (final LargeSimpleSV event : events) {
                    if (event.getType() == SimpleSVType.TYPES.DEL) {
                        calledEventTree.put(event.getInterval(), event);
                    }
                }
                if (progressMeter != null) {
                    progressMeter.update(SVIntervalUtils.convertToSimpleInterval(interval, dictionary));
                }
            }
        }

        //Search breakpoint pairs
        for (final IntrachromosomalBreakpointPair breakpoints : pairedBreakpoints) {
            final SVInterval leftInterval = new SVInterval(breakpoints.getContig(), breakpoints.getInterval().getStart(), breakpoints.getInterval().getStart());
            final SVInterval rightInterval = new SVInterval(breakpoints.getContig(), breakpoints.getInterval().getEnd(), breakpoints.getInterval().getEnd());
            final Collection<LargeSimpleSV> events = getEventsOnInterval(leftInterval, rightInterval, breakpoints.getInterval(), breakpoints, arguments.breakpointPadding);
            for (final LargeSimpleSV event : events) {
                if (event.getType() == SimpleSVType.TYPES.DUP_TAND) {
                    calledEventTree.put(event.getInterval(), event);
                }
            }
            if (progressMeter != null) {
                progressMeter.update(SVIntervalUtils.convertToSimpleInterval(breakpoints.getInterval(), dictionary));
            }
        }

        //Search links
        final Iterator<SVIntervalTree.Entry<EvidenceTargetLink>> linkIter = intrachromosomalLinkTree.iterator();
        while (linkIter.hasNext()) {
            final EvidenceTargetLink link = linkIter.next().getValue();
            final SVInterval leftInterval = link.getPairedStrandedIntervals().getLeft().getInterval();
            final SVInterval rightInterval = link.getPairedStrandedIntervals().getRight().getInterval();
            if (leftInterval.getEnd() < rightInterval.getStart()) {
                final SVInterval callInterval = new SVInterval(leftInterval.getContig(), leftInterval.getEnd(), rightInterval.getStart());
                final Collection<LargeSimpleSV> events = getEventsOnInterval(leftInterval, rightInterval, callInterval, null, arguments.evidenceTargetLinkPadding);
                for (final LargeSimpleSV event : events) {
                    if (event.getType() == SimpleSVType.TYPES.DUP_TAND) {
                        calledEventTree.put(event.getInterval(), event);
                    }
                }
                if (progressMeter != null) {
                    progressMeter.update(SVIntervalUtils.convertToSimpleInterval(callInterval, dictionary));
                }
            }
        }


        final Set<LargeSimpleSV> callsToRemove = new HashSet<>();
        final Iterator<SVIntervalTree.Entry<LargeSimpleSV>> iterator = calledEventTree.iterator();
        while (iterator.hasNext()) {
            final SVIntervalTree.Entry<LargeSimpleSV> entry = iterator.next();
            final SVInterval interval = entry.getInterval();
            if (interval.overlaps(new SVInterval(2, 37710188, 37712149))) {
                int x = 0;
            }
            final LargeSimpleSV largeSimpleSV = entry.getValue();
            final Set<EvidenceTargetLink> supportingEvidence = new HashSet<>(largeSimpleSV.getSupportingEvidence());
            final List<LargeSimpleSV> conflictingCalls = Utils.stream(calledEventTree.overlappers(interval)).map(SVIntervalTree.Entry::getValue)
                    .filter(overlapper -> overlapper.getSupportingEvidence().stream().anyMatch(supportingEvidence::contains))
                    .filter(overlapper -> !callsToRemove.contains(overlapper))
                    .collect(Collectors.toList());
            if (conflictingCalls.size() < 2) continue;
            final double maxScore = conflictingCalls.stream().mapToDouble(call -> call.getReadPairEvidence()).max().getAsDouble();
            boolean bFoundMax = false;
            for (final LargeSimpleSV conflictingCall : conflictingCalls) {
                final double conflictingCallScore = conflictingCall.getReadPairEvidence();
                if (conflictingCallScore < maxScore) {
                    callsToRemove.add(conflictingCall);
                } else if (conflictingCallScore == maxScore) {
                    if (!bFoundMax) {
                        bFoundMax = true;
                    } else {
                        callsToRemove.add(conflictingCall);
                    }
                }
            }
        }

        final SVIntervalTree<LargeSimpleSV> filteredCalls = new SVIntervalTree<>();
        Utils.stream(calledEventTree.iterator())
                .filter(entry -> entry.getInterval().getContig() == 0)  //TODO
                .filter(entry -> !callsToRemove.contains(entry.getValue()))
                .filter(entry -> {
                    final Set<CalledCopyRatioSegment> overlappingSegments = copyRatioSegmentOverlapDetector.getOverlaps(SVIntervalUtils.convertToSimpleInterval(entry.getInterval(), dictionary));
                    final CalledCopyRatioSegment.Call expectedCall = entry.getValue().getType() == SimpleSVType.TYPES.DEL ? CalledCopyRatioSegment.Call.DELETION : CalledCopyRatioSegment.Call.AMPLIFICATION;
                    final int callOverlap = overlappingSegments.stream().filter(segment -> segment.getCall() == expectedCall).mapToInt(segment -> SVIntervalUtils.convertToSVInterval(segment.getInterval(), dictionary).overlapLen(entry.getInterval())).sum();
                    final double fractionOverlap = callOverlap / (double) entry.getInterval().getLength();
                    return fractionOverlap > 0.5;
                }).forEach(event -> filteredCalls.put(event.getInterval(), event.getValue()));

        logger.info("Running read depth model on " + filteredCalls.size() + " events");
        final ReadDepthModel readDepthModel = new ReadDepthModel(filteredCalls, copyRatioSegmentOverlapDetector, dictionary);
        final Tuple2<Double,List<ReadDepthEvent>> result = readDepthModel.solve(ctx);
        final Collection<ReadDepthEvent> finalResult = result._2;
        return new Tuple2<>(finalResult, Collections.emptyList());
    }

/*
    private final static class GradientDescentSolver {
        public final double[] gradientR;
        public final double[] deltaR;
        public final double[] mtR;
        public final double[] vtR;
        public final double[] r; //prob exists

        public final double[] gradientQ;
        public final double[] deltaQ;
        public final double[] mtQ;
        public final double[] vtQ;
        public final double[] q; //zygosity

        private final ReadDepthModelParameters modelParameters;

        public GradientDescentSolver(final ReadDepthModelParameters modelParameters, final int size) {
            this.modelParameters = modelParameters;
            gradientQ = new double[size];
            deltaQ = new double[size];
            mtQ = new double[size];
            vtQ = new double[size];
            q = new double[size];
            Arrays.fill(q, 1e-3);
            gradientR = new double[size];
            deltaR = new double[size];
            mtR = new double[size];
            vtR = new double[size];
            r = new double[size];
            Arrays.fill(r, 1e-3);
        }

        public void initialize(final double[] initR, final double[] initQ) {
            for (int i = 0; i < r.length; i++) {
                r[i] = initR[i];
            }
            for (int i = 0; i < q.length; i++) {
                q[i] = initQ[i];
            }
        }

        public double computeGradient(final List<LargeSimpleSV> eventsList,
                                    final List<Tuple2<Integer,Integer>> nearestCallDistances,
                                    final ArrayList<Tuple2<List<OverlapInfo>,Double>> copyNumberInfo) {
            return computeGradient(q, gradientQ, mtQ, eventsList, nearestCallDistances, copyNumberInfo);
                //+ computeGradient(r, gradientR, mtR, eventsList, nearestCallDistances, copyNumberInfo);
        }

        private double computeGradient(final double[] x,
                                     final double[] gradientX,
                                     final double[] mtX,
                                     final List<LargeSimpleSV> eventsList,
                                     final List<Tuple2<Integer,Integer>> nearestCallDistances,
                                     final ArrayList<Tuple2<List<OverlapInfo>,Double>> copyNumberInfo) {
            double total = 0;
            for (int i = 0; i < x.length; i++) {
                double xTemp = x[i];
                x[i] += modelParameters.gradientDelta; //mtX[i]*0.9 +
                final double testLogPosterior1 = computeLogPosterior(eventsList, nearestCallDistances, copyNumberInfo);
                x[i] = xTemp - modelParameters.gradientDelta; //+ mtX[i]*0.9
                final double testLogPosterior2 = computeLogPosterior(eventsList, nearestCallDistances, copyNumberInfo);
                x[i] = xTemp;
                gradientX[i] = (testLogPosterior1 - testLogPosterior2) / (2 * modelParameters.gradientDelta);
                total += gradientX[i]*gradientX[i];
            }
            return total;
        }

        public void computeDelta(final int numIter) {
            //computeDelta(gradientR, mtR, vtR, deltaR, numIter);
            computeDelta(gradientQ, mtQ, vtQ, deltaQ, numIter);
        }

        private void computeDelta(final double[] gradientX, final double[] mtX, final double[] vtX, final double[] deltaX, final int numIter) {
            for (int i = 0; i < gradientX.length; i++) {

                final double beta1 = 0.9; //0.99;
                final double beta2 = 0.999; //0.9999;
                if (numIter == 0) {
                    //mtX[i] = gradientX[i];
                    //vtX[i] = gradientX[i] * gradientX[i];
                }

                //Basic
                //deltaX[i] = gradientX[i] * learningRate;
                //if (Math.abs(deltaX[i]) > maxStepSize) {
                //    deltaX[i] = maxStepSize * Math.signum(deltaX[i]);
                //}

                //Simple momentum
                //deltaX[i] = mtX[i] * 0.99 + gradientX[i] * learningRate;
                //if (Math.abs(deltaX[i]) > maxStepSize) {
                //    deltaX[i] = maxStepSize * Math.signum(deltaX[i]);
                //}
                //mtX[i] = deltaX[i];

                //ADAM
                mtX[i] = (beta1 * mtX[i] +  (1 - beta1) * gradientX[i]); // / (1 - Math.pow(beta1, numIter + 1));
                vtX[i] = (beta2 * vtX[i] +  (1 - beta2) * gradientX[i] * gradientX[i]); // / (1 - Math.pow(beta2, numIter + 1));
                deltaX[i] = mtX[i] * modelParameters.learningRate / (Math.sqrt(vtX[i]) + 1e-8);

                //ADAMAX
                //mtX[i] = (beta1 * mtX[i] +  (1 - beta1) * gradientX[i]) / (1 - Math.pow(beta1, numIter + 1));
                //vtX[i] = Math.max(beta2 * vtX[i], Math.abs(gradientX[i]));
                //deltaX[i] = mtX[i] * learningRate / (vtX[i] + 1e-8);

            }
        }
        public double step() {
            return Math.sqrt(step(r, deltaR) + step(q, deltaQ));
        }

        private double step(final double[] x, final double[] deltaX) {
            double deltaRiTotal = 0;
            for (int i = 0; i < x.length; i++) {
                x[i] = x[i] + deltaX[i];
                deltaRiTotal += deltaX[i]*deltaX[i];
            }
            return deltaRiTotal;
        }

        public double computeLogPosterior(final List<LargeSimpleSV> eventsList,
                                        final List<Tuple2<Integer,Integer>> nearestCallDistances,
                                        final ArrayList<Tuple2<List<OverlapInfo>,Double>> copyNumberInfo) {
            return LargeSimpleSVCaller.computeLogPosterior(r, q, modelParameters, copyNumberInfo, eventsList, nearestCallDistances);
        }
    }
*/


    private Collection<LargeSimpleSV> testReadDepth(final SVIntervalTree<LargeSimpleSV> callTree, final ProgressMeter progressMeter) {
        logger.info("Evaluating calls for read depth...");
        final Collection<LargeSimpleSV> calledEvents = new ArrayList<>();
        for (int i = 0; i < dictionary.size(); i++) {
            calledEvents.addAll(testReadDepthOnContig(callTree, i, progressMeter));
        }
        return calledEvents;
    }

    private Collection<LargeSimpleSV> testReadDepthOnContig(final SVIntervalTree<LargeSimpleSV> callTree, final int contig, final ProgressMeter progressMeter) {
/*        final Set<LargeSimpleSV> visited = new HashSet<>(SVUtils.hashMapCapacity(callTree.size()));

        final Set<LargeSimpleSV> eventsToFilter = new HashSet<>(SVUtils.hashMapCapacity(callTree.size()));
        final List<LargeSimpleSV> callList = Utils.stream(callTree.iterator()).map(SVIntervalTree.Entry::getValue).collect(Collectors.toList());
        Collections.sort(callList, Comparator.comparingInt(event -> -event.getSize()));
        for (final LargeSimpleSV call : callList) {
            if (!eventsToFilter.contains(call)) {
                final List<LargeSimpleSV> overlappers = Utils.stream(callTree.overlappers(call.getInterval())).map(SVIntervalTree.Entry::getValue)
                        .filter(event -> !eventsToFilter.contains(event))
                        .collect(Collectors.toList());
                if (overlappers.size() > 1) {
                    eventsToFilter.add(call);
                }
            }
        }
*/

        final List<LargeSimpleSV> filteredCalls = Utils.stream(callTree.iterator()).filter(entry -> entry.getInterval().getContig() == contig)
                .map(SVIntervalTree.Entry::getValue).collect(Collectors.toList()); // callList.stream().filter(call -> !eventsToFilter.contains(call)).collect(Collectors.toList());
        final SVIntervalTree<LargeSimpleSV> filteredCallTree = callTree; // new SVIntervalTree<>();
        /*
        for (final LargeSimpleSV event : filteredCalls) {
            filteredCallTree.put(event.getInterval(), event);
        }
        */
        final SAMSequenceRecord contigRecord = dictionary.getSequence(contig);
        if (contigRecord == null) {
            throw new IllegalArgumentException("Could not find contig with index " + contig + " in dictionary");
        }
        final SVIntervalTree<Float> copyRatioTree = new SVIntervalTree<>();
        for (final SVCopyRatio copyRatio : copyRatios.get(contig)) {
            copyRatioTree.put(copyRatio.getInterval(), copyRatio.getLog2CopyRatio());
        }

        final Map<SimpleSVType.TYPES, CalledCopyRatioSegment.Call> expectedCallMap = new HashMap<>(SVUtils.hashMapCapacity(2));
        expectedCallMap.put(SimpleSVType.TYPES.DEL, CalledCopyRatioSegment.Call.DELETION);
        expectedCallMap.put(SimpleSVType.TYPES.DUP_TAND, CalledCopyRatioSegment.Call.AMPLIFICATION);

        final Map<SimpleSVType.TYPES, Set<Integer>> validHMMStatesMap = new HashMap<>(SVUtils.hashMapCapacity(2));
        validHMMStatesMap.put(SimpleSVType.TYPES.DEL, IntStream.range(0, 2).boxed().collect(Collectors.toSet()));
        validHMMStatesMap.put(SimpleSVType.TYPES.DUP_TAND, IntStream.range(3, arguments.hmmMaxStates).boxed().collect(Collectors.toSet()));

        final Collection<LargeSimpleSV> supportedCalls = new ArrayList<>();
        for (final LargeSimpleSV largeSimpleSV : filteredCalls) {
            //if (!visited.contains(largeSimpleSV)) {
                /*
                SVInterval setInterval = largeSimpleSV.getInterval();
                final Set<LargeSimpleSV> overlappers = new HashSet<>();
                overlappers.add(largeSimpleSV);
                int lastSize;
                do {
                    lastSize = overlappers.size();
                    final SVInterval finalInterval = setInterval;
                    overlappers.addAll(Utils.stream(filteredCallTree.overlappers(setInterval))
                            .filter(overlapper -> SVIntervalUtils.hasReciprocalOverlap(finalInterval, overlapper.getInterval(), 0.2))
                            .map(SVIntervalTree.Entry::getValue)
                            .collect(Collectors.toList()));
                    final int newStart = overlappers.stream().mapToInt(LargeSimpleSV::getStart).min().getAsInt();
                    final int newEnd = overlappers.stream().mapToInt(LargeSimpleSV::getEnd).max().getAsInt();
                    setInterval = new SVInterval(setInterval.getContig(), newStart, newEnd);
                } while (lastSize < overlappers.size());
                visited.addAll(overlappers);

                if (overlappers.size() > 1) continue;
                if (!expectedCalls.containsKey(largeSimpleSV.getType())) continue;
                */

            final SVInterval eventSVInterval = largeSimpleSV.getInterval();

            final SVInterval leftEnd = new SVInterval(eventSVInterval.getContig(), eventSVInterval.getStart(), eventSVInterval.getStart());
            final SVInterval rightEnd = new SVInterval(eventSVInterval.getContig(), eventSVInterval.getEnd(), eventSVInterval.getEnd());
            if (highCoverageIntervalTree.hasOverlapper(leftEnd) || highCoverageIntervalTree.hasOverlapper(rightEnd)) {
                continue;
            }
            final int highCoverageLength = Utils.stream(highCoverageIntervalTree.overlappers(eventSVInterval))
                    .mapToInt(entry -> entry.getInterval().overlapLen(eventSVInterval)).sum();
            if (highCoverageLength > 0.5 * eventSVInterval.getLength()) {
                continue;
            }

            final SimpleInterval eventSimpleInterval = SVIntervalUtils.convertToSimpleInterval(eventSVInterval, dictionary);
            final SimpleSVType.TYPES eventType = largeSimpleSV.getType();

            final List<CalledCopyRatioSegment> overlappingSegments = new ArrayList<>(copyRatioSegmentOverlapDetector.getOverlaps(eventSimpleInterval));
            /*
            Collections.sort(overlappingSegments, IntervalUtils.getDictionaryOrderComparator(dictionary));
            List<CalledCopyRatioSegment> mergedOverlappingSegments = new ArrayList<>(overlappingSegments.size());
            CalledCopyRatioSegment lastSegment = null;
            for (final CalledCopyRatioSegment segment : overlappingSegments) {
                if (lastSegment == null) {
                    lastSegment = segment;
                } else if (segment.getCall() == lastSegment.getCall() && lastSegment.getEnd() + 1 == segment.getStart()) {
                    final SimpleInterval mergedInterval = new SimpleInterval(segment.getContig(), lastSegment.getStart(), segment.getEnd());
                    final CopyRatioSegment mergedSegment = new CopyRatioSegment(mergedInterval, lastSegment.getNumPoints() + segment.getNumPoints(), 0);
                    lastSegment = new CalledCopyRatioSegment(mergedSegment, lastSegment.getCall());
                } else {
                    mergedOverlappingSegments.add(lastSegment);
                    lastSegment = segment;
                }
            }
            if (lastSegment != null) {
                mergedOverlappingSegments.add(lastSegment);
            }
            */
            /*
            if (mergedOverlappingSegments.size() > 1) {
                int lastNumIntervals;
                List<CalledCopyRatioSegment> newMergedSegments;
                do {
                    newMergedSegments = new ArrayList<>(mergedOverlappingSegments.size());
                    for (int i = 0; i < mergedOverlappingSegments.size() - 1; i++) {
                        final CalledCopyRatioSegment a = mergedOverlappingSegments.get(i);
                        final CalledCopyRatioSegment b = mergedOverlappingSegments.get(i + 1);
                        if (a.getEnd() + 1 == b.getStart() && a.getCall() == b.getCall()) {
                            final SimpleInterval mergedInterval = new SimpleInterval(a.getContig(), a.getStart(), b.getEnd());
                            final CopyRatioSegment mergedSegment = new CopyRatioSegment(mergedInterval, a.getNumPoints() + b.getNumPoints(), 0);
                            newMergedSegments.add(new CalledCopyRatioSegment(mergedSegment, a.getCall()));
                            i++;
                        } else {
                            newMergedSegments.add(a);
                            if (i == mergedOverlappingSegments.size() - 2) {
                                newMergedSegments.add(b);
                            }
                        }
                    }
                    lastNumIntervals = mergedOverlappingSegments.size();
                    mergedOverlappingSegments = newMergedSegments;
                } while (lastNumIntervals > mergedOverlappingSegments.size());
            }*/

            final Collection<CalledCopyRatioSegment> supportingSegments = overlappingSegments.stream().filter(segment -> segment.getCall() == expectedCallMap.get(eventType)
                    && SVIntervalUtils.hasReciprocalOverlap(SVIntervalUtils.convertToSVInterval(segment.getInterval(), dictionary), eventSVInterval, 0.5)
                    && Math.abs(segment.getInterval().getStart() - eventSVInterval.getStart()) < 10000
                    && Math.abs(segment.getInterval().getEnd() - eventSVInterval.getEnd()) < 10000)
                    .collect(Collectors.toList());
            if (expectedCallMap.containsKey(eventType) && !supportingSegments.isEmpty()) {
                //largeSimpleSV.setReadDepthSupportType("SEGMENTS");
                supportedCalls.add(largeSimpleSV);
            } else if (eventSVInterval.getLength() < 10000) {
                final List<SVCopyRatio> eventBins = getCopyRatiosOnInterval(eventSVInterval, copyRatioTree, arguments.copyRatioBinTrimming, dictionary);
                if (validHMMStatesMap.containsKey(eventType) && supportedByHMM(eventBins, validHMMStatesMap.get(eventType))) {
                    //largeSimpleSV.setReadDepthSupportType("HMM");
                    supportedCalls.add(largeSimpleSV);
                }
            } else if (eventType == SimpleSVType.TYPES.DEL && largeSimpleSV.getSupportingEvidence().size() >= 3) {
                final List<SVCopyRatio> eventBins = getCopyRatiosOnInterval(eventSVInterval, copyRatioTree, arguments.copyRatioBinTrimming, dictionary);
                if (fractionEmpty(eventSVInterval, eventBins) > 0.8) {
                    //largeSimpleSV.setReadDepthSupportType("DEL_RESCUE");
                    supportedCalls.add(largeSimpleSV);
                }
            }
            //}
            if (progressMeter != null) {
                progressMeter.setRecordLabel(largeSimpleSV.getInterval() + ":" + largeSimpleSV.getType());
            }
        }
        return supportedCalls;
    }

    /**
     * Gets the event with the highest score on the interval defned by leftInterval and rightInterval. For evidence
     * target links to count as evidence, their paired intervals must overlap leftInterval and rightInterval.
     *
     * @param leftInterval        The left interval
     * @param rightInterval       The right interval
     * @param callInterval        Calls will be made using this interval
     * @param disallowedIntervals No calls will be made if they sufficiently overlap intervals in this tree
     * @param breakpoints         Breakpoints associated with the interval (may be null)
     * @return Optional containing the called event, if any
     */
    private Optional<LargeSimpleSV> getHighestScoringEventOnInterval(final SVInterval leftInterval, final SVInterval rightInterval, final SVInterval callInterval, final SVIntervalTree<LargeSimpleSV> disallowedIntervals, final IntrachromosomalBreakpointPair breakpoints, final int evidencePadding) {
        if (SVIntervalUtils.hasReciprocalOverlapInTree(callInterval, disallowedIntervals, arguments.maxCallReciprocalOverlap))
            return Optional.empty();
        final Stream<LargeSimpleSV> candidateEvents = getEventsOnInterval(leftInterval, rightInterval, callInterval, breakpoints, evidencePadding).stream();
        return candidateEvents.max(Comparator.comparingDouble(event -> event.getScore(arguments.counterEvidencePseudocount)));
    }

    /**
     * Gets collection of valid events on the interval. Evidence target link stranded intervals must overlap leftInterval and rightInterval.
     *
     * @param leftInterval  The left interval
     * @param rightInterval The right interval
     * @param callInterval  The interval to use for the call
     * @param breakpoints   Associated breakpoints (may be null)
     * @return Collection of called events
     */
    private Collection<LargeSimpleSV> getEventsOnInterval(final SVInterval leftInterval,
                                                          final SVInterval rightInterval,
                                                          final SVInterval callInterval,
                                                          final IntrachromosomalBreakpointPair breakpoints,
                                                          final int evidencePadding) {

        if (leftInterval.getContig() != rightInterval.getContig()) return Collections.emptyList();
        if (callInterval.getLength() < arguments.minEventSize)
            return Collections.emptyList();

        final Collection<LargeSimpleSV> events = new ArrayList<>();

        final LargeSimpleSV tandemDuplication = tandemDuplicationFactory.call(leftInterval, rightInterval, callInterval, breakpoints, evidencePadding);
        if (tandemDuplication != null) events.add(tandemDuplication);

        final LargeSimpleSV deletion = largeDeletionFactory.call(leftInterval, rightInterval, callInterval, breakpoints, evidencePadding);
        if (deletion != null) events.add(deletion);

        return events;
    }

    /**
     * Tests whether the copy ratios support an event with an HMM
     *
     * @param copyRatioBins Copy ratios over the interval
     * @return True if supported
     */
    protected boolean supportedByHMM(final List<SVCopyRatio> copyRatioBins, final Set<Integer> validStates) {
        if (copyRatioBins.isEmpty()) {
            return false;
        }
        final List<Double> copyRatios = copyRatioBins.stream().map(SVCopyRatio::getLog2CopyRatio).map(val -> Math.pow(2.0, val)).collect(Collectors.toList());
        final int numStates = Math.min(arguments.hmmMaxStates, copyRatios.stream().mapToInt(val -> (int) (2 * val)).max().getAsInt() + 1);
        final RealVector copyNumberPriors = CopyNumberHMM.uniformPrior(numStates);
        final List<Integer> positionsList = CopyNumberHMM.positionsList(copyRatios.size());
        final CopyNumberHMM copyNumberHMM = new CopyNumberHMM(copyNumberPriors, arguments.hmmTransitionProb);
        final List<Integer> copyNumberStates = ViterbiAlgorithm.apply(copyRatios, positionsList, copyNumberHMM);
        return testHMMState(copyNumberStates, validStates);
    }

    /**
     * Tests if the HMM state path contains a sufficient proportion of valid states
     *
     * @param states      State path
     * @param validStates Valid HMM states
     * @return True if the threshold is met
     */
    private boolean testHMMState(final List<Integer> states, final Set<Integer> validStates) {
        return validStateFrequency(states, validStates) >= arguments.hmmValidStatesMinFraction * states.size();
    }

    /**
     * Counts the number of states in the path that are one of the valid states
     */
    private int validStateFrequency(final List<Integer> states, final Set<Integer> validStates) {
        return (int) states.stream().filter(validStates::contains).count();
    }

}
