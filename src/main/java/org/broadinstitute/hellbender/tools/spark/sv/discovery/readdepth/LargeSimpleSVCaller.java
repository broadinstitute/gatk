package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ProgressMeter;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledCopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.LargeSimpleSV;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.Utils;
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

    private static final int MAX_COPY_RATIO_EVENT_SIZE = 100000;

    private final OverlapDetector<CopyRatio> copyRatioOverlapDetector;
    private final SVIntervalTree<VariantContext> structuralVariantCallTree;
    private final SVIntervalTree<EvidenceTargetLink> intrachromosomalLinkTree;
    private final SVIntervalTree<EvidenceTargetLink> interchromosomalLinkTree;
    private final SVIntervalTree<GATKRead> contigTree;
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
                               final SAMSequenceDictionary dictionary,
                               final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromReadDepthArgumentCollection arguments) {
        Utils.nonNull(breakpoints, "Breakpoint collection cannot be null");
        Utils.nonNull(assembledContigs, "Contig collection cannot be null");
        Utils.nonNull(evidenceTargetLinks, "Evidence target link collection cannot be null");
        Utils.nonNull(copyRatios, "Copy ratio collection cannot be null");
        Utils.nonNull(copyRatioSegments, "Copy ratio segments collection cannot be null");
        Utils.nonNull(dictionary, "Dictionary cannot be null");
        Utils.nonNull(arguments, "Parameter arguments collection cannot be null");

        this.dictionary = dictionary;
        this.arguments = arguments;

        logger.info("Building interval trees...");

        pairedBreakpoints = getIntrachromosomalBreakpointPairs(breakpoints);
        structuralVariantCallTree = buildVariantIntervalTree(structuralVariantCalls);

        final Collection<EvidenceTargetLink> filteredEvidenceTargetLinks = new ArrayList<>(evidenceTargetLinks);

        final Collection<EvidenceTargetLink> intrachromosomalEvidenceTargetLinks = getIntrachromosomalLinks(filteredEvidenceTargetLinks);
        final Collection<EvidenceTargetLink> interchromosomalEvidenceTargetLinks = getInterchromosomalLinks(filteredEvidenceTargetLinks);
        intrachromosomalLinkTree = buildEvidenceIntervalTree(intrachromosomalEvidenceTargetLinks, 0, false);
        interchromosomalLinkTree = buildEvidenceIntervalTree(interchromosomalEvidenceTargetLinks, 0, true);

        contigTree = buildReadIntervalTree(assembledContigs);
        copyRatioOverlapDetector = getMinimalCopyRatioCollection(copyRatios, copyRatios.getMetadata(),
                filteredEvidenceTargetLinks, pairedBreakpoints, dictionary, arguments.minEventSize, MAX_COPY_RATIO_EVENT_SIZE,
                arguments.breakpointPadding + arguments.hmmPadding,
                arguments.evidenceTargetLinkPadding + arguments.hmmPadding).getOverlapDetector();
        copyRatioSegmentOverlapDetector = copyRatioSegments.getOverlapDetector();

        logger.info("Initializing event factories...");
        tandemDuplicationFactory = new LargeTandemDuplicationFactory(intrachromosomalLinkTree, interchromosomalLinkTree,
                structuralVariantCallTree, contigTree, arguments, copyRatioSegmentOverlapDetector, copyRatioOverlapDetector, dictionary);
        largeDeletionFactory = new LargeDeletionFactory(intrachromosomalLinkTree, interchromosomalLinkTree,
                structuralVariantCallTree, contigTree, arguments, copyRatioSegmentOverlapDetector, copyRatioOverlapDetector, dictionary);
        dispersedDuplicationFactory = new DispersedDuplicationFactory(intrachromosomalLinkTree, interchromosomalLinkTree,
                structuralVariantCallTree, contigTree, arguments, copyRatioSegmentOverlapDetector, copyRatioOverlapDetector, dictionary);

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
        final SVIntervalTree mergedIntervalTree = new SVIntervalTree();
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
    public Tuple2<Collection<LargeSimpleSV>, List<VariantContext>> callEvents(final ProgressMeter progressMeter) {

        if (progressMeter != null) {
            progressMeter.setRecordLabel("intervals");
        }

        //In order to prevent duplicate calls, keep calls in a tree and filter subsequent intervals with sufficient overlap
        final SVIntervalTree<LargeSimpleSV> calledEventTree = new SVIntervalTree<>();

        //Filter unsupported existing calls
        final List<VariantContext> filteredCalls = new ArrayList<>();
        for (final SVIntervalTree.Entry<VariantContext> entry : structuralVariantCallTree) {
            final VariantContext variantContext = entry.getValue();
            final SVInterval interval = entry.getInterval();
            if (variantContext.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "").equals(GATKSVVCFConstants.SYMB_ALT_ALLELE_DEL)
                    && interval.getLength() >= arguments.minEventSize) {
                final SVInterval leftInterval = new SVInterval(interval.getContig(), interval.getStart(), interval.getStart());
                final SVInterval rightInterval = new SVInterval(interval.getContig(), interval.getEnd(), interval.getEnd());
                final Collection<LargeSimpleSV> events = getEventsOnInterval(leftInterval, rightInterval, interval, null, arguments.breakpointPadding);
                for (final LargeSimpleSV event : events) {
                    calledEventTree.put(event.getInterval(), event);
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
            //final Optional<LargeSimpleSV> event = getHighestScoringEventOnInterval(leftBreakpointInterval, rightBreakpointInterval, breakpoints.getInterval(), calledEventTree, breakpoints, arguments.breakpointPadding);
            final Collection<LargeSimpleSV> events = getEventsOnInterval(leftInterval, rightInterval, breakpoints.getInterval(), breakpoints, arguments.breakpointPadding);
            for (final LargeSimpleSV event : events) {
                calledEventTree.put(event.getInterval(), event);
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
            final int callEnd1 = (leftInterval.getStart() + leftInterval.getEnd()) / 2;
            final int callEnd2 = (rightInterval.getStart() + rightInterval.getEnd()) / 2;
            final SVInterval callInterval = new SVInterval(leftInterval.getContig(), Math.min(callEnd1, callEnd2), Math.max(callEnd1, callEnd2));
            //final Optional<LargeSimpleSV> event = getHighestScoringEventOnInterval(leftInterval, rightInterval, callInterval, calledEventTree, null, arguments.evidenceTargetLinkPadding);
            final Collection<LargeSimpleSV> events = getEventsOnInterval(leftInterval, rightInterval, callInterval, null, arguments.evidenceTargetLinkPadding);
            for (final LargeSimpleSV event : events) {
                calledEventTree.put(event.getInterval(), event);
            }
            if (progressMeter != null) {
                progressMeter.update(SVIntervalUtils.convertToSimpleInterval(callInterval, dictionary));
            }
        }

        final Set<LargeSimpleSV> callsToRemove = new HashSet<>();
        final Iterator<SVIntervalTree.Entry<LargeSimpleSV>> iterator = calledEventTree.iterator();
        while (iterator.hasNext()) {
            final SVIntervalTree.Entry<LargeSimpleSV> entry = iterator.next();
            final SVInterval interval = entry.getInterval();
            final LargeSimpleSV largeSimpleSV = entry.getValue();
            final Set<EvidenceTargetLink> supportingEvidence = new HashSet<>(largeSimpleSV.getSupportingEvidence());
            final List<LargeSimpleSV> conflictingCalls = Utils.stream(calledEventTree.overlappers(interval)).map(SVIntervalTree.Entry::getValue)
                    .filter(overlapper -> overlapper.getSupportingEvidence().stream().anyMatch(supportingEvidence::contains))
                    .filter(overlapper -> !callsToRemove.contains(overlapper))
                    .collect(Collectors.toList());
            if (conflictingCalls.size() < 2) continue;
            final double maxScore = conflictingCalls.stream().mapToDouble(call -> call.getScore(arguments.counterEvidencePseudocount)).max().getAsDouble();
            boolean bFoundMax = false;
            for (final LargeSimpleSV conflictingCall : conflictingCalls) {
                final double conflictingCallScore = conflictingCall.getScore(arguments.counterEvidencePseudocount);
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
        final SVIntervalTree<LargeSimpleSV> uniqueCalledEvents = new SVIntervalTree<>();
        final Iterator<SVIntervalTree.Entry<LargeSimpleSV>> iterator2 = calledEventTree.iterator();
        while (iterator2.hasNext()) {
            final SVIntervalTree.Entry<LargeSimpleSV> entry = iterator2.next();
            if (!callsToRemove.contains(entry.getValue())) {
                uniqueCalledEvents.put(entry.getInterval(), entry.getValue());
            }
        }

        return new Tuple2<>(testReadDepth(uniqueCalledEvents), filteredCalls);
    }

    private Collection<LargeSimpleSV> testReadDepth(final SVIntervalTree<LargeSimpleSV> callTree) {
        final Set<LargeSimpleSV> visited = new HashSet<>(SVUtils.hashMapCapacity(callTree.size()));
        final Iterator<SVIntervalTree.Entry<LargeSimpleSV>> iter = callTree.iterator();
        final Collection<LargeSimpleSV> supportedCalls = new ArrayList<>();
        while (iter.hasNext()) {
            final SVIntervalTree.Entry<LargeSimpleSV> entry = iter.next();
            final LargeSimpleSV largeSimpleSV = entry.getValue();
            if (!visited.contains(largeSimpleSV)) {
                SVInterval setInterval = entry.getInterval();
                final Set<LargeSimpleSV> overlappers = new HashSet<>();
                overlappers.add(largeSimpleSV);
                int lastSize;
                do {
                    lastSize = overlappers.size();
                    final SVInterval finalInterval = setInterval;
                    overlappers.addAll(Utils.stream(callTree.overlappers(setInterval))
                            .filter(overlapper -> SVIntervalUtils.hasReciprocalOverlap(finalInterval, overlapper.getInterval(), 0.2))
                            .map(SVIntervalTree.Entry::getValue)
                            .collect(Collectors.toList()));
                    final int newStart = overlappers.stream().mapToInt(LargeSimpleSV::getStart).min().getAsInt();
                    final int newEnd = overlappers.stream().mapToInt(LargeSimpleSV::getEnd).max().getAsInt();
                    setInterval = new SVInterval(setInterval.getContig(), newStart, newEnd);
                } while (lastSize < overlappers.size());
                visited.addAll(overlappers);
                if (setInterval.overlaps(new SVInterval(5, 118691685, 118691686))) {
                    int x = 0;
                }
                if (overlappers.size() > 2) continue;

                final List<LargeSimpleSV> overlappersList = new ArrayList<>(overlappers);
                double bestScore = 0;
                List<LargeSimpleSV> bestOverlapperCombination = Collections.emptyList();
                final int setSize = overlappers.size();
                {
                    final Iterator<int[]> combinationIter = CombinatoricsUtils.combinationsIterator(overlappersList.size(), setSize);
                    while (combinationIter.hasNext()) {
                        final int[] combinationIndices = combinationIter.next();
                        final List<LargeSimpleSV> overlapperCombination = new ArrayList<>(combinationIndices.length);
                        for (int k = 0; k < combinationIndices.length; k++) {
                            overlapperCombination.add(overlappersList.get(combinationIndices[k]));
                        }

                        final int newStart = overlapperCombination.stream().mapToInt(LargeSimpleSV::getStart).min().getAsInt();
                        final int newEnd = overlapperCombination.stream().mapToInt(LargeSimpleSV::getEnd).max().getAsInt();
                        setInterval = new SVInterval(setInterval.getContig(), newStart, newEnd);

                        final int totalNumDUP = (int) overlapperCombination.stream().filter(event -> event.getType() == SimpleSVType.TYPES.DUP_TAND).count();
                        final int totalNumDEL = (int) overlapperCombination.stream().filter(event -> event.getType() == SimpleSVType.TYPES.DEL).count();

                        final int binSize = 50;
                        final int numBins = 1 + (setInterval.getLength() / binSize);
                        int numMatchingBinCalls = 0;
                        final int intervalStart = setInterval.getStart();
                        for (int i = 0; i < numBins; i++) {
                            final SVInterval binInterval = new SVInterval(setInterval.getContig(), intervalStart + i * binSize, intervalStart + (i + 1) * binSize);

                            final List<LargeSimpleSV> binOverlappingEvents = overlapperCombination.stream().filter(event -> event.getInterval().overlaps(binInterval)).collect(Collectors.toList());
                            final int numDUP = (int) binOverlappingEvents.stream().filter(event -> event.getType() == SimpleSVType.TYPES.DUP_TAND).count();
                            final int numDEL = (int) binOverlappingEvents.stream().filter(event -> event.getType() == SimpleSVType.TYPES.DEL).count();

                            final Set<CalledCopyRatioSegment> segments = copyRatioSegmentOverlapDetector.getOverlaps(SVIntervalUtils.convertToSimpleInterval(binInterval, dictionary));
                            if (numDUP == 0 && numDEL == 1) { //1 deletion
                                for (final CalledCopyRatioSegment segment : segments) {
                                    if (segment.getCall() == CalledCopyRatioSegment.Call.DELETION
                                            || (segment.getCall() == CalledCopyRatioSegment.Call.NEUTRAL && totalNumDEL == 1 && totalNumDUP == 1)) {
                                        numMatchingBinCalls++;
                                    }
                                }
                            } else if (numDUP == 0 && numDEL == 2) { //2 deletions
                                for (final CalledCopyRatioSegment segment : segments) {
                                    if (segment.getCall() == CalledCopyRatioSegment.Call.DELETION) {
                                        numMatchingBinCalls++;
                                    }
                                }
                            } else if (numDUP == 1 && numDEL == 0) { //1 tandem duplication
                                for (final CalledCopyRatioSegment segment : segments) {
                                    if (segment.getCall() == CalledCopyRatioSegment.Call.AMPLIFICATION
                                            || (segment.getCall() == CalledCopyRatioSegment.Call.NEUTRAL && totalNumDEL == 1 && totalNumDUP == 1)) {
                                        numMatchingBinCalls++;
                                    }
                                }
                            } else if (numDUP == 2 && numDEL == 0) { //2 tandem duplications
                                for (final CalledCopyRatioSegment segment : segments) {
                                    if (segment.getCall() == CalledCopyRatioSegment.Call.AMPLIFICATION) {
                                        numMatchingBinCalls++;
                                    }
                                }
                            } else if (numDUP == 1 && numDEL == 1) { //1 dispersed duplication or 1 deletion/1 tandem duplication
                                for (final CalledCopyRatioSegment segment : segments) {
                                    if (segment.getCall() == CalledCopyRatioSegment.Call.NEUTRAL) {
                                        numMatchingBinCalls++;
                                    }
                                }
                            } else {
                                continue;
                            }
                        }
                        final double score = numMatchingBinCalls / (double) numBins;
                        if (score > 0.8) {
                            if (score * overlapperCombination.size() > bestScore) {
                                bestScore = score * overlapperCombination.size();
                                bestOverlapperCombination = overlapperCombination;
                            }
                        }
                    }
                }
                if (bestScore > 0) {
                    supportedCalls.addAll(bestOverlapperCombination);
                    if (bestOverlapperCombination.size() > 1) {
                        int x = 0;
                    }
                }
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

}
