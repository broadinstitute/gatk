package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.ValidationStringency;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.ProgressMeter;
import org.broadinstitute.hellbender.engine.TraversalParameters;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
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
public final class StructuralVariantIntervalFinder {

    private static final Logger logger = LogManager.getLogger(LargeSimpleSVCaller.class);

    private final SVIntervalTree<Object> highCoverageIntervalTree;
    private final SVIntervalTree<Object> mappableIntervalTree;
    private final SVIntervalTree<Object> blacklistTree;
    private final SVIntervalTree<VariantContext> structuralVariantCallTree;
    private final SVIntervalTree<EvidenceTargetLink> intrachromosomalLinkTree;
    private final SVIntervalTree<EvidenceTargetLink> interchromosomalLinkTree;
    private final SVIntervalTree<GATKRead> contigTree;
    private final LargeSimpleSVFactory tandemDuplicationFactory;
    private final LargeDeletionFactory largeDeletionFactory;
    private final SAMSequenceDictionary dictionary;
    private final Collection<IntrachromosomalBreakpointPair> pairedBreakpoints;
    private final StructuralVariationDiscoveryArgumentCollection.StructuralVariantIntervalsForCNV arguments;

    public StructuralVariantIntervalFinder(final Collection<VariantContext> breakpoints,
                               final Collection<VariantContext> structuralVariantCalls,
                               final Collection<GATKRead> assembledContigs,
                               final Collection<EvidenceTargetLink> evidenceTargetLinks,
                               final List<SVInterval> highCoverageIntervals,
                               final List<SVInterval> mappableIntervals,
                               final List<SVInterval> blacklist,
                               final SAMSequenceDictionary dictionary,
                               final StructuralVariationDiscoveryArgumentCollection.StructuralVariantIntervalsForCNV arguments) {
        Utils.nonNull(breakpoints, "Breakpoint collection cannot be null");
        Utils.nonNull(assembledContigs, "Contig collection cannot be null");
        Utils.nonNull(evidenceTargetLinks, "Evidence target link collection cannot be null");
        Utils.nonNull(dictionary, "Dictionary cannot be null");
        Utils.nonNull(arguments, "Parameter arguments collection cannot be null");
        Utils.nonNull(highCoverageIntervals, "High coverage intervals list cannot be null");

        this.dictionary = dictionary;
        this.arguments = arguments;

        logger.info("Building interval trees...");

        pairedBreakpoints = getIntrachromosomalBreakpointPairs(breakpoints);
        structuralVariantCallTree = buildVariantIntervalTree(structuralVariantCalls);
        highCoverageIntervalTree = buildSVIntervalTree(highCoverageIntervals);
        mappableIntervalTree = buildSVIntervalTree(mappableIntervals);
        blacklistTree = buildSVIntervalTree(blacklist);

        final Collection<EvidenceTargetLink> filteredEvidenceTargetLinks = new ArrayList<>(evidenceTargetLinks);

        final Collection<EvidenceTargetLink> intrachromosomalEvidenceTargetLinks = getIntrachromosomalLinks(filteredEvidenceTargetLinks);
        final Collection<EvidenceTargetLink> interchromosomalEvidenceTargetLinks = getInterchromosomalLinks(filteredEvidenceTargetLinks);
        intrachromosomalLinkTree = buildEvidenceIntervalTree(intrachromosomalEvidenceTargetLinks, 0, false);
        interchromosomalLinkTree = buildEvidenceIntervalTree(interchromosomalEvidenceTargetLinks, 0, true);

        contigTree = buildReadIntervalTree(assembledContigs);

        logger.info("Initializing event factories...");
        tandemDuplicationFactory = new LargeTandemDuplicationFactory(intrachromosomalLinkTree, interchromosomalLinkTree,
                structuralVariantCallTree, contigTree, arguments, dictionary);
        largeDeletionFactory = new LargeDeletionFactory(intrachromosomalLinkTree, interchromosomalLinkTree,
                structuralVariantCallTree, contigTree, arguments, dictionary);
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
        if (callInterval.getLength() < arguments.smallEventSize)
            return Collections.emptyList();

        final Collection<LargeSimpleSV> events = new ArrayList<>();

        final LargeSimpleSV tandemDuplication = tandemDuplicationFactory.call(leftInterval, rightInterval, callInterval, breakpoints, evidencePadding);
        if (tandemDuplication != null) events.add(tandemDuplication);

        final LargeSimpleSV deletion = largeDeletionFactory.call(leftInterval, rightInterval, callInterval, breakpoints, evidencePadding);
        if (deletion != null) events.add(deletion);

        return events;
    }

    /**
     * Returns all events. Searches by iterating over the breakpoint pairs and then the evidence target links.
     */
    public Map<String,List<SimpleInterval>> getIntervals(final ProgressMeter progressMeter) {

        if (progressMeter != null) {
            progressMeter.setRecordLabel("intervals");
        }

        final List<SVInterval> intervals = new ArrayList<>();

        final List<SVInterval> deletionCallIntervals = Utils.stream(structuralVariantCallTree.iterator())
                .filter(entry -> entry.getValue().getStructuralVariantType() == StructuralVariantType.DEL)
                .filter(entry -> entry.getInterval().getLength() >= arguments.smallEventSize)
                .map(SVIntervalTree.Entry::getInterval)
                .collect(Collectors.toList());
        intervals.addAll(deletionCallIntervals);

        //Search breakpoint pairs for tandem duplications
        final SVIntervalTree<LargeSimpleSV> duplicationEventTree = new SVIntervalTree<>();
        for (final IntrachromosomalBreakpointPair breakpoints : pairedBreakpoints) {
            final SVInterval leftInterval = new SVInterval(breakpoints.getContig(), breakpoints.getInterval().getStart(), breakpoints.getInterval().getStart());
            final SVInterval rightInterval = new SVInterval(breakpoints.getContig(), breakpoints.getInterval().getEnd(), breakpoints.getInterval().getEnd());
            final Collection<LargeSimpleSV> events = getEventsOnInterval(leftInterval, rightInterval, breakpoints.getInterval(), breakpoints, arguments.breakpointPadding);
            for (final LargeSimpleSV event : events) {
                if (event.getEventType() == SimpleSVType.TYPES.DUP_TANDEM) {
                    duplicationEventTree.put(event.getInterval(), event);
                }
            }
        }

        //Search links for tandem duplications
        final Iterator<SVIntervalTree.Entry<EvidenceTargetLink>> linkIter = intrachromosomalLinkTree.iterator();
        while (linkIter.hasNext()) {
            final EvidenceTargetLink link = linkIter.next().getValue();
            final SVInterval leftInterval = link.getPairedStrandedIntervals().getLeft().getInterval();
            final SVInterval rightInterval = link.getPairedStrandedIntervals().getRight().getInterval();
            if (leftInterval.getEnd() < rightInterval.getStart()) {
                final SVInterval callInterval = new SVInterval(leftInterval.getContig(), leftInterval.getEnd(), rightInterval.getStart());
                final Collection<LargeSimpleSV> events = getEventsOnInterval(leftInterval, rightInterval, callInterval, null, arguments.evidenceTargetLinkPadding);
                for (final LargeSimpleSV event : events) {
                    if (event.getEventType() == SimpleSVType.TYPES.DUP_TANDEM) {
                        duplicationEventTree.put(event.getInterval(), event);
                    }
                }
            }
        }


        /*
        //Remove conflicting calls
        final Set<LargeSimpleSV> callsToRemove = new HashSet<>();
        final Iterator<SVIntervalTree.Entry<LargeSimpleSV>> iterator = duplicationEventTree.iterator();
        while (iterator.hasNext()) {
            final SVIntervalTree.Entry<LargeSimpleSV> entry = iterator.next();
            final SVInterval interval = entry.getInterval();
            final LargeSimpleSV largeSimpleSV = entry.getValue();
            final Set<EvidenceTargetLink> supportingEvidence = new HashSet<>(largeSimpleSV.getSupportingEvidence());
            final List<LargeSimpleSV> conflictingCalls = Utils.stream(duplicationEventTree.overlappers(interval)).map(SVIntervalTree.Entry::getValue)
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
        */

        //Filter calls
        final SVIntervalTree<LargeSimpleSV> filteredCalls = new SVIntervalTree<>();
        Utils.stream(duplicationEventTree.iterator())
                //.filter(entry -> entry.getInterval().getContig() >= 0)  //TODO
                //.filter(entry -> !callsToRemove.contains(entry.getValue()))
                .filter(entry -> entry.getValue().getReadPairEvidence() > 0)
                .filter(entry -> entry.getValue().getReadPairEvidence() > Math.log10(entry.getInterval().getLength()) - 2)
                .filter(entry -> !blacklistTree.hasOverlapper(entry.getInterval()))
                .filter(entry -> {
                    final SVInterval interval = entry.getInterval();
                    final SVInterval start = new SVInterval(interval.getContig(), interval.getStart(), interval.getStart());
                    final SVInterval end = new SVInterval(interval.getContig(), interval.getEnd(), interval.getEnd());
                    return mappableIntervalTree.hasOverlapper(start) && mappableIntervalTree.hasOverlapper(end);
                }).filter(entry -> {
            final SVInterval interval = entry.getInterval();
            final SVInterval start = new SVInterval(interval.getContig(), interval.getStart(), interval.getStart());
            final SVInterval end = new SVInterval(interval.getContig(), interval.getEnd(), interval.getEnd());
            final SVInterval paddedStart = SVIntervalUtils.getPaddedInterval(start, 500, dictionary);
            final SVInterval paddedEnd = SVIntervalUtils.getPaddedInterval(end, 500, dictionary);
            return !highCoverageIntervalTree.hasOverlapper(paddedStart) && !highCoverageIntervalTree.hasOverlapper(paddedEnd);
        }).forEach(event -> filteredCalls.put(event.getInterval(), event.getValue()));

        final List<SVInterval> eventIntervals = Utils.stream(filteredCalls.iterator())
                .map(SVIntervalTree.Entry::getInterval)
                .collect(Collectors.toList());
        /*
        final List<GenomeLoc> eventIntervals = Utils.stream(filteredCalls.iterator()).flatMap(call -> {
            final int padding = 1000;
            final SVInterval leftInterval = SVIntervalUtils.getPaddedInterval(new SVInterval(call.getInterval().getContig(), call.getInterval().getStart(), call.getInterval().getStart()), padding, dictionary);
            final SVInterval rightInterval = SVIntervalUtils.getPaddedInterval(new SVInterval(call.getInterval().getContig(), call.getInterval().getEnd(), call.getInterval().getEnd()), padding, dictionary);
            return Stream.of(SVIntervalUtils.convertToGenomeLoc(leftInterval, dictionary), SVIntervalUtils.convertToGenomeLoc(rightInterval, dictionary));
        }).collect(Collectors.toList());
        */
        intervals.addAll(eventIntervals);

        final List<SimpleInterval> smallIntervals = stratifyIntervals(intervals, arguments.smallEventSize, arguments.mediumEventSize, dictionary);
        final List<SimpleInterval> mediumIntervals = stratifyIntervals(intervals, arguments.mediumEventSize, arguments.largeEventSize,dictionary);
        final List<SimpleInterval> largeIntervals = stratifyIntervals(intervals, arguments.largeEventSize, arguments.xlargeEventSize, dictionary);
        final List<SimpleInterval> extraLargeIntervals = stratifyIntervals(intervals, arguments.xlargeEventSize, Integer.MAX_VALUE, dictionary);

        final Map<String,List<SimpleInterval>> intervalsMap = new HashMap<>(SVUtils.hashMapCapacity(3));
        intervalsMap.put("S", smallIntervals);
        intervalsMap.put("M", mediumIntervals);
        intervalsMap.put("L", largeIntervals);
        intervalsMap.put("XL", extraLargeIntervals);

        for (final String key : intervalsMap.keySet()) {
            logger.info(key + " events: " + intervalsMap.get(key).size() + " intervals, " + intervalsMap.get(key).stream().mapToInt(SimpleInterval::getLengthOnReference).sum() + " bp");
        }

        return intervalsMap;
    }

    private static int getNumBins(final List<SVInterval> intervals, final int mediumIndex, final int largeIndex, final int xlargeIndex, final SAMSequenceDictionary dictionary) {
        final int numSmall = mergeIntervals(intervals.subList(0, mediumIndex).stream(), dictionary).mapToInt(interval -> interval.getLengthOnReference()).sum() / 10;
        final int numMedium =  mergeIntervals(intervals.subList(mediumIndex, largeIndex).stream(), dictionary).mapToInt(interval -> interval.getLengthOnReference()).sum() / 100;
        final int numLarge =  mergeIntervals(intervals.subList(largeIndex, xlargeIndex).stream(), dictionary).mapToInt(interval -> interval.getLengthOnReference()).sum() / 1000;
        final int numXLarge =  mergeIntervals(intervals.subList(xlargeIndex, intervals.size()).stream(), dictionary).mapToInt(interval -> interval.getLengthOnReference()).sum() / 10000;
        return numSmall + numMedium + numLarge + numXLarge;
    }

    private static List<SimpleInterval> stratifyIntervals(final Collection<SVInterval> intervals, final int minSize, final int maxSize, final SAMSequenceDictionary dictionary) {
        final Stream<SVInterval> filteredIntervals = getIntervalsInSizeRange(intervals, minSize, maxSize);
        return mergeIntervals(filteredIntervals, dictionary).collect(Collectors.toList());
    }

    private static Stream<SimpleInterval> mergeIntervals(final Stream<SVInterval> intervals, final SAMSequenceDictionary dictionary) {
        final List<GenomeLoc> genomeLocList = intervals.map(interval -> SVIntervalUtils.convertToGenomeLoc(interval, dictionary)).collect(Collectors.toList());
        return IntervalUtils.mergeIntervalLocations(genomeLocList, IntervalMergingRule.ALL).stream()
                .map(loc -> SVIntervalUtils.convertToSimpleInterval(loc))
                .sorted(IntervalUtils.getDictionaryOrderComparator(dictionary));
    }

    private static Stream<SVInterval> getIntervalsInSizeRange(final Collection<SVInterval> intervals, final int minSize, final int maxSize) {
        return intervals.stream().filter(interval -> interval.getLength() >= minSize && interval.getLength() < maxSize);
    }

    public static List<Tuple2<SVInterval,GATKRead>> getReads(final String inputPath, final JavaSparkContext ctx, final List<SimpleInterval> intervals, final SAMSequenceDictionary dictionary) {
        final TraversalParameters traversalParameters = new TraversalParameters(intervals, false);
        final ReadsSparkSource readsSource = new ReadsSparkSource(ctx, ValidationStringency.DEFAULT_STRINGENCY);
        final JavaRDD<GATKRead> reads = readsSource.getParallelReads(inputPath, null, traversalParameters, 8000000);
        return reads.filter(read -> !read.isUnmapped()).map(read -> new Tuple2<>(SVIntervalUtils.convertToSVInterval(new SimpleInterval(read.getContig(), read.getUnclippedStart(), read.getUnclippedEnd()), dictionary), read)).collect();
    }

}
