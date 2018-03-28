package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.commons.math3.linear.RealVector;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.LargeSimpleSV;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.utils.IntrachromosomalBreakpointPair;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.ViterbiAlgorithm;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Base class for calling large (i.e. > 500 bp) non-complex structural variants with breakpoints on the same chromosome.
 * It requires the putative event interval to be known and makes calls by integrating evidence target links, copy ratios,
 * and copy ratio segments.
 */
public abstract class LargeSimpleSVFactory {

    protected final SVIntervalTree<EvidenceTargetLink> intrachromosomalLinkTree;
    protected final SVIntervalTree<EvidenceTargetLink> interchromosomalLinkTree;
    protected final SVIntervalTree<GATKRead> contigTree;
    protected final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromReadDepth arguments;
    protected final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector;
    protected final OverlapDetector<CopyRatio> copyRatioOverlapDetector;
    protected final SAMSequenceDictionary dictionary;

    public LargeSimpleSVFactory(final SVIntervalTree<EvidenceTargetLink> intrachromosomalLinkTree,
                                final SVIntervalTree<EvidenceTargetLink> interchromosomalLinkTree,
                                final SVIntervalTree<GATKRead> contigTree,
                                final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromReadDepth arguments,
                                final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector,
                                final OverlapDetector<CopyRatio> copyRatioOverlapDetector,
                                final SAMSequenceDictionary dictionary) {
        Utils.nonNull(intrachromosomalLinkTree, "Intrachromosomal link tree cannot be null");
        Utils.nonNull(interchromosomalLinkTree, "Interchromosomal link tree cannot be null");
        Utils.nonNull(contigTree, "Contig tree cannot be null");
        Utils.nonNull(arguments, "Arguments cannot be null");
        Utils.nonNull(copyRatioSegmentOverlapDetector, "Copy ratio segment overlap detector cannot be null");
        Utils.nonNull(copyRatioOverlapDetector, "Copy ratio overlap detector cannot be null");
        Utils.nonNull(dictionary, "Sequence dictionary cannot be null");
        this.intrachromosomalLinkTree = intrachromosomalLinkTree;
        this.interchromosomalLinkTree = interchromosomalLinkTree;
        this.contigTree = contigTree;
        this.arguments = arguments;
        this.copyRatioSegmentOverlapDetector = copyRatioSegmentOverlapDetector;
        this.copyRatioOverlapDetector = copyRatioOverlapDetector;
        this.dictionary = dictionary;
    }

    /**
     * Returns list of ordered CopyRatio objects on the given interval
     *
     * @param interval        Interval over which to retrieve bins
     * @param overlapDetector Copy ratio overlap detector
     * @param binsToTrim      Number of bins to trim from either side
     * @param dictionary      Sequence dictionary
     * @return List of copy ratios
     */
    private static List<CopyRatio> getCopyRatiosOnInterval(final SVInterval interval, final OverlapDetector<CopyRatio> overlapDetector,
                                                           final int binsToTrim, final SAMSequenceDictionary dictionary) {
        final int start = interval.getStart();
        final int end = interval.getEnd();
        final SAMSequenceRecord sequence = dictionary.getSequence(interval.getContig());
        if (sequence == null) {
            throw new IllegalArgumentException("Could not find contig with index " + interval.getContig() + " in the sequence dictionary");
        }
        final SimpleInterval simpleInterval = new SimpleInterval(sequence.getSequenceName(), start + 1, end + 1);
        if (simpleInterval.size() == 0) {
            return Collections.emptyList();
        }
        final List<CopyRatio> copyRatios = overlapDetector.getOverlaps(simpleInterval).stream().collect(Collectors.toList());

        if (copyRatios.size() <= 2 * binsToTrim) return Collections.emptyList();
        Collections.sort(copyRatios, Comparator.comparing(CopyRatio::getStart));
        return copyRatios.subList(binsToTrim, copyRatios.size() - binsToTrim);
    }

    /**
     * Counts the number of valid unique counter-evidence read templates that are not supporting evidence
     *
     * @param counterEvidenceLinks Collection of counter-evidence
     * @param evidenceLinks        Collection of supporting evidence
     * @param minEvidenceCount     Counter-evidence links with fewer than this many reads will be filtered out
     * @param evidenceTypeGetter   Function that returns the type of evidence (i.e. split reads or read pairs)
     * @return Counter-evidence count
     */
    private static int countUniqueCounterEvidence(final Collection<EvidenceTargetLink> counterEvidenceLinks, final Collection<EvidenceTargetLink> evidenceLinks,
                                                  final int minEvidenceCount, final Function<EvidenceTargetLink, Set<String>> evidenceTypeGetter) {
        final Collection<String> counterEvidenceTemplates = counterEvidenceLinks.stream()
                .map(link -> evidenceTypeGetter.apply(link))
                .filter(set -> set.size() >= minEvidenceCount)
                .flatMap(Set::stream)
                .collect(Collectors.toList());
        final Set<String> evidenceTemplates = evidenceLinks.stream()
                .map(link -> Stream.concat(link.getReadPairTemplateNames().stream(), link.getSplitReadTemplateNames().stream()))
                .flatMap(Function.identity()).collect(Collectors.toSet());
        return (int) counterEvidenceTemplates.stream().filter(evidenceTemplates::contains).count();
    }

    private static int countUniqueCounterEvidenceReadPairs(final Collection<EvidenceTargetLink> counterEvidenceLinks, final Collection<EvidenceTargetLink> evidenceLinks, final int minEvidenceCount) {
        return countUniqueCounterEvidence(counterEvidenceLinks, evidenceLinks, minEvidenceCount, EvidenceTargetLink::getReadPairTemplateNames);
    }

    private static int countUniqueCounterEvidenceSplitReads(final Collection<EvidenceTargetLink> counterEvidenceLinks, final Collection<EvidenceTargetLink> evidenceLinks, final int minEvidenceCount) {
        return countUniqueCounterEvidence(counterEvidenceLinks, evidenceLinks, minEvidenceCount, EvidenceTargetLink::getSplitReadTemplateNames);
    }

    /**
     * Returns true if the link is of proper orientation to support the event
     */
    protected abstract boolean isEvidenceOrientation(final EvidenceTargetLink link);

    /**
     * Returns true if this event should be filtered based on the copy ratio bins found in the interval
     */
    protected abstract boolean isInvalidCoverage(final List<CopyRatio> copyRatios);

    /**
     * Gets valid copy ratio states that would support the presence of an event
     */
    protected abstract Set<Integer> getValidHMMCopyStates(final int numStates);

    /**
     * Determines if the event is supported by copy ratio segment calls
     *
     * @param interval            Event interval
     * @param overlappingSegments Segments overlapping the interval
     * @param dictionary          Sequence dictionary
     * @return True if the event is supported by the overlapping segments
     */
    protected abstract boolean supportedBySegmentCalls(final SVInterval interval, final Set<CalledCopyRatioSegment> overlappingSegments, final SAMSequenceDictionary dictionary);

    /**
     * Gets a new call object corresponding to the factory's type
     */
    protected abstract LargeSimpleSV getNewSV(final int start,
                                              final int end,
                                              final int contigId,
                                              final String contig,
                                              final int readPairEvidence,
                                              final int splitReadEvidence,
                                              final int readPairCounterEvidence,
                                              final int splitReadCounterEvidence,
                                              final List<CopyRatio> coverage,
                                              final List<Integer> copyNumberStates,
                                              final IntrachromosomalBreakpointPair breakpoints);

    /**
     * Generates a call for an event with breakpoint supporting evidence on the given intervals
     *
     * @param leftInterval    Interval overlapping left breakpoint
     * @param rightInterval   Interval overlapping right breakpoint
     * @param callInterval    Interval to call if successful
     * @param breakpoints     Associated breakpoint pair, if any
     * @param evidencePadding Amount to pad the left and right intervals when searching for supporting evidence
     * @return An event call or null if unsuccessful
     */
    public LargeSimpleSV create(final SVInterval leftInterval,
                                final SVInterval rightInterval,
                                final SVInterval callInterval,
                                final IntrachromosomalBreakpointPair breakpoints,
                                final int evidencePadding) {

        Utils.nonNull(leftInterval, "Left interval cannot be null");
        Utils.nonNull(rightInterval, "Right interval cannot be null");
        if (leftInterval.getContig() != rightInterval.getContig()) {
            throw new IllegalArgumentException("Left and right intervals must be on the same contig");
        }
        if (breakpoints != null && leftInterval.getContig() != breakpoints.getContig()) {
            throw new IllegalArgumentException("Intervals must be on the same contig as the breakpoints");
        }
        final int contigId = leftInterval.getContig();
        final SAMSequenceRecord sequence = dictionary.getSequence(contigId);
        if (sequence == null) {
            throw new IllegalArgumentException("Could not find interval contig with index " + contigId + " in the sequence dictionary");
        }

        //Get evidence links whose left and right intervals overlap with the padded intervals and have proper strandedness for the event type
        final SVInterval paddedLeftInterval = SVIntervalUtils.getPaddedInterval(leftInterval, evidencePadding, dictionary);
        final SVInterval paddedRightInterval = SVIntervalUtils.getPaddedInterval(rightInterval, evidencePadding, dictionary);
        final Collection<EvidenceTargetLink> overlappingLinks = getOverlappingLinks(paddedLeftInterval, paddedRightInterval, intrachromosomalLinkTree, dictionary);
        final Collection<EvidenceTargetLink> evidenceLinks = getLinksWithEvidenceOrientation(overlappingLinks);
        if (evidenceLinks.isEmpty()) return null;

        //Get "outer" and "inner" intervals
        final SVInterval outerInterval = new SVInterval(contigId, leftInterval.getStart(), rightInterval.getEnd());
        final SVInterval innerInterval;
        if (rightInterval.getStart() > leftInterval.getEnd()) {
            innerInterval = new SVInterval(contigId, leftInterval.getEnd(), rightInterval.getStart());
        } else {
            innerInterval = new SVInterval(contigId, leftInterval.getEnd(), leftInterval.getEnd() + 1);
        }

        //Get overlapping counterevidence links that suggest a more complex signature
        final Collection<EvidenceTargetLink> counterEvidenceLinks = localOverlappingLinks(outerInterval, intrachromosomalLinkTree, arguments.SPANNING_COUNTEREVIDENCE_RANGE, dictionary);
        counterEvidenceLinks.addAll(SVIntervalUtils.getOverlappingLinksOnInterval(outerInterval, interchromosomalLinkTree));
        counterEvidenceLinks.removeAll(evidenceLinks);

        //Tally evidence and counterevidence
        final int readPairEvidence = readPairEvidence(evidenceLinks);
        final int splitReadEvidence = splitReadEvidence(evidenceLinks);
        final int readPairCounterEvidence = countUniqueCounterEvidenceReadPairs(counterEvidenceLinks, evidenceLinks, arguments.MIN_LINK_COUNTEREVIDENCE);
        final int splitReadCounterEvidence = countUniqueCounterEvidenceSplitReads(counterEvidenceLinks, evidenceLinks, arguments.MIN_LINK_COUNTEREVIDENCE);

        //Score the event and reject if too small
        final double evidenceToCounterEvidenceRatio = LargeSimpleSV.computeScore(readPairEvidence, splitReadEvidence, readPairCounterEvidence, splitReadCounterEvidence, arguments.COUNTEREVIDENCE_PSEUDOCOUNT);
        if (evidenceToCounterEvidenceRatio < arguments.MIN_EVIDENCE_TO_COUNTEREVIDENCE_RATIO) return null;

        //Test if the event matches a model segments call
        final String contigName = sequence.getSequenceName();
        final Set<CalledCopyRatioSegment> overlappingSegments = copyRatioSegmentOverlapDetector.getOverlaps(new SimpleInterval(contigName, outerInterval.getStart(), outerInterval.getEnd()));
        boolean supportedBySegmentCalls = supportedBySegmentCalls(innerInterval, overlappingSegments, dictionary);
        final SVInterval hmmInterval = SVIntervalUtils.getPaddedInterval(innerInterval, arguments.HMM_PADDING, dictionary);
        final List<CopyRatio> copyRatioBins = getCopyRatiosOnInterval(hmmInterval, copyRatioOverlapDetector, arguments.COPY_NUMBER_BIN_TRIMMING, dictionary);
        if (isInvalidCoverage(copyRatioBins)) return null;
        if (supportedBySegmentCalls) {
            return getNewSV(callInterval.getStart(), callInterval.getEnd(), contigId, contigName, readPairEvidence, splitReadEvidence, readPairCounterEvidence, splitReadCounterEvidence,
                    Collections.emptyList(), Collections.emptyList(), breakpoints);
        }

        //Run copy number state HMM over copy ratios and test if the states (a la Viterbi) match the valid states for the event type
        if (copyRatioBins.isEmpty()) return null;
        final List<Double> copyRatios = copyRatioBins.stream().map(CopyRatio::getLog2CopyRatioValue).collect(Collectors.toList());
        final int numStates = Math.min(2 * arguments.MAX_COPY_RATIO_STATE + 1, copyRatios.stream().mapToInt(val -> (int) Math.pow(2.0, val)).max().getAsInt() + 2);
        final RealVector copyNumberPriors = CopyNumberHMM.getCopyNumberPrior(numStates);
        final List<Integer> positionsList = CopyNumberHMM.getCopyNumberHMMPositions(copyRatios.size());
        final CopyNumberHMM copyNumberHMM = new CopyNumberHMM(copyNumberPriors, arguments.COPY_NUMBER_HMM_ALPHA);
        final List<Integer> copyNumberStates = ViterbiAlgorithm.apply(copyRatios, positionsList, copyNumberHMM);
        if (testHMMState(copyNumberStates, numStates, arguments.MIN_EVENT_HMM_COVERAGE)) {
            return getNewSV(callInterval.getStart(), callInterval.getEnd(), contigId, contigName, readPairEvidence, splitReadEvidence, readPairCounterEvidence, splitReadCounterEvidence, copyRatioBins, copyNumberStates, breakpoints);
        }

        return null;
    }

    /**
     * Counts the number of unique split reads
     */
    protected int splitReadEvidence(final Collection<EvidenceTargetLink> links) {
        return (int) links.stream().flatMap(link -> link.getSplitReadTemplateNames().stream()).distinct().count();
    }

    /**
     * Counts the number of unique read pairs
     */
    protected int readPairEvidence(final Collection<EvidenceTargetLink> links) {
        return (int) links.stream().filter(link -> isEvidenceOrientation(link)).flatMap(link -> link.getReadPairTemplateNames().stream()).distinct().count();
    }

    /**
     * Returns collection of the links with the proper orientation that supports the event
     */
    protected Collection<EvidenceTargetLink> getLinksWithEvidenceOrientation(final Collection<EvidenceTargetLink> links) {
        return links.stream().filter(link -> isEvidenceOrientation(link)).collect(Collectors.toList());
    }

    /**
     * Gets collection of evidence links whose target intervals align with the given left and right interval.
     *
     * @param leftInterval  Left interval that must overlap the link's left interval
     * @param rightInterval Right interval that much overlap the link's right interval
     * @param tree          Tree of evidence
     * @param dictionary    Sequence dictionary
     * @return Collection of overlapping links
     */
    private Collection<EvidenceTargetLink> getOverlappingLinks(final SVInterval leftInterval, final SVInterval rightInterval, final SVIntervalTree<EvidenceTargetLink> tree, final SAMSequenceDictionary dictionary) {
        final Collection<EvidenceTargetLink> leftOverlappingLinks = SVIntervalUtils.getOverlappingLinksOnInterval(leftInterval, tree);
        return leftOverlappingLinks.stream().filter(link -> link.getPairedStrandedIntervals().getRight().getInterval().overlaps(rightInterval)).collect(Collectors.toList());
    }

    /**
     * Gets collection of evidence links that overlap the given interval "locally." This means that at least one of the
     * evidence left or right intervals must overlap the interval after being padded by some amount (i.e. the evidence cannot
     * completely span the padded interval). This approach is useful for finding counter-evidence on the interval without
     * pulling in evidence from irrelevant distant events.
     *
     * @param interval   Interval to retrieve evidence on
     * @param tree       Tree of evidence
     * @param localRange Number of bases to pad the given interval
     * @param dictionary Sequence dictionary
     * @return Collection of locally overlapping links
     */
    private Collection<EvidenceTargetLink> localOverlappingLinks(final SVInterval interval, final SVIntervalTree<EvidenceTargetLink> tree, final int localRange, final SAMSequenceDictionary dictionary) {
        final Collection<EvidenceTargetLink> overlappingLinks = SVIntervalUtils.getOverlappingLinksOnInterval(interval, tree);
        final SVInterval localInterval = SVIntervalUtils.getPaddedInterval(interval, localRange, dictionary);
        return overlappingLinks.stream().filter(link -> {
            final SVInterval linkInterval = SVIntervalUtils.getOuterIntrachromosomalLinkInterval(link);
            final SVInterval leftInterval = link.getPairedStrandedIntervals().getLeft().getInterval();
            final SVInterval rightInterval = link.getPairedStrandedIntervals().getRight().getInterval();
            return SVIntervalUtils.containsInterval(localInterval, linkInterval) || interval.overlaps(leftInterval) || interval.overlaps(rightInterval);
        }).collect(Collectors.toList());
    }

    /**
     * Tests if the HMM state path contains a sufficient proportion of valid states
     *
     * @param states              State path
     * @param numStates           Number of HMM states
     * @param minEventHMMCoverage Minimum proportion of valid states
     * @return True if the threshold is met
     */
    private boolean testHMMState(final List<Integer> states, final int numStates, final double minEventHMMCoverage) {
        return validStateFrequency(states, getValidHMMCopyStates(numStates)) >= minEventHMMCoverage * states.size();
    }

    /**
     * Counts the number of states in the path that are one of the valid states
     */
    private int validStateFrequency(final List<Integer> states, final Set<Integer> validStates) {
        return (int) states.stream().filter(validStates::contains).count();
    }
}
