package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.VariantContext;
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
 * Base class for calling large (i.e. > 500 bp) non-complex structural variants a breakpoint pair on the same chromosome.
 * It requires the putative event interval to be known and makes calls by integrating evidence target links, copy ratios,
 * and copy ratio segments.
 */
public abstract class LargeSimpleSVFactory {

    protected final SVIntervalTree<EvidenceTargetLink> intrachromosomalLinkTree;
    protected final SVIntervalTree<EvidenceTargetLink> interchromosomalLinkTree;
    protected final SVIntervalTree<VariantContext> structuralVariantCallTree;
    protected final SVIntervalTree<GATKRead> contigTree;
    protected final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromReadDepthArgumentCollection arguments;
    protected final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector;
    protected final OverlapDetector<CopyRatio> copyRatioOverlapDetector;
    protected final SAMSequenceDictionary dictionary;

    public LargeSimpleSVFactory(final SVIntervalTree<EvidenceTargetLink> intrachromosomalLinkTree,
                                final SVIntervalTree<EvidenceTargetLink> interchromosomalLinkTree,
                                final SVIntervalTree<VariantContext> structuralVariantCallTree,
                                final SVIntervalTree<GATKRead> contigTree,
                                final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromReadDepthArgumentCollection arguments,
                                final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector,
                                final OverlapDetector<CopyRatio> copyRatioOverlapDetector,
                                final SAMSequenceDictionary dictionary) {
        Utils.nonNull(intrachromosomalLinkTree, "Intrachromosomal link tree cannot be null");
        Utils.nonNull(interchromosomalLinkTree, "Interchromosomal link tree cannot be null");
        Utils.nonNull(structuralVariantCallTree, "Structural variant call tree cannot be null");
        Utils.nonNull(contigTree, "Contig tree cannot be null");
        Utils.nonNull(arguments, "Arguments cannot be null");
        Utils.nonNull(copyRatioSegmentOverlapDetector, "Copy ratio segment overlap detector cannot be null");
        Utils.nonNull(copyRatioOverlapDetector, "Copy ratio overlap detector cannot be null");
        Utils.nonNull(dictionary, "Sequence dictionary cannot be null");
        this.intrachromosomalLinkTree = intrachromosomalLinkTree;
        this.interchromosomalLinkTree = interchromosomalLinkTree;
        this.structuralVariantCallTree = structuralVariantCallTree;
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
    @VisibleForTesting
    static List<CopyRatio> getCopyRatiosOnInterval(final SVInterval interval, final OverlapDetector<CopyRatio> overlapDetector,
                                                             final int binsToTrim, final SAMSequenceDictionary dictionary) {
        final SimpleInterval simpleInterval = SVIntervalUtils.convertToSimpleInterval(interval, dictionary);
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
    protected static int countUniqueCounterEvidence(final Collection<EvidenceTargetLink> counterEvidenceLinks, final Collection<EvidenceTargetLink> evidenceLinks,
                                                    final int minEvidenceCount, final Function<EvidenceTargetLink, Set<String>> evidenceTypeGetter) {
        final Set<String> counterEvidenceTemplates = counterEvidenceLinks.stream()
                .map(link -> evidenceTypeGetter.apply(link))
                .filter(set -> set.size() >= minEvidenceCount)
                .flatMap(Set::stream)
                .collect(Collectors.toSet());
        final Set<String> evidenceTemplates = evidenceLinks.stream()
                .map(link -> Stream.concat(link.getReadPairTemplateNames().stream(), link.getSplitReadTemplateNames().stream()))
                .flatMap(Function.identity()).collect(Collectors.toSet());
        return (int) counterEvidenceTemplates.stream().filter(template -> !evidenceTemplates.contains(template)).count();
    }

    private static int countUniqueCounterEvidenceReadPairs(final Collection<EvidenceTargetLink> counterEvidenceLinks, final Collection<EvidenceTargetLink> evidenceLinks, final int minEvidenceCount) {
        return countUniqueCounterEvidence(counterEvidenceLinks, evidenceLinks, minEvidenceCount, EvidenceTargetLink::getReadPairTemplateNames);
    }

    private static int countUniqueCounterEvidenceSplitReads(final Collection<EvidenceTargetLink> counterEvidenceLinks, final Collection<EvidenceTargetLink> evidenceLinks, final int minEvidenceCount) {
        return countUniqueCounterEvidence(counterEvidenceLinks, evidenceLinks, minEvidenceCount, EvidenceTargetLink::getSplitReadTemplateNames);
    }

    /**
     * Counts the number of unique split reads
     */
    protected static int countSplitReads(final Collection<EvidenceTargetLink> links) {
        return (int) links.stream().flatMap(link -> link.getSplitReadTemplateNames().stream()).distinct().count();
    }

    /**
     * Counts the number of unique read pairs
     */
    protected static int countReadPairs(final Collection<EvidenceTargetLink> links) {
        return (int) links.stream().flatMap(link -> link.getReadPairTemplateNames().stream()).distinct().count();
    }

    /**
     * Gets collection of evidence links whose target intervals align with the given left and right interval.
     *
     * @param leftInterval  Left interval that must overlap the link's left interval
     * @param rightInterval Right interval that much overlap the link's right interval
     * @param tree          Tree of evidence
     * @return Collection of overlapping links
     */
    protected static Collection<EvidenceTargetLink> getMatchingLinks(final SVInterval leftInterval, final SVInterval rightInterval, final SVIntervalTree<EvidenceTargetLink> tree) {
        final Collection<EvidenceTargetLink> leftOverlappingLinks = SVIntervalUtils.getOverlappingValuesOnInterval(leftInterval, tree);
        return leftOverlappingLinks.stream()
                .filter(link -> link.getPairedStrandedIntervals().getLeft().getInterval().overlaps(leftInterval)
                        && link.getPairedStrandedIntervals().getRight().getInterval().overlaps(rightInterval))
                .collect(Collectors.toList());
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
    protected static Collection<EvidenceTargetLink> localOverlappingLinks(final SVInterval interval, final SVIntervalTree<EvidenceTargetLink> tree, final int localRange, final SAMSequenceDictionary dictionary) {
        final SVInterval localInterval = SVIntervalUtils.getPaddedInterval(interval, localRange, dictionary);
        final Collection<EvidenceTargetLink> overlappingLinks = SVIntervalUtils.getOverlappingValuesOnInterval(localInterval, tree);
        return overlappingLinks.stream().filter(link -> {
            final SVInterval linkInterval = SVIntervalUtils.getOuterIntrachromosomalLinkInterval(link);
            final SVInterval leftInterval = link.getPairedStrandedIntervals().getLeft().getInterval();
            final SVInterval rightInterval = link.getPairedStrandedIntervals().getRight().getInterval();
            return linkInterval.overlaps(interval) && (localInterval.overlaps(leftInterval) || localInterval.overlaps(rightInterval));
        }).collect(Collectors.toList());
    }

    /**
     * Returns true if the link is of proper orientation to support the event
     */
    protected abstract boolean hasSupportingEvidenceOrientation(final EvidenceTargetLink link);

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
     * Returns number of supporting evidence read pairs
     */
    protected abstract int countSupportingEvidenceReadPairs(final Collection<EvidenceTargetLink> links);

    /**
     * Returns number of supporting evidence split reads
     */
    protected abstract int countSupportingEvidenceSplitReads(final Collection<EvidenceTargetLink> links);

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
                                              final IntrachromosomalBreakpointPair breakpoints,
                                              final Collection<EvidenceTargetLink> supportingEvidence);

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
    public LargeSimpleSV call(final SVInterval leftInterval,
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
        final Collection<EvidenceTargetLink> overlappingLinks = getMatchingLinks(paddedLeftInterval, paddedRightInterval, intrachromosomalLinkTree);
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
        final Collection<EvidenceTargetLink> counterEvidenceLinks = localOverlappingLinks(outerInterval, intrachromosomalLinkTree, arguments.localCounterevidenceRange, dictionary);
        counterEvidenceLinks.addAll(SVIntervalUtils.getOverlappingValuesOnInterval(outerInterval, interchromosomalLinkTree));
        counterEvidenceLinks.removeAll(evidenceLinks);

        //Tally evidence and counterevidence
        final int readPairEvidence = countSupportingEvidenceReadPairs(evidenceLinks);
        if (readPairEvidence == 0) return null; //TODO
        final int splitReadEvidence = countSupportingEvidenceSplitReads(evidenceLinks);
        final int readPairCounterEvidence = countUniqueCounterEvidenceReadPairs(counterEvidenceLinks, evidenceLinks, arguments.minCountervidenceClusterSize);
        final int splitReadCounterEvidence = countUniqueCounterEvidenceSplitReads(counterEvidenceLinks, evidenceLinks, arguments.minCountervidenceClusterSize);

        //Score the event and reject if too small
        final double eventScore = LargeSimpleSV.computeScore(readPairEvidence, splitReadEvidence, readPairCounterEvidence, splitReadCounterEvidence, arguments.counterEvidencePseudocount, breakpoints != null);
        if (eventScore < arguments.minScore) return null;

        final String contigName = sequence.getSequenceName();
        return getNewSV(callInterval.getStart(), callInterval.getEnd(), contigId, contigName, readPairEvidence, splitReadEvidence, readPairCounterEvidence, splitReadCounterEvidence, breakpoints, evidenceLinks);

        /*
        //Test if the copy ratio bins indicate we should reject this call
        final List<CopyRatio> eventBins = getCopyRatiosOnInterval(innerInterval, copyRatioOverlapDetector, arguments.copyRatioBinTrimming, dictionary);
        if (isInvalidCoverage(eventBins)) return null;

        //Test if the event matches a model segments call
        final String contigName = sequence.getSequenceName();
        final Set<CalledCopyRatioSegment> overlappingSegments = copyRatioSegmentOverlapDetector.getOverlaps(SVIntervalUtils.convertToSimpleInterval(innerInterval, dictionary));
        final boolean eventSegmentSupport = supportedBySegmentCalls(innerInterval, overlappingSegments, dictionary);
        if (eventSegmentSupport) {
            return getNewSV(callInterval.getStart(), callInterval.getEnd(), contigId, contigName, readPairEvidence,
                    splitReadEvidence, readPairCounterEvidence, splitReadCounterEvidence, breakpoints);
        }

        //Run copy number state HMM over copy ratios and test if the states (a la Viterbi) match the valid states for the event type
        if (supportedByHMM(eventBins)) {
            return getNewSV(callInterval.getStart(), callInterval.getEnd(), contigId, contigName, readPairEvidence, splitReadEvidence, readPairCounterEvidence, splitReadCounterEvidence, breakpoints);
        }

        return null;
        */
    }

    private boolean checkFlankingSegments(final SVInterval interval, final Set<CalledCopyRatioSegment> overlappingSegments, final SAMSequenceDictionary dictionary) {
        final int neutralBases = overlappingSegments.stream().filter(segment -> segment.getCall() == CalledCopyRatioSegment.Call.NEUTRAL)
                .mapToInt(segment -> SVIntervalUtils.convertToSVInterval(segment.getInterval(), dictionary).overlapLen(interval)).sum();
        return neutralBases / (double) interval.getLength() >= arguments.minSegmentOverlap;
    }

    /**
     * Returns collection of the links with the proper orientation that supports the event
     */
    protected Collection<EvidenceTargetLink> getLinksWithEvidenceOrientation(final Collection<EvidenceTargetLink> links) {
        return links.stream().filter(link -> hasSupportingEvidenceOrientation(link)).collect(Collectors.toList());
    }

    /**
     * Tests whether the copy ratios support an event with an HMM
     *
     * @param copyRatioBins Copy ratios over the interval
     * @return True if supported
     */
    protected boolean supportedByHMM(final List<CopyRatio> copyRatioBins) {
        if (copyRatioBins.isEmpty()) return false;
        final List<Double> copyRatios = copyRatioBins.stream().map(CopyRatio::getLog2CopyRatioValue).map(val -> Math.pow(2.0, val)).collect(Collectors.toList());
        final int numStates = Math.min(arguments.hmmMaxStates, copyRatios.stream().mapToInt(val -> (int) (2 * val)).max().getAsInt() + 1);
        final RealVector copyNumberPriors = CopyNumberHMM.uniformPrior(numStates);
        final List<Integer> positionsList = CopyNumberHMM.positionsList(copyRatios.size());
        final CopyNumberHMM copyNumberHMM = new CopyNumberHMM(copyNumberPriors, arguments.hmmTransitionProb);
        final List<Integer> copyNumberStates = ViterbiAlgorithm.apply(copyRatios, positionsList, copyNumberHMM);
        return testHMMState(copyNumberStates, numStates);
    }

    /**
     * Tests if the HMM state path contains a sufficient proportion of valid states
     *
     * @param states    State path
     * @param numStates Number of HMM states
     * @return True if the threshold is met
     */
    private boolean testHMMState(final List<Integer> states, final int numStates) {
        return validStateFrequency(states, numStates) >= arguments.hmmValidStatesMinFraction * states.size();
    }

    /**
     * Counts the number of states in the path that are one of the valid states
     */
    private int validStateFrequency(final List<Integer> states, final int numStates) {
        final Set<Integer> validStates = getValidHMMCopyStates(numStates);
        return (int) states.stream().filter(validStates::contains).count();
    }
}
