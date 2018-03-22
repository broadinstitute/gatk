package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.math3.linear.RealVector;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.IntrachromosomalBreakpointPair;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype.CopyNumberHMM;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.hmm.ViterbiAlgorithm;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.stream.Collectors;

public abstract class SimpleSVFactory {

    protected final SVIntervalTree<EvidenceTargetLink> intrachromosomalLinkTree;
    protected final SVIntervalTree<EvidenceTargetLink> interchromosomalLinkTree;
    protected final SVIntervalTree<GATKRead> contigTree;
    protected final StructuralVariationDiscoveryArgumentCollection.SimpleVariantDiscoveryArgumentCollection arguments;
    protected final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector;
    protected final OverlapDetector<CopyRatio> readDepthOverlapDetector;
    protected final SAMSequenceDictionary dictionary;

    public SimpleSVFactory(final SVIntervalTree<EvidenceTargetLink> intrachromosomalLinkTree,
                           final SVIntervalTree<EvidenceTargetLink> interchromosomalLinkTree,
                           final SVIntervalTree<GATKRead> contigTree,
                           final StructuralVariationDiscoveryArgumentCollection.SimpleVariantDiscoveryArgumentCollection arguments,
                           final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector,
                           final OverlapDetector<CopyRatio> readDepthOverlapDetector,
                           final SAMSequenceDictionary dictionary) {
        this.intrachromosomalLinkTree = intrachromosomalLinkTree;
        this.interchromosomalLinkTree = interchromosomalLinkTree;
        this.contigTree = contigTree;
        this.arguments = arguments;
        this.copyRatioSegmentOverlapDetector = copyRatioSegmentOverlapDetector;
        this.readDepthOverlapDetector = readDepthOverlapDetector;
        this.dictionary = dictionary;
    }

    protected static int readPairCounterEvidence(final List<EvidenceTargetLink> counterEvidenceLinks, final List<EvidenceTargetLink> evidenceLinks, final int minEvidenceCount) {
        final List<String> counterEvidenceTemplates = counterEvidenceLinks.stream().filter(link -> link.getReadPairs() >= minEvidenceCount).flatMap(link -> link.getReadPairTemplateNames().stream()).collect(Collectors.toList());
        final Set<String> evidenceTemplates = evidenceLinks.stream().map(link -> {
            final Collection<String> linkTemplates = new ArrayList<>(link.getReadPairTemplateNames());
            linkTemplates.addAll(link.getSplitReadTemplateNames());
            return linkTemplates;
        }).flatMap(Collection::stream).collect(Collectors.toSet());
        counterEvidenceTemplates.removeAll(evidenceTemplates);
        return counterEvidenceTemplates.size();
    }

    protected abstract boolean isEvidenceOrientation(final EvidenceTargetLink link);

    protected abstract boolean isInvalidCoverage(final List<CopyRatio> copyRatios);

    protected abstract Set<Integer> getValidHMMCopyStates(final int numStates);

    protected abstract boolean supportedBySegmentCalls(final SVInterval interval, final Set<CalledCopyRatioSegment> overlappingSegments);

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

    public LargeSimpleSV create(final SVInterval leftInterval,
                                final SVInterval rightInterval,
                                final SVInterval callInterval,
                                final IntrachromosomalBreakpointPair breakpoints) {

        //Get evidence links whose left and right intervals overlap with the input intervals and have proper strandedness for the event type
        final List<EvidenceTargetLink> overlappingLinks = overlappingLinks(leftInterval, rightInterval, intrachromosomalLinkTree, dictionary);
        final List<EvidenceTargetLink> evidenceLinks = getLinksWithEvidenceOrientation(overlappingLinks);
        if (evidenceLinks.isEmpty()) return null;

        final int contigId = leftInterval.getContig();
        final String contig = dictionary.getSequence(contigId).getSequenceName();

        //Get "outer" and "inner" intervals
        final SVInterval outerInterval = new SVInterval(contigId, leftInterval.getStart(), rightInterval.getEnd());
        final SVInterval innerInterval;
        if (rightInterval.getStart() > leftInterval.getEnd()) {
            innerInterval = new SVInterval(contigId, leftInterval.getEnd(), rightInterval.getStart());
        } else {
            innerInterval = new SVInterval(contigId, leftInterval.getEnd(), leftInterval.getEnd() + 1);
        }

        //Get overlapping counterevidence links that suggest a more complex signature
        final List<EvidenceTargetLink> localOverlappingLinks = localOverlappingLinks(outerInterval, intrachromosomalLinkTree, arguments.SPANNING_COUNTEREVIDENCE_RANGE, dictionary);
        localOverlappingLinks.addAll(IntervalUtils.getOverlappingLinksOnInterval(outerInterval, interchromosomalLinkTree));
        final List<EvidenceTargetLink> counterEvidenceLinks = new ArrayList<>(CollectionUtils.subtract(localOverlappingLinks, evidenceLinks));

        //Tally evidence and counterevidence
        final int readPairEvidence = readPairEvidence(evidenceLinks);
        final int splitReadEvidence = splitReadEvidence(evidenceLinks);
        final int readPairCounterEvidence = readPairCounterEvidence(counterEvidenceLinks, evidenceLinks, arguments.MIN_LINK_COUNTEREVIDENCE);
        final int splitReadCounterEvidence = splitReadCounterEvidence(counterEvidenceLinks, evidenceLinks, arguments.MIN_LINK_COUNTEREVIDENCE);

        //Score the event and reject if too small
        final double evidenceToCounterEvidenceRatio = LargeSimpleSV.computeScore(readPairEvidence, splitReadEvidence, readPairCounterEvidence, splitReadCounterEvidence, arguments.COUNTEREVIDENCE_PSEUDOCOUNT);
        if (evidenceToCounterEvidenceRatio < arguments.MIN_EVIDENCE_TO_COUNTEREVIDENCE_RATIO) return null;

        //Test if the event matches a model segments call
        final Set<CalledCopyRatioSegment> overlappingSegments = copyRatioSegmentOverlapDetector.getOverlaps(new SimpleInterval(contig, outerInterval.getStart(), outerInterval.getEnd()));
        boolean supportedBySegmentCalls = supportedBySegmentCalls(outerInterval, overlappingSegments);
        final SVInterval hmmInterval = IntervalUtils.getPaddedInterval(innerInterval, arguments.HMM_PADDING, dictionary);
        final List<CopyRatio> copyRatioBins = IntervalUtils.getCopyRatiosOnInterval(hmmInterval, readDepthOverlapDetector, arguments.COPY_NUMBER_BIN_TRIMMING, contigId, dictionary);
        if (isInvalidCoverage(copyRatioBins)) return null;
        if (supportedBySegmentCalls) {
            return getNewSV(callInterval.getStart(), callInterval.getEnd(), contigId, contig, readPairEvidence, splitReadEvidence, readPairCounterEvidence, splitReadCounterEvidence,
                    Collections.emptyList(), Collections.emptyList(), breakpoints);
        }

        //Run copy number state HMM over copy ratios and test if the states (a la Viterbi) match the valid states for the event type
        if (copyRatioBins.isEmpty()) return null;
        final List<Double> copyRatios = copyRatioBins.stream().map(CopyRatio::getLog2CopyRatioValue).collect(Collectors.toList());
        final int numStates = Math.min(2 * arguments.highDepthCoveragePeakFactor + 1, copyRatios.stream().mapToInt(val -> (int) Math.pow(2.0, val)).max().getAsInt() + 2);
        final RealVector copyNumberPriors = CopyNumberHMM.getCopyNumberPrior(numStates);
        final List<Integer> positionsList = CopyNumberHMM.getCopyNumberHMMPositions(copyRatios.size());
        final CopyNumberHMM copyNumberHMM = new CopyNumberHMM(copyNumberPriors, arguments.COPY_NUMBER_HMM_ALPHA);
        final List<Integer> copyNumberStates = ViterbiAlgorithm.apply(copyRatios, positionsList, copyNumberHMM);
        if (testHMMState(copyNumberStates, numStates, arguments.MIN_EVENT_HMM_COVERAGE)) {
            return getNewSV(callInterval.getStart(), callInterval.getEnd(), contigId, contig, readPairEvidence, splitReadEvidence, readPairCounterEvidence, splitReadCounterEvidence, copyRatioBins, copyNumberStates, breakpoints);
        }

        return null;
    }

    protected int splitReadEvidence(final List<EvidenceTargetLink> links) {
        return (int) links.stream().flatMap(link -> link.getSplitReadTemplateNames().stream()).distinct().count();
    }

    protected int readPairEvidence(final List<EvidenceTargetLink> links) {
        return (int) links.stream().filter(link -> isEvidenceOrientation(link)).flatMap(link -> link.getReadPairTemplateNames().stream()).distinct().count();
    }

    protected List<EvidenceTargetLink> getLinksWithEvidenceOrientation(final List<EvidenceTargetLink> links) {
        return links.stream().filter(link -> isEvidenceOrientation(link)).collect(Collectors.toList());
    }

    private List<EvidenceTargetLink> overlappingLinks(final SVInterval leftBreakpointInterval, final SVInterval rightBreakpointInterval, final SVIntervalTree<EvidenceTargetLink> tree, final SAMSequenceDictionary dictionary) {
        final List<EvidenceTargetLink> leftOverlappingLinks = IntervalUtils.getOverlappingLinksOnInterval(leftBreakpointInterval, tree);
        return leftOverlappingLinks.stream().filter(link -> IntervalUtils.linkEndsOverlapIntervals(leftBreakpointInterval, rightBreakpointInterval, link)).collect(Collectors.toList());
    }

    private List<EvidenceTargetLink> localOverlappingLinks(final SVInterval breakpointInterval, final SVIntervalTree<EvidenceTargetLink> tree, final int localRange, final SAMSequenceDictionary dictionary) {
        final List<EvidenceTargetLink> overlappingLinks = IntervalUtils.getOverlappingLinksOnInterval(breakpointInterval, tree);
        final SVInterval localInterval = IntervalUtils.getPaddedInterval(breakpointInterval, localRange, dictionary);
        return overlappingLinks.stream().filter(link -> {
            final SVInterval linkInterval = IntervalUtils.getOuterLinkInterval(link);
            final SVInterval leftInterval = link.getPairedStrandedIntervals().getLeft().getInterval();
            final SVInterval rightInterval = link.getPairedStrandedIntervals().getRight().getInterval();
            return IntervalUtils.containsInterval(localInterval, linkInterval) || breakpointInterval.overlaps(leftInterval) || breakpointInterval.overlaps(rightInterval);
        }).collect(Collectors.toList());
    }

    private int splitReadCounterEvidence(final List<EvidenceTargetLink> counterEvidenceLinks, final List<EvidenceTargetLink> evidenceLinks, final int minEvidenceCount) {
        final List<String> counterEvidenceTemplates = counterEvidenceLinks.stream().filter(link -> link.getSplitReads() >= minEvidenceCount).flatMap(link -> link.getSplitReadTemplateNames().stream()).collect(Collectors.toList());
        final Set<String> evidenceTemplates = evidenceLinks.stream().map(link -> {
            final Collection<String> linkTemplates = new ArrayList<>(link.getReadPairTemplateNames());
            linkTemplates.addAll(link.getSplitReadTemplateNames());
            return linkTemplates;
        }).flatMap(Collection::stream).collect(Collectors.toSet());
        counterEvidenceTemplates.removeAll(evidenceTemplates);
        return counterEvidenceTemplates.size();
    }

    private boolean testHMMState(final List<Integer> states, final int numStates, final double minEventHMMCoverage) {
        return validStateFrequency(states, getValidHMMCopyStates(numStates)) >= minEventHMMCoverage * states.size();
    }

    private int validStateFrequency(final List<Integer> states, final Set<Integer> validStates) {
        return (int) states.stream().filter(validStates::contains).count();
    }

    protected static boolean segmentCallMatchesInterval(final CalledCopyRatioSegment segment, final SVInterval interval, final CalledCopyRatioSegment.Call callType, final SAMSequenceDictionary dictionary, final double minSegmentReciprocalOverlap) {
        final SVInterval segmentInterval = new SVInterval(dictionary.getSequenceIndex(segment.getContig()), segment.getStart(), segment.getEnd());
        return segment.getCall() == callType && IntervalUtils.hasReciprocalOverlap(segmentInterval, interval, minSegmentReciprocalOverlap);
    }
}
