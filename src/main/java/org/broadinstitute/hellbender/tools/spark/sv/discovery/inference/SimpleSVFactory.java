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
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
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
    protected final ReadMetadata readMetadata;
    protected final SAMSequenceDictionary dict;

    public SimpleSVFactory(final SVIntervalTree<EvidenceTargetLink> intrachromosomalLinkTree,
                           final SVIntervalTree<EvidenceTargetLink> interchromosomalLinkTree,
                           final SVIntervalTree<GATKRead> contigTree,
                           final StructuralVariationDiscoveryArgumentCollection.SimpleVariantDiscoveryArgumentCollection arguments,
                           final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector,
                           final OverlapDetector<CopyRatio> readDepthOverlapDetector,
                           final ReadMetadata readMetadata,
                           final SAMSequenceDictionary dict) {
        this.intrachromosomalLinkTree = intrachromosomalLinkTree;
        this.interchromosomalLinkTree = interchromosomalLinkTree;
        this.contigTree = contigTree;
        this.arguments = arguments;
        this.copyRatioSegmentOverlapDetector = copyRatioSegmentOverlapDetector;
        this.readDepthOverlapDetector = readDepthOverlapDetector;
        this.readMetadata = readMetadata;
        this.dict = dict;
    }

    protected abstract boolean isEvidenceOrientation(final EvidenceTargetLink link);

    protected abstract boolean isInvalidCoverage(final List<CopyRatio> copyRatios);

    protected abstract Set<Integer> getValidHMMCopyStates(final int numStates);

    protected abstract boolean supportedBySegmentCalls(final SVInterval interval, final Set<CalledCopyRatioSegment> overlappingSegments);

    protected abstract int filterEvidence(final List<EvidenceTargetLink> links);

    protected abstract SimpleSV getNewSV(final int start,
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

    public SimpleSV create(final SVInterval leftInterval,
                           final SVInterval rightInterval,
                           final SVInterval callInterval,
                           final boolean applyFilter,
                           final IntrachromosomalBreakpointPair breakpoints) {

        final List<EvidenceTargetLink> breakpointOverlappingLinks = breakpointOverlappingLinks(leftInterval, rightInterval, intrachromosomalLinkTree, dict);
        final List<EvidenceTargetLink> evidenceLinks = getLinksWithEvidenceOrientation(breakpointOverlappingLinks);
        if (evidenceLinks.isEmpty()) return null;

        final int contigId = leftInterval.getContig();
        final String contig = readMetadata.getContigName(contigId);

        final SVInterval outerInterval = new SVInterval(contigId, leftInterval.getStart(), rightInterval.getEnd());
        final SVInterval innerInterval;
        if (rightInterval.getStart() > leftInterval.getEnd()) {
            innerInterval = new SVInterval(contigId, leftInterval.getEnd(), rightInterval.getStart());
        } else {
            innerInterval = new SVInterval(contigId, leftInterval.getEnd(), leftInterval.getEnd() + 1);
        }
        final List<EvidenceTargetLink> localOverlappingLinks = localOverlappingLinks(outerInterval, intrachromosomalLinkTree, arguments.SPANNING_COUNTEREVIDENCE_RANGE, dict);
        localOverlappingLinks.addAll(IntervalUtils.getOverlappingLinksOnInterval(outerInterval, interchromosomalLinkTree));
        final List<EvidenceTargetLink> counterEvidenceLinks = new ArrayList<>(CollectionUtils.subtract(localOverlappingLinks, evidenceLinks));

        final int filteredEvidence = filterEvidence(evidenceLinks);
        final int readPairEvidence = filteredEvidence; //readPairEvidence(filteredEvidenceLinks);
        final int splitReadEvidence = 0; //splitReadEvidence(filteredEvidenceLinks);
        final int readPairCounterEvidence = readPairCounterEvidence(counterEvidenceLinks, evidenceLinks, arguments.MIN_LINK_COUNTEREVIDENCE);
        final int splitReadCounterEvidence = splitReadCounterEvidence(counterEvidenceLinks, evidenceLinks, arguments.MIN_LINK_COUNTEREVIDENCE);
        final double evidenceToCounterEvidenceRatio = SimpleSV.computeScore(readPairEvidence, splitReadEvidence, readPairCounterEvidence, splitReadCounterEvidence, arguments.COUNTEREVIDENCE_PSEUDOCOUNT);

        if (evidenceToCounterEvidenceRatio < arguments.MIN_EVIDENCE_TO_COUNTEREVIDENCE_RATIO && applyFilter)
            return null;

        final Set<CalledCopyRatioSegment> overlappingSegments = copyRatioSegmentOverlapDetector.getOverlaps(new SimpleInterval(contig, outerInterval.getStart(), outerInterval.getEnd()));
        boolean supportedBySegmentCalls = supportedBySegmentCalls(outerInterval, overlappingSegments);
        final SVInterval hmmInterval = IntervalUtils.getPaddedInterval(innerInterval, arguments.HMM_PADDING, dict);
        final List<CopyRatio> copyRatioBins = IntervalUtils.getCopyRatiosOnInterval(hmmInterval, readDepthOverlapDetector, arguments.COPY_NUMBER_BIN_TRIMMING, contigId, dict);
        if (isInvalidCoverage(copyRatioBins) && applyFilter) return null;
        if (supportedBySegmentCalls && applyFilter) {
            return getNewSV(callInterval.getStart(), callInterval.getEnd(), contigId, contig, readPairEvidence, splitReadEvidence, readPairCounterEvidence, splitReadCounterEvidence,
                    Collections.emptyList(), Collections.emptyList(), breakpoints);
        }
        if (copyRatioBins.isEmpty()) {
            if (applyFilter) {
                return null;
            }
            return getNewSV(callInterval.getStart(), callInterval.getEnd(), contigId, contig, readPairEvidence, splitReadEvidence, readPairCounterEvidence, splitReadCounterEvidence,
                    Collections.emptyList(), Collections.emptyList(), breakpoints);
        }
        final List<Double> copyRatios = copyRatioBins.stream().map(CopyRatio::getLog2CopyRatioValue).collect(Collectors.toList());
        final int numStates = Math.min(2 * arguments.highDepthCoveragePeakFactor + 1, copyRatios.stream().mapToInt(val -> (int) Math.pow(2.0, val)).max().getAsInt() + 2);
        final RealVector copyNumberPriors = CopyNumberHMM.getCopyNumberPrior(numStates);
        final List<Integer> positionsList = CopyNumberHMM.getCopyNumberHMMPositions(copyRatios.size());
        final CopyNumberHMM copyNumberHMM = new CopyNumberHMM(copyNumberPriors, arguments.COPY_NUMBER_HMM_ALPHA);
        final List<Integer> copyNumberStates = ViterbiAlgorithm.apply(copyRatios, positionsList, copyNumberHMM);
        if (!applyFilter || testHMMState(copyNumberStates, numStates, arguments.MIN_EVENT_HMM_COVERAGE)) {
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

    private List<EvidenceTargetLink> breakpointOverlappingLinks(final SVInterval leftBreakpointInterval, final SVInterval rightBreakpointInterval, final SVIntervalTree<EvidenceTargetLink> tree, final SAMSequenceDictionary dict) {
        final List<EvidenceTargetLink> leftOverlappingLinks = IntervalUtils.getOverlappingLinksOnInterval(leftBreakpointInterval, tree);
        return leftOverlappingLinks.stream().filter(link -> IntervalUtils.linkEndsOverlapIntervals(leftBreakpointInterval, rightBreakpointInterval, link)).collect(Collectors.toList());
    }

    private List<EvidenceTargetLink> localOverlappingLinks(final SVInterval breakpointInterval, final SVIntervalTree<EvidenceTargetLink> tree, final int localRange, final SAMSequenceDictionary dict) {
        final List<EvidenceTargetLink> overlappingLinks = IntervalUtils.getOverlappingLinksOnInterval(breakpointInterval, tree);
        final SVInterval localInterval = IntervalUtils.getPaddedInterval(breakpointInterval, localRange, dict);
        return overlappingLinks.stream().filter(link -> {
            final SVInterval linkInterval = IntervalUtils.getLargeLinkEventInterval(link);
            final SVInterval leftInterval = link.getPairedStrandedIntervals().getLeft().getInterval();
            final SVInterval rightInterval = link.getPairedStrandedIntervals().getRight().getInterval();
            return IntervalUtils.containsInterval(localInterval, linkInterval) || breakpointInterval.overlaps(leftInterval) || breakpointInterval.overlaps(rightInterval);
        }).collect(Collectors.toList());
    }

    static int readPairCounterEvidence(final List<EvidenceTargetLink> counterEvidenceLinks, final List<EvidenceTargetLink> evidenceLinks, final int minEvidenceCount) {
        final List<String> counterEvidenceTemplates = counterEvidenceLinks.stream().filter(link -> link.getReadPairs() >= minEvidenceCount).flatMap(link -> link.getReadPairTemplateNames().stream()).collect(Collectors.toList());
        final Set<String> evidenceTemplates = evidenceLinks.stream().map(link -> {
            final Collection<String> linkTemplates = new ArrayList<>(link.getReadPairTemplateNames());
            linkTemplates.addAll(link.getSplitReadTemplateNames());
            return linkTemplates;
        }).flatMap(Collection::stream).collect(Collectors.toSet());
        counterEvidenceTemplates.removeAll(evidenceTemplates);
        return counterEvidenceTemplates.size();
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
}
