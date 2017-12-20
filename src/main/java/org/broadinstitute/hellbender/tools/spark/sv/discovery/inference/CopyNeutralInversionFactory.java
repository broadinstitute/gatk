package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.IntrachromosomalBreakpointPair;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.utils.PairedStrandedIntervals;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class CopyNeutralInversionFactory extends SimpleSVFactory {

    public CopyNeutralInversionFactory(final SVIntervalTree<EvidenceTargetLink> intrachromosomalLinkTree,
                                       final SVIntervalTree<EvidenceTargetLink> interchromosomalLinkTree,
                                       final SVIntervalTree<GATKRead> contigTree,
                                       final StructuralVariationDiscoveryArgumentCollection.SimpleVariantDiscoveryArgumentCollection arguments,
                                       final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector,
                                       final OverlapDetector<CopyRatio> readDepthOverlapDetector,
                                       final ReadMetadata readMetadata,
                                       final SAMSequenceDictionary dict) {
        super(intrachromosomalLinkTree, interchromosomalLinkTree, contigTree, arguments, copyRatioSegmentOverlapDetector, readDepthOverlapDetector, readMetadata, dict);
    }

    @Override
    protected SimpleSV getNewSV(final int start,
                                final int end,
                                final int contigId,
                                final String contig,
                                final int readPairEvidence,
                                final int splitReadEvidence,
                                final int readPairCounterEvidence,
                                final int splitReadCounterEvidence,
                                final List<CopyRatio> coverage,
                                final List<Integer> copyNumberStates,
                                final IntrachromosomalBreakpointPair breakpoints) {
        return new SimpleSV(SimpleSVType.TYPES.INV, start, end, contigId, contig, readPairEvidence, splitReadEvidence, readPairCounterEvidence, splitReadCounterEvidence,
                coverage, copyNumberStates, breakpoints);
    }

    @Override
    protected boolean isEvidenceOrientation(final EvidenceTargetLink link) {
        final PairedStrandedIntervals intervals = link.getPairedStrandedIntervals();
        return intervals.getLeft().getStrand() == intervals.getRight().getStrand();
    }

    @Override
    protected Set<Integer> getValidHMMCopyStates(final int numStates) {
        return IntStream.range(2, 3).boxed().collect(Collectors.toSet());
    }

    @Override
    protected boolean supportedBySegmentCalls(final SVInterval interval, final Set<CalledCopyRatioSegment> overlappingSegments) {
        return overlappingSegments.stream().anyMatch(segment -> segmentCallMatchesInterval(segment, interval));
    }

    @Override
    protected boolean isInvalidCoverage(final List<CopyRatio> copyRatios) {
        if (!(copyRatios.stream().filter(ratio -> ratio.getLog2CopyRatioValue() < arguments.NON_DELETION_INVALID_COPY_RATIO_THRESHOLD).count() > arguments.NON_DELETION_INVALID_BIN_FRACTION * copyRatios.size())) {
            return false;
        }
        return copyRatios.stream().anyMatch(ratio -> Math.pow(2.0, ratio.getLog2CopyRatioValue())  >= arguments.highDepthCoveragePeakFactor);
    }

    public boolean segmentCallMatchesInterval(final CalledCopyRatioSegment segment, final SVInterval interval) {
        final SVInterval segmentInterval = new SVInterval(readMetadata.getContigID(segment.getContig()), segment.getStart(), segment.getEnd());
        return segment.getCall() == CalledCopyRatioSegment.Call.NEUTRAL && interval.overlapLen(segmentInterval) >= arguments.MIN_COPY_NEUTRAL_SV_SEGMENT_OVERLAP * interval.getLength();
    }

    @Override
    protected int filterEvidence(final List<EvidenceTargetLink> links) {
        final List<EvidenceTargetLink> righties = links.stream().filter(link -> link.getPairedStrandedIntervals().getLeft().getStrand() && link.getPairedStrandedIntervals().getRight().getStrand()).collect(Collectors.toList());
        final List<EvidenceTargetLink> lefties = links.stream().filter(link -> !link.getPairedStrandedIntervals().getLeft().getStrand() && !link.getPairedStrandedIntervals().getRight().getStrand()).collect(Collectors.toList());
        final int numRightyEvidence = righties.stream().mapToInt(link -> link.getSplitReads() + link.getReadPairs()).sum();
        final int numLeftyEvidence = lefties.stream().mapToInt(link -> link.getSplitReads() + link.getReadPairs()).sum();
        final int minEvidenceType = Math.min(numRightyEvidence, numLeftyEvidence);
        return minEvidenceType * 2;
    }
}
