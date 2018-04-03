package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.utils.IntrachromosomalBreakpointPair;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalUtils;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.LargeSimpleSV;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.utils.PairedStrandedIntervals;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Calls large tandem duplication variants
 */
public class LargeTandemDuplicationFactory extends LargeSimpleSVFactory {

    public LargeTandemDuplicationFactory(final SVIntervalTree<EvidenceTargetLink> intrachromosomalLinkTree,
                                         final SVIntervalTree<EvidenceTargetLink> interchromosomalLinkTree,
                                         final SVIntervalTree<GATKRead> contigTree,
                                         final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromReadDepthArgumentCollection arguments,
                                         final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector,
                                         final OverlapDetector<CopyRatio> copyRatioOverlapDetector,
                                         final SAMSequenceDictionary dictionary) {
        super(intrachromosomalLinkTree, interchromosomalLinkTree, contigTree, arguments, copyRatioSegmentOverlapDetector, copyRatioOverlapDetector, dictionary);
    }

    @Override
    protected LargeSimpleSV getNewSV(final int start,
                                     final int end,
                                     final int contigId,
                                     final String contig,
                                     final int readPairEvidence,
                                     final int splitReadEvidence,
                                     final int readPairCounterEvidence,
                                     final int splitReadCounterEvidence,
                                     final IntrachromosomalBreakpointPair breakpoints) {
        return new LargeSimpleSV(SimpleSVType.TYPES.DUP, start, end, contigId, readPairEvidence, splitReadEvidence, readPairCounterEvidence, splitReadCounterEvidence, breakpoints);
    }

    /**
     * Checks if the link matches "outie" read pair orientation, i.e. -/+
     */
    @Override
    protected boolean hasSupportingEvidenceOrientation(final EvidenceTargetLink link) {
        final PairedStrandedIntervals intervals = link.getPairedStrandedIntervals();
        return !intervals.getLeft().getStrand() && intervals.getRight().getStrand();
    }

    /**
     * Returns elevated copy number states (3, 4, ...)
     */
    @Override
    protected Set<Integer> getValidHMMCopyStates(final int numStates) {
        return IntStream.range(3, numStates).boxed().collect(Collectors.toSet());
    }

    /**
     * Tests if an amplification call matches the interval
     */
    @Override
    protected boolean supportedBySegmentCalls(final SVInterval interval, final Set<CalledCopyRatioSegment> overlappingSegments, final SAMSequenceDictionary dictionary) {
        //return overlappingSegments.stream().filter(call -> call.getCall() == CalledCopyRatioSegment.Call.AMPLIFICATION).map(CalledCopyRatioSegment::getInterval).map(segmentInterval -> SVIntervalUtils.convertToSVInterval(segmentInterval, dictionary)).anyMatch(segment -> SVIntervalUtils.hasReciprocalOverlap(interval, segment, arguments.minSegmentOverlap));
        //if (overlappingSegments.stream().anyMatch(segment -> segment.getMeanLog2CopyRatio() > 3)) return false;
        final int amplifiedBases = overlappingSegments.stream().filter(segment -> segment.getCall() == CalledCopyRatioSegment.Call.AMPLIFICATION)
                .mapToInt(segment -> SVIntervalUtils.convertToSVInterval(segment.getInterval(), dictionary).overlapLen(interval)).sum();
        //final int overlappingSegmentsLength = overlappingSegments.stream().filter(segment -> segment.getCall() == CalledCopyRatioSegment.Call.AMPLIFICATION).mapToInt(call -> call.getLengthOnReference()).sum();
        return amplifiedBases / (double) interval.getLength() >= arguments.minSegmentOverlap; // && amplifiedBases / (double) overlappingSegmentsLength >= arguments.minSegmentOverlap;
    }

    /**
     * Tests if there is a suspicious number of copy ratio bins that are too high or too low
     */
    @Override
    protected boolean isInvalidCoverage(final List<CopyRatio> copyRatios) {
        if (copyRatios.stream().anyMatch(copyRatio -> copyRatio.getLog2CopyRatioValue() > Math.log(6.0) / Math.log(2.0))) return false;
        if (!(copyRatios.stream().filter(ratio -> ratio.getLog2CopyRatioValue() < arguments.tandemDuplicationInvalidLog2CopyRatioThreshold).count() > arguments.tandemDuplicationInvalidBinFraction * copyRatios.size())) {
            return false;
        }
        return copyRatios.stream().anyMatch(ratio -> Math.pow(2.0, ratio.getLog2CopyRatioValue())  >= arguments.hmmMaxStates);
    }
}
