package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.LargeSimpleSV;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.utils.IntrachromosomalBreakpointPair;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Collection;
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
                                         final SVIntervalTree<VariantContext> structuralVariantCallTree,
                                         final SVIntervalTree<GATKRead> contigTree,
                                         final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromReadDepthArgumentCollection arguments,
                                         final SAMSequenceDictionary dictionary) {
        super(intrachromosomalLinkTree, interchromosomalLinkTree, structuralVariantCallTree, contigTree, arguments, dictionary);
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
                                     final IntrachromosomalBreakpointPair breakpoints,
                                     final Collection<EvidenceTargetLink> supportingEvidence) {
        return new LargeSimpleSV(SimpleSVType.TYPES.DUP_TAND, start, end, contigId, readPairEvidence, splitReadEvidence, readPairCounterEvidence, splitReadCounterEvidence, breakpoints, supportingEvidence);
    }

    @Override
    protected int countSupportingEvidenceReadPairs(final Collection<EvidenceTargetLink> links) {
        return countReadPairs(links);
    }

    @Override
    protected int countSupportingEvidenceSplitReads(final Collection<EvidenceTargetLink> links) {
        return countSplitReads(links);
    }

    /**
     * Checks if the link matches "outie" read pair orientation, i.e. -/+
     */
    @Override
    protected boolean hasSupportingEvidenceOrientation(final EvidenceTargetLink link) {
        return hasOutieOrientation(link);
    }

    @Override
    protected boolean isCounterEvidenceOrientation(final EvidenceTargetLink link) {
        return !hasInnieOrientation(link) && !hasOutieOrientation(link);
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
        final int amplifiedBases = overlappingSegments.stream().filter(segment -> segment.getCall() == CalledCopyRatioSegment.Call.AMPLIFICATION)
                .mapToInt(segment -> SVIntervalUtils.convertToSVInterval(segment.getInterval(), dictionary).overlapLen(interval)).sum();
        return amplifiedBases / (double) interval.getLength() >= arguments.minSegmentOverlap;
    }

    /**
     * Tests if there is a suspicious number of copy ratio bins that are too high or too low
     */
    @Override
    protected boolean isInvalidCoverage(final List<CopyRatio> copyRatios) {
        if (!(copyRatios.stream().filter(ratio -> ratio.getLog2CopyRatioValue() < arguments.tandemDuplicationInvalidLog2CopyRatioThreshold).count() > arguments.tandemDuplicationInvalidBinFraction * copyRatios.size())) {
            return false;
        }
        return copyRatios.stream().anyMatch(ratio -> Math.pow(2.0, ratio.getLog2CopyRatioValue())  >= arguments.hmmMaxStates);
    }
}
