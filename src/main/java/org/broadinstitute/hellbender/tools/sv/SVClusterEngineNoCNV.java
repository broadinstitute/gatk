package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

public class SVClusterEngineNoCNV extends SVClusterEngine {

    public SVClusterEngineNoCNV(final SAMSequenceDictionary dictionary) {
        super(dictionary);
    }

    public SVClusterEngineNoCNV(final SAMSequenceDictionary dictionary, boolean depthOnly, BreakpointSummaryStrategy strategy) {
        super(dictionary, depthOnly, strategy);
    }

    @Override
    protected boolean clusterTogether(final SVCallRecord a, final SVCallRecord b) {
        if (!a.getType().equals(b.getType())) {
            return false;
        }
        final boolean depthOnlyA = isDepthOnlyCall(a);
        final boolean depthOnlyB = isDepthOnlyCall(b);
        if (depthOnlyA && depthOnlyB) {
            return clusterTogetherBothDepthOnly(a, b);
        } else if (depthOnlyA != depthOnlyB) {
            return clusterTogetherMixedEvidence(a, b);
        } else {
            return clusterTogetherBothWithEvidence(a, b);
        }
    }


    // TODO optimize intervals
    @Override
    protected SimpleInterval getClusteringInterval(final SVCallRecord call, final SimpleInterval clusterMinStartInterval) {
        final int padding = (int) Math.ceil(Math.max(Math.max(getEndpointClusteringPadding(call), call.getLength() * MIN_RECIPROCAL_OVERLAP_DEPTH), MIXED_CLUSTERING_WINDOW));
        final int minStart = call.getPositionA() - padding;
        final int maxStart = call.getPositionA() + padding;
        final String currentContig = getCurrentContig();
        if (clusterMinStartInterval == null) {
            return IntervalUtils.trimIntervalToContig(currentContig, minStart, maxStart, dictionary.getSequence(currentContig).getSequenceLength());
        }
        //NOTE: this is an approximation -- best method would back calculate cluster bounds, then rederive start and end based on call + cluster
        final int newMinStart = Math.min(minStart, clusterMinStartInterval.getStart());
        final int newMaxStart = Math.max(maxStart, clusterMinStartInterval.getEnd());
        return IntervalUtils.trimIntervalToContig(currentContig, newMinStart, newMaxStart, dictionary.getSequence(currentContig).getSequenceLength());
    }
}
