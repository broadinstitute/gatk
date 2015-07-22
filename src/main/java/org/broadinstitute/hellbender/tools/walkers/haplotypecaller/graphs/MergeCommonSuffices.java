package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.Utils;

/**
 * Merge headless configurations:
 *
 * Performs the transformation:
 *
 * { x + S_i -> y -> Z }
 *
 * goes to:
 *
 * { x -> S_i -> y + Z }
 *
 * for all nodes that match this configuration.
 */
final class MergeCommonSuffices extends VertexBasedTransformer {
    MergeCommonSuffices(final SeqGraph graph) {
        super(graph);
    }

    @Override
    boolean tryToTransform(final SeqVertex bottom) {
        Utils.nonNull(bottom);
        return SharedSequenceMerger.merge(getGraph(), bottom);
    }
}
