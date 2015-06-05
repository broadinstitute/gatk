package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collection;
import java.util.HashSet;

/**
 * Performs the transformation:
 *
 * { x + S_i + y -> Z }
 *
 * goes to:
 *
 * { x -> S_i -> y -> Z }
 *
 * for all nodes that match this configuration.
 *
 * Differs from the diamond transform in that no top node is required
 */
final class SplitCommonSuffices extends VertexBasedTransformer {
    private final Collection<SeqVertex> alreadySplit = new HashSet<>();

    SplitCommonSuffices(final SeqGraph graph) {
        super(graph);
    }

    @Override
    boolean tryToTransform(final SeqVertex bottom) {
        Utils.nonNull(bottom);

        if (alreadySplit.contains(bottom)) {
            return false;
        } else {
            alreadySplit.add(bottom);
            return CommonSuffixSplitter.split(getGraph(), bottom);
        }
    }
}
