package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Set;

/**
 * Merge tail configurations:
 *
 * Performs the transformation:
 *
 * { A -> x + S_i + y }
 *
 * goes to:
 *
 * { A -> x -> S_i -> y }
 *
 * for all nodes that match this configuration.
 *
 * Differs from the diamond transform in that no bottom node is required
 */
final class MergeTails extends VertexBasedTransformer {

    /**
     * The minimum number of common bp from the prefix (head merging) or suffix (tail merging)
     * required before we'll merge in such configurations.  A large value here is critical to avoid
     * merging inappropriate head or tail nodes, which introduces large insertion / deletion events
     * as the merge operation creates a link among the non-linked sink / source vertices
     */
    static final int MIN_COMMON_SEQUENCE_TO_MERGE_SOURCE_SINK_VERTICES = 10;

    MergeTails(final SeqGraph graph) {
        super(graph);
    }

    @Override
    protected boolean tryToTransform(final SeqVertex top) {
        Utils.nonNull(top);

        final Set<SeqVertex> tails = getGraph().outgoingVerticesOf(top);
        if ( tails.size() <= 1 ) {
            return false;
        }

        for ( final SeqVertex t : tails ) {
            if (!getGraph().isSink(t) || getGraph().inDegreeOf(t) > 1) {
                return false;
            }
        }

        if ( dontModifyGraphEvenIfPossible() ) {
            return true;
        }

        final SharedVertexSequenceSplitter splitter = new SharedVertexSequenceSplitter(getGraph(), tails);

        return splitter.meetsMinMergableSequenceForSuffix(MIN_COMMON_SEQUENCE_TO_MERGE_SOURCE_SINK_VERTICES) && splitter.splitAndUpdate(top, null);
    }
}