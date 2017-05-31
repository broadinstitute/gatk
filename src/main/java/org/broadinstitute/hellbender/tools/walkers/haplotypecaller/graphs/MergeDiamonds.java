package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Set;

/**
 * Merge diamond configurations:
 *
 * Performance the transformation:
 *
 * { A -> x + S_i + y -> Z }
 *
 * goes to:
 *
 * { A -> x -> S_i -> y -> Z }
 *
 * for all nodes that match this configuration.
 */
final class MergeDiamonds extends VertexBasedTransformer {
    MergeDiamonds(final SeqGraph graph) {
        super(graph);
    }

    @Override
    protected boolean tryToTransform(final SeqVertex top) {
        Utils.nonNull(top);
        final Set<SeqVertex> middles = getGraph().outgoingVerticesOf(top);
        if ( middles.size() <= 1 )
        // we can only merge if there's at least two middle nodes
        {
            return false;
        }

        SeqVertex bottom = null;
        for ( final SeqVertex mi : middles ) {
            // all nodes must have at least 1 connection
            if ( getGraph().outDegreeOf(mi) < 1 ) {
                return false;
            }

            // can only have 1 incoming node, the root vertex
            if ( getGraph().inDegreeOf(mi) != 1 ) {
                return false;
            }

            // make sure that all outgoing vertices of mi go only to the bottom node
            for ( final SeqVertex mt : getGraph().outgoingVerticesOf(mi) ) {
                if ( bottom == null ) {
                    bottom = mt;
                } else if ( ! bottom.equals(mt) ) {
                    return false;
                }
            }
        }

        // bottom has some connections coming in from other nodes, don't allow
        if ( getGraph().inDegreeOf(bottom) != middles.size() ) {
            return false;
        }

        if ( dontModifyGraphEvenIfPossible() ) {
            return true;
        }

        // actually do the merging, returning true if at least 1 base was successfully split
        final SharedVertexSequenceSplitter splitter = new SharedVertexSequenceSplitter(getGraph(), middles);
        return splitter.meetsMinMergableSequenceForEitherPrefixOrSuffix(1) && splitter.splitAndUpdate(top, bottom);
    }
}