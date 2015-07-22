package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.Utils;

/**
 * Base class for transformation operations that need to iterate over proposed vertices, where
 * each proposed vertex is a seed vertex for a potential transformation.
 *
 * transformUntilComplete will iteratively apply the tryToTransform function on each vertex in the graph
 * until no vertex can be found that can be transformed.
 *
 * Note that in order to eventually terminate tryToTransform must transform the graph such that eventually
 * no vertices are candidates for further transformations.
 */
abstract class VertexBasedTransformer {
    /**
     * For testing purposes we sometimes want to test that can be transformed capabilities are working
     * without actually modifying the graph */
    private boolean dontModifyGraphEvenIfPossible = false;

    public boolean dontModifyGraphEvenIfPossible() { return dontModifyGraphEvenIfPossible; }
    public void setDontModifyGraphEvenIfPossible() { this.dontModifyGraphEvenIfPossible = true; }

    private final SeqGraph graph;

    VertexBasedTransformer(final SeqGraph graph){
        Utils.nonNull(graph);
        this.graph= graph;
    }

    SeqGraph getGraph() {
        return graph;
    }

    /**
     * Merge until the graph has no vertices that are candidates for merging
     */
    public boolean transformUntilComplete() {
        boolean didAtLeastOneTransform = false;
        boolean foundNodesToMerge = true;
        while( foundNodesToMerge ) {
            foundNodesToMerge = false;

            for( final SeqVertex v : graph.vertexSet() ) {
                foundNodesToMerge = tryToTransform(v);
                if ( foundNodesToMerge ) {
                    didAtLeastOneTransform = true;
                    break;
                }
            }
        }

        return didAtLeastOneTransform;
    }

    /**
     * Merge, if possible, seeded on the vertex v
     * @param v the proposed seed vertex to merge
     * @return true if some useful merging happened, false otherwise
     */
    abstract boolean tryToTransform(final SeqVertex v);
}
