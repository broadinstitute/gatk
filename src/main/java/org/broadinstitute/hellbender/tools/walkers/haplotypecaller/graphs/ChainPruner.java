package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import com.google.common.annotations.VisibleForTesting;

import java.util.*;


public abstract class ChainPruner<V extends BaseVertex, E extends BaseEdge> {
    public ChainPruner() { }

    public void pruneLowWeightChains(final BaseGraph<V,E> graph) {
        final List<Path<V, E>> chains = findAllChains(graph);
        final Collection<Path<V, E>> chainsToRemove = chainsToRemove(chains);
        chainsToRemove.forEach(c -> graph.removeAllEdges(c.getEdges()));
        graph.removeSingletonOrphanVertices();
    }

    @VisibleForTesting
    List<Path<V, E>> findAllChains(BaseGraph<V, E> graph) {
        final Deque<V> chainStarts = new LinkedList<>(graph.getSources());
        final List<Path<V,E>> chains = new LinkedList<>();
        final Set<V> alreadySeen = new HashSet<>(chainStarts);

        while(!chainStarts.isEmpty()) {
            final V chainStart = chainStarts.pop();
            for ( final E outEdge : graph.outgoingEdgesOf(chainStart) ) {
                final Path<V,E> chain = findChain(outEdge, graph);
                chains.add(chain);
                final V chainEnd = chain.getLastVertex();
                if (!alreadySeen.contains(chainEnd)) {
                    chainStarts.add(chainEnd);
                    alreadySeen.add(chainEnd);
                }
            }
        }
        return chains;
    }

    /**
     *
     * @return a fully extended linear path
     */
    private Path<V,E> findChain(final E startEdge, final BaseGraph<V, E> graph) {
        final List<E> edges = new ArrayList<>();
        edges.add(startEdge);
        final V firstVertex = graph.getEdgeSource(startEdge);
        V lastVertex = graph.getEdgeTarget(startEdge);

        while(true) {
            // chain ends if: 1) no out edges; 2) multiple out edges; 3) multiple in edges; 4) cycle back to start of chain
            final Set<E> outEdges = graph.outgoingEdgesOf(lastVertex);
            if (outEdges.size() != 1 || graph.inDegreeOf(lastVertex) > 1 || lastVertex.equals(firstVertex)) {
                break;
            }
            final E nextEdge = outEdges.iterator().next();
            edges.add(nextEdge);
            lastVertex = graph.getEdgeTarget(nextEdge);
        }
        return new Path<V, E>(edges, lastVertex, graph);
    }

    protected abstract Collection<Path<V,E>> chainsToRemove(final List<Path<V, E>> chains);
}
