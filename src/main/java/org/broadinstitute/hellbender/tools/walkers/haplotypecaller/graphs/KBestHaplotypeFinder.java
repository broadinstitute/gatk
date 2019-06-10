package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;
import java.util.Set;

/**
 * A common interface for the different KBestHaplotypeFinder implementations to conform to
 */
public abstract class KBestHaplotypeFinder<V extends BaseVertex, E extends BaseEdge> {
    protected final BaseGraph<V, E> graph;
    final Set<V> sinks;
    final Set<V> sources;

    public KBestHaplotypeFinder(final Set<V> sinks, final Set<V> sources, final BaseGraph<V, E> graph) {
        Utils.nonNull(graph, "graph cannot be null");
        Utils.nonNull(sources, "sources cannot be null");
        Utils.nonNull(sinks, "sinks cannot be null");
        Utils.validateArg(graph.containsAllVertices(sources), "source does not belong to the graph");
        Utils.validateArg(graph.containsAllVertices(sinks), "sink does not belong to the graph");

        this.sinks = sinks;
        this.sources = sources;
        //TODO dealing with cycles here due to a bug in some of the graph transformations that produces cycles.
        //TODO Once that is solve, the if-else below should be substituted by a throw if there is any cycles,
        //TODO just the line commented out below if you want to trade early-bug-fail for speed.
        this.graph = removeCyclesIfNecessary(graph, sources, sinks);
    }

    protected abstract BaseGraph<V, E> removeCyclesIfNecessary(BaseGraph<V, E> graph, Set<V> sources, Set<V> sinks);

    public abstract List<KBestHaplotype<V, E>> findBestHaplotypes(int maxNumberOfHaplotypes);

    @VisibleForTesting
    public List<KBestHaplotype<V, E>> findBestHaplotypes() {
        return findBestHaplotypes(Integer.MAX_VALUE);
    }
}
