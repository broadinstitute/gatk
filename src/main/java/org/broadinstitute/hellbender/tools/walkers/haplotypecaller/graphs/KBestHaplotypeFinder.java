package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.jgrapht.alg.CycleDetector;

import java.util.*;

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

    private BaseGraph<V, E> removeCyclesIfNecessary(BaseGraph<V, E> graph, Set<V> sources, Set<V> sinks) {
        if (keepCycles()) {
            return graph;
        } else {
            return new CycleDetector<>(graph).detectCycles() ? removeCyclesAndVerticesThatDontLeadToSinks(graph, sources, sinks) : graph;
        }
    }

    /**
     * Removes edges that produces cycles and also dead vertices that do not lead to any sink vertex.
     * @return never {@code null}.
     */
    private BaseGraph<V, E> removeCyclesAndVerticesThatDontLeadToSinks(final BaseGraph<V, E> original, final Collection<V> sources, final Set<V> sinks) {
        final Set<E> edgesToRemove = new HashSet<>(original.edgeSet().size());
        final Set<V> vertexToRemove = new HashSet<>(original.vertexSet().size());

        boolean foundSomePath = false;
        for (final V source : sources) {
            final Set<V> parentVertices = new HashSet<>(original.vertexSet().size());
            foundSomePath = findGuiltyVerticesAndEdgesToRemoveCycles(original, source, sinks, edgesToRemove, vertexToRemove, parentVertices) || foundSomePath;
        }

        Utils.validate(foundSomePath, () -> "could not find any path from the source vertex to the sink vertex after removing cycles: "
                + Arrays.toString(sources.toArray()) + " => " + Arrays.toString(sinks.toArray()));

        Utils.validate(!(edgesToRemove.isEmpty() && vertexToRemove.isEmpty()), "cannot find a way to remove the cycles");

        final BaseGraph<V, E> result = original.clone();
        result.removeAllEdges(edgesToRemove);
        result.removeAllVertices(vertexToRemove);
        return result;
    }

    /**
     * Recursive call that looks for edges and vertices that need to be removed to get rid of cycles.
     *
     * @param graph the original graph.
     * @param currentVertex current search vertex.
     * @param sinks considered sink vertices.
     * @param edgesToRemove collection  of edges that need to be removed in order to get rid of cycles.
     * @param verticesToRemove collection of vertices that can be removed.
     * @param parentVertices collection of vertices that preceded the {@code currentVertex}; i.e. the it can be
     *                       reached from those vertices using edges existing in {@code graph}.
     *
     * @return {@code true} to indicate that the some sink vertex is reachable by {@code currentVertex},
     *  {@code false} otherwise.
     */
    private boolean findGuiltyVerticesAndEdgesToRemoveCycles(final BaseGraph<V, E> graph,
                                                             final V currentVertex,
                                                             final Set<V> sinks,
                                                             final Set<E> edgesToRemove,
                                                             final Set<V> verticesToRemove,
                                                             final Set<V> parentVertices) {
        if (sinks.contains(currentVertex)) {
            return true;
        }

        final Set<E> outgoingEdges = graph.outgoingEdgesOf(currentVertex);
        parentVertices.add(currentVertex);

        boolean reachesSink = false;
        for (final E edge : outgoingEdges) {
            final V child = graph.getEdgeTarget(edge);
            if (parentVertices.contains(child)) {
                edgesToRemove.add(edge);
            } else {
                final boolean childReachSink = findGuiltyVerticesAndEdgesToRemoveCycles(graph, child, sinks, edgesToRemove, verticesToRemove, parentVertices);
                reachesSink = reachesSink || childReachSink;
            }
        }
        if (!reachesSink) {
            verticesToRemove.add(currentVertex);
        }
        return reachesSink;
    }

    // Switch to be used in deciding whether or not to alter the graph for cycle safety.
    public abstract boolean keepCycles();

    public abstract List<KBestHaplotype<V, E>> findBestHaplotypes(int maxNumberOfHaplotypes);

    @VisibleForTesting
    public List<KBestHaplotype<V, E>> findBestHaplotypes() {
        return findBestHaplotypes(Integer.MAX_VALUE);
    }
}
