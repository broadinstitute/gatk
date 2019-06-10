package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.apache.commons.lang3.mutable.MutableInt;
import org.broadinstitute.hellbender.utils.Utils;
import org.jgrapht.alg.CycleDetector;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Efficient algorithm to obtain the list of best haplotypes given the {@link SeqGraph instace}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class GraphBasedKBestHaplotypeFinder<V extends BaseVertex, E extends BaseEdge> extends KBestHaplotypeFinder<V, E> {

    /**
     * Constructs a new best haplotypes finder.
     *
     * @param graph the seq-graph to search.
     * @param sources source vertices for all haplotypes.
     * @param sinks sink vertices for all haplotypes.
     *
     * @throws IllegalArgumentException if <ul>
     *     <li>any of {@code graph}, {@code sources} or {@code sinks} is {@code null} or</li>
     *     <li>any of {@code sources}' or any {@code sinks}' member is not a vertex in {@code graph}.</li>
     * </ul>
     */
    public GraphBasedKBestHaplotypeFinder(final BaseGraph<V, E> graph, final Set<V> sources, final Set<V> sinks) {
        super(sinks, sources, graph);
    }

    @Override
    protected BaseGraph<V, E> removeCyclesIfNecessary(BaseGraph<V, E> graph, Set<V> sources, Set<V> sinks) {
        return new CycleDetector<>(graph).detectCycles() ? removeCyclesAndVerticesThatDontLeadToSinks(graph,sources,sinks) : graph;
    }

    /**
     * Constructor for the special case of a single source and sink
     */
    public GraphBasedKBestHaplotypeFinder(final BaseGraph<V, E> graph, final V source, final V sink) {
        this(graph, Collections.singleton(source), Collections.singleton(sink));
    }

    /**
     * Constructor for the default case of all sources and sinks
     */
    public GraphBasedKBestHaplotypeFinder(final BaseGraph<V, E> graph) {
        this(graph, graph.getSources(), graph.getSinks());
    }

    /**
     * Implement Dijkstra's algorithm as described in https://en.wikipedia.org/wiki/K_shortest_path_routing
     */
    @Override
    public List<KBestHaplotype<V, E>> findBestHaplotypes(final int maxNumberOfHaplotypes) {
        final List<KBestHaplotype<V, E>> result = new ArrayList<>();
        final PriorityQueue<KBestHaplotype<V, E>> queue = new PriorityQueue<>(Comparator.comparingDouble(KBestHaplotype<V, E>::score).reversed());
        sources.forEach(source -> queue.add(new KBestHaplotype<>(source, graph)));

        final Map<V, MutableInt> vertexCounts = graph.vertexSet().stream()
                .collect(Collectors.toMap(v -> v, v -> new MutableInt(0)));

        while (!queue.isEmpty() && result.size() < maxNumberOfHaplotypes) {
            final KBestHaplotype<V, E> pathToExtend = queue.poll();
            final V vertexToExtend = pathToExtend.getLastVertex();
            if (sinks.contains(vertexToExtend)) {
                result.add(pathToExtend);
            } else {
                if (vertexCounts.get(vertexToExtend).getAndIncrement() < maxNumberOfHaplotypes) {
                    final Set<E> outgoingEdges = graph.outgoingEdgesOf(vertexToExtend);
                    int totalOutgoingMultiplicity = 0;
                    for (final BaseEdge edge : outgoingEdges) {
                        totalOutgoingMultiplicity += edge.getMultiplicity();
                    }

                    for (final E edge : outgoingEdges) {
                        final V targetVertex = graph.getEdgeTarget(edge);
                        queue.add(new KBestHaplotype<>(pathToExtend, edge, totalOutgoingMultiplicity));
                    }
                }
            }
        }
        return result;
    }

    /**
     * Removes edges that produces cycles and also dead vertices that do not lead to any sink vertex.
     * @return never {@code null}.
     */
    protected BaseGraph<V, E> removeCyclesAndVerticesThatDontLeadToSinks(final BaseGraph<V, E> original, final Collection<V> sources, final Set<V> sinks) {
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
}
