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
public final class KBestHaplotypeFinder {

    private final SeqGraph graph;
    final Set<SeqVertex> sinks;
    final Set<SeqVertex> sources;

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
    public KBestHaplotypeFinder(final SeqGraph graph, final Set<SeqVertex> sources, final Set<SeqVertex> sinks) {
        Utils.nonNull(graph, "graph cannot be null");
        Utils.nonNull(sources, "sources cannot be null");
        Utils.nonNull(sinks, "sinks cannot be null");
        Utils.validateArg(graph.containsAllVertices(sources), "source does not belong to the graph");
        Utils.validateArg(graph.containsAllVertices(sinks), "sink does not belong to the graph");

        //TODO dealing with cycles here due to a bug in some of the graph transformations that produces cycles.
        //TODO Once that is solve, the if-else below should be substituted by a throw if there is any cycles,
        //TODO just the line commented out below if you want to trade early-bug-fail for speed.
        //this.graph = graph;
        this.graph = new CycleDetector<>(graph).detectCycles() ? removeCyclesAndVerticesThatDontLeadToSinks(graph,sources,sinks) : graph;

        this.sinks = sinks;
        this.sources = sources;
    }

    /**
     * Constructor for the special case of a single source and sink
     */
    public KBestHaplotypeFinder(final SeqGraph graph, final SeqVertex source, final SeqVertex sink) {
        this(graph, Collections.singleton(source), Collections.singleton(sink));
    }

    /**
     * Constructor for the default case of all sources and sinks
     */
    public KBestHaplotypeFinder(final SeqGraph graph) {
        this(graph, graph.getSources(), graph.getSinks());
    }

    /**
     * Implement Dijkstra's algorithm as described in https://en.wikipedia.org/wiki/K_shortest_path_routing
     */
    public List<KBestHaplotype> findBestHaplotypes(final int maxNumberOfHaplotypes) {
        final List<KBestHaplotype> result = new ArrayList<>();
        final PriorityQueue<KBestHaplotype> queue = new PriorityQueue<>(Comparator.comparingDouble(KBestHaplotype::score).reversed());
        sources.forEach(source -> queue.add(new KBestHaplotype(source, graph)));

        final Map<SeqVertex, MutableInt> vertexCounts = graph.vertexSet().stream()
                .collect(Collectors.toMap(v -> v, v -> new MutableInt(0)));

        while (!queue.isEmpty() && result.size() < maxNumberOfHaplotypes) {
            final KBestHaplotype pathToExtend = queue.poll();
            final SeqVertex vertexToExtend = pathToExtend.getLastVertex();
            if (sinks.contains(vertexToExtend)) {
                result.add(pathToExtend);
            } else {
                final Set<BaseEdge> outgoingEdges = graph.outgoingEdgesOf(vertexToExtend);
                int totalOutgoingMultiplicity = 0;
                for (final BaseEdge edge : outgoingEdges) {
                    totalOutgoingMultiplicity += edge.getMultiplicity();
                }

                for (final BaseEdge edge : outgoingEdges) {
                    final SeqVertex targetVertex = graph.getEdgeTarget(edge);
                    if (vertexCounts.get(targetVertex).getAndIncrement() < maxNumberOfHaplotypes) {
                        queue.add(new KBestHaplotype(pathToExtend, edge, totalOutgoingMultiplicity));
                    }
                }
            }
        }
        return result;
    }

    public List<KBestHaplotype> findBestHaplotypes() {
       return findBestHaplotypes(Integer.MAX_VALUE);
    }

    /**
     * Removes edges that produces cycles and also dead vertices that do not lead to any sink vertex.
     * @return never {@code null}.
     */
    private static SeqGraph removeCyclesAndVerticesThatDontLeadToSinks(final SeqGraph original, final Collection<SeqVertex> sources, final Set<SeqVertex> sinks) {
        final Set<BaseEdge> edgesToRemove = new HashSet<>(original.edgeSet().size());
        final Set<SeqVertex> vertexToRemove = new HashSet<>(original.vertexSet().size());

        boolean foundSomePath = false;
        for (final SeqVertex source : sources) {
            final Set<SeqVertex> parentVertices = new HashSet<>(original.vertexSet().size());
            foundSomePath = findGuiltyVerticesAndEdgesToRemoveCycles(original, source, sinks, edgesToRemove, vertexToRemove, parentVertices) || foundSomePath;
        }

        Utils.validate(foundSomePath, () -> "could not find any path from the source vertex to the sink vertex after removing cycles: "
                    + Arrays.toString(sources.toArray()) + " => " + Arrays.toString(sinks.toArray()));

        Utils.validate(!(edgesToRemove.isEmpty() && vertexToRemove.isEmpty()), "cannot find a way to remove the cycles");

        final SeqGraph result = original.clone();
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
    private static boolean findGuiltyVerticesAndEdgesToRemoveCycles(final SeqGraph graph,
                                                                    final SeqVertex currentVertex,
                                                                    final Set<SeqVertex> sinks,
                                                                    final Set<BaseEdge> edgesToRemove,
                                                                    final Set<SeqVertex> verticesToRemove,
                                                                    final Set<SeqVertex> parentVertices) {
        if (sinks.contains(currentVertex)) {
            return true;
        }

        final Set<BaseEdge> outgoingEdges = graph.outgoingEdgesOf(currentVertex);
        parentVertices.add(currentVertex);

        boolean reachesSink = false;
        for (final BaseEdge edge : outgoingEdges) {
            final SeqVertex child = graph.getEdgeTarget(edge);
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
