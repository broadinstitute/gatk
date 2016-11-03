package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import htsjdk.samtools.Cigar;
import joptsimple.internal.Strings;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A path thought a BaseGraph
 *
 * class to keep track of paths
 *
 */
public final class Path<T extends BaseVertex, E extends BaseEdge> {

    // the last vertex seen in the path
    private final T lastVertex;

    // the list of edges comprising the path
    private final List<E> edgesInOrder;

    // the scores for the path
    private final int totalScore;

    // the graph from which this path originated
    private final BaseGraph<T, E> graph;

    /**
     * Create a new Path containing no edges and starting at initialVertex
     * @param initialVertex the starting vertex of the path
     * @param graph the graph this path will follow through
     */
    public Path(final T initialVertex, final BaseGraph<T, E> graph) {
        Utils.nonNull(initialVertex, "initialVertex cannot be null");
        Utils.nonNull(graph, "graph cannot be null");
        Utils.validateArg(graph.containsVertex(initialVertex), () -> "Vertex " + initialVertex + " must be part of graph " + graph);

        lastVertex = initialVertex;
        edgesInOrder = new ArrayList<>(0);
        totalScore = 0;
        this.graph = graph;
    }

    /**
     * Create a new Path extending p with edge
     *
     * @param p the path to extend.
     * @param edge the edge to extend path with.
     *
     * @throws IllegalArgumentException if {@code p} or {@code edge} are {@code null}, or {@code edge} is
     * not part of {@code p}'s graph, or {@code edge} does not have as a source the last vertex in {@code p}.
     */
    public Path(final Path<T,E> p, final E edge) {
        Utils.nonNull(p, "Path cannot be null");
        Utils.nonNull(edge, "Edge cannot be null");
        Utils.validateArg(p.graph.containsEdge(edge), () -> "Graph must contain edge " + edge + " but it doesn't");
        if ( ! p.graph.getEdgeSource(edge).equals(p.lastVertex) ) { throw new IllegalStateException("Edges added to path must be contiguous."); }

        graph = p.graph;
        lastVertex = p.graph.getEdgeTarget(edge);
        edgesInOrder = new ArrayList<>(p.length() + 1);
        edgesInOrder.addAll(p.edgesInOrder);
        edgesInOrder.add(edge);
        totalScore = p.totalScore + edge.getMultiplicity();
    }

    /**
     * Length of the path in edges.
     *
     * @return {@code 0} or greater.
     */
    public int length() {
        return edgesInOrder.size();
    }

    /**
     * Prepend a path with an edge.
     *
     * @param edge the extending edge.
     * @param p the original path.
     *
     * @throws IllegalArgumentException if {@code p} or {@code edge} are {@code null}, or {@code edge} is
     * not part of {@code p}'s graph, or {@code edge} does not have as a target the first vertex in {@code p}.
     */
    public Path(final E edge, final Path<T,E> p) {
        Utils.nonNull(p, "Path cannot be null");
        Utils.nonNull(edge, "Edge cannot be null");
        Utils.validateArg(p.graph.containsEdge(edge), () -> "Graph must contain edge " + edge + " but it doesn't");
        if ( ! p.graph.getEdgeTarget(edge).equals(p.getFirstVertex())) { throw new IllegalStateException("Edges added to path must be contiguous."); }
        graph = p.graph;
        lastVertex = p.lastVertex;
        edgesInOrder = new ArrayList<>(p.length() + 1);
        edgesInOrder.add(edge);
        edgesInOrder.addAll(p.getEdges());
        totalScore = p.totalScore + edge.getMultiplicity();
    }

    /**
     * Check that two paths have the same edges and total score
     * @param path the other path we might be the same as
     * @return true if this and path are the same
     */
    public boolean pathsAreTheSame(final Path<T,E> path) {
        return totalScore == path.totalScore && edgesInOrder.equals(path.edgesInOrder);
    }

    /**
     * Does this path contain the given vertex?
     *
     * @param v a non-null vertex
     * @return true if v occurs within this path, false otherwise
     */
    public boolean containsVertex(final T v) {
        Utils.nonNull(v, "Vertex cannot be null");
        return v.equals(getFirstVertex()) || edgesInOrder.stream().map(graph::getEdgeTarget).anyMatch(v::equals);
    }

    @Override
    public String toString() {
        final String joinedPath = Strings.join(getVertices().stream().map(v -> v.getSequenceString()).collect(Collectors.toList()), "->");
        return String.format("Path{score=%d, path=%s}", totalScore, joinedPath);
    }

    /**
     * Get the graph of this path
     * @return a non-null graph
     */
    public BaseGraph<T, E> getGraph() {
        return graph;
    }

    /**
     * Get the edges of this path in order.
     * Returns an unmodifiable view of the underlying list
     * @return a non-null list of edges
     */
    public List<E> getEdges() { return Collections.unmodifiableList(edgesInOrder); }

    /**
     * Get the list of vertices in this path in order defined by the edges of the path
     * @return a non-null, non-empty list of vertices
     */
    public List<T> getVertices() {
        final List<T> result = new ArrayList<>(edgesInOrder.size()+1);
        result.add(getFirstVertex());
        result.addAll(edgesInOrder.stream().map(graph::getEdgeTarget).collect(Collectors.toList()));
        return result;
    }

    /**
     * Get the total score of this path (bigger is better)
     * @return a positive integer
     */
    public int getScore() { return totalScore; }

    /**
     * Get the final vertex of the path
     * @return a non-null vertex
     */
    public T getLastVertex() { return lastVertex; }

    /**
     * Get the first vertex in this path
     * @return a non-null vertex
     */
    public T getFirstVertex() {
        if (edgesInOrder.isEmpty()) {
            return lastVertex;
        } else {
            return getGraph().getEdgeSource(edgesInOrder.get(0));
        }
    }

    /**
     * The base sequence for this path. Pull the full sequence for source nodes and then the suffix for all subsequent nodes
     * @return  non-null sequence of bases corresponding to this path
     */
    public byte[] getBases() {
        if( getEdges().isEmpty() ) { return graph.getAdditionalSequence(lastVertex); }

        byte[] bases = graph.getAdditionalSequence(graph.getEdgeSource(edgesInOrder.get(0)));
        for( final E e : edgesInOrder ) {
            bases = ArrayUtils.addAll(bases, graph.getAdditionalSequence(graph.getEdgeTarget(e)));
        }
        return bases;
    }


    /**
     * Calculate the cigar elements for this path against the reference sequence
     *
     * @param refSeq the reference sequence that all of the bases in this path should align to
     * @return a Cigar mapping this path to refSeq, or null if no reasonable alignment could be found
     */
    public  Cigar calculateCigar(final byte[] refSeq) {
        //Note: CigarUtils.calculateCigar already checks for null
        return CigarUtils.calculateCigar(refSeq,getBases());
    }

}
