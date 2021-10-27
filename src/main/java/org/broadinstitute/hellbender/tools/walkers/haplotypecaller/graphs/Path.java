package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import joptsimple.internal.Strings;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;

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
public class Path<V extends BaseVertex, E extends BaseEdge> {

    // the last vertex seen in the path
    private final V lastVertex;

    // the list of edges comprising the path
    private final List<E> edgesInOrder;

    // the graph from which this path originated
    private final BaseGraph<V, E> graph;

    /**
     * Create a new Path containing no edges and starting at initialVertex
     * @param initialVertex the starting vertex of the path
     * @param graph the graph this path will follow through
     */
    public Path(final V initialVertex, final BaseGraph<V, E> graph) {
        lastVertex = Utils.nonNull(initialVertex, "initialVertex cannot be null");
        this.graph = Utils.nonNull(graph, "graph cannot be null");
        Utils.validateArg(graph.containsVertex(initialVertex), () -> "Vertex " + initialVertex + " must be part of graph " + graph);

        edgesInOrder = new ArrayList<>(0);
    }

    /**
     * Constructor that does not check arguments' validity i.e. doesn't check that edges are in order
     */
    public Path(final List<E> edgesInOrder, final V lastVertex, final BaseGraph<V,E> graph) {
        this.lastVertex = lastVertex;
        this.graph = graph;
        this.edgesInOrder = edgesInOrder;
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
    public Path(final Path<V,E> p, final E edge) {
        Utils.nonNull(p, "Path cannot be null");
        Utils.nonNull(edge, "Edge cannot be null");
        Utils.validateArg(p.graph.containsEdge(edge), () -> "Graph must contain edge " + edge + " but it doesn't");
        Utils.validate( p.graph.getEdgeSource(edge).equals(p.lastVertex), "Edges added to path must be contiguous.");

        graph = p.graph;
        lastVertex = p.graph.getEdgeTarget(edge);
        edgesInOrder = new ArrayList<>(p.length() + 1);
        edgesInOrder.addAll(p.edgesInOrder);
        edgesInOrder.add(edge);
    }

    /**
     * Create a new Path extending p with edge
     *
     * @param p the path to extend.
     * @param edges list of edges to extend. Does not check arguments' validity i.e. doesn't check that edges are in order
     *
     * @throws IllegalArgumentException if {@code p} or {@code edges} are {@code null} or empty, or {@code edges} is
     * not part of {@code p}'s graph, or {@code edges} does not have as a source the last vertex in {@code p}.
     */
    public Path(final Path<V,E> p, final List<E> edges) {
        Utils.nonNull(p, "Path cannot be null");
        Utils.nonEmpty(edges, "Edge cannot be null");
        edges.forEach(edge -> Utils.validateArg(p.graph.containsEdge(edge), () -> "Graph must contain edge " + edge + " but it doesn't"));
        // Sanity check that the provided path is contiguous
        V tmpVertex = p.lastVertex;
        for (int i = 0; i < edges.size(); i++) {
            if ( ! p.graph.getEdgeSource(edges.get(i)).equals(tmpVertex) ) {
                throw new IllegalStateException("Edges added to path must be contiguous.");
            }
            tmpVertex = p.graph.getEdgeTarget(edges.get(i));
        }

        graph = p.graph;
        lastVertex = tmpVertex;
        edgesInOrder = new ArrayList<>(p.length() + 1);
        edgesInOrder.addAll(p.edgesInOrder);
        edgesInOrder.addAll(edges);
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
    public Path(final E edge, final Path<V,E> p) {
        Utils.nonNull(p, "Path cannot be null");
        Utils.nonNull(edge, "Edge cannot be null");
        Utils.validateArg(p.graph.containsEdge(edge), () -> "Graph must contain edge " + edge + " but it doesn't");
        if ( ! p.graph.getEdgeTarget(edge).equals(p.getFirstVertex())) { throw new IllegalStateException("Edges added to path must be contiguous."); }
        graph = p.graph;
        lastVertex = p.lastVertex;
        edgesInOrder = new ArrayList<>(p.length() + 1);
        edgesInOrder.add(edge);
        edgesInOrder.addAll(p.getEdges());
    }

    @VisibleForTesting
    boolean pathsAreTheSame(final Path<V,E> path) {
        return edgesInOrder.equals(path.edgesInOrder);
    }

    /**
     * Does this path contain the given vertex?
     *
     * @param v a non-null vertex
     * @return true if v occurs within this path, false otherwise
     */
    public boolean containsVertex(final V v) {
        Utils.nonNull(v, "Vertex cannot be null");
        return v.equals(getFirstVertex()) || edgesInOrder.stream().map(graph::getEdgeTarget).anyMatch(v::equals);
    }

    @Override
    public String toString() {
        final String joinedPath = Strings.join(getVertices().stream().map(v -> v.getSequenceString()).collect(Collectors.toList()), "->");
        return String.format("Path{path=%s}", joinedPath);
    }

    /**
     * Get the graph of this path
     * @return a non-null graph
     */
    public BaseGraph<V, E> getGraph() {
        return graph;
    }

    /**
     * Get the edges of this path in order.
     * Returns an unmodifiable view of the underlying list
     * @return a non-null list of edges
     */
    public List<E> getEdges() { return Collections.unmodifiableList(edgesInOrder); }

    public E getLastEdge() { return edgesInOrder.get(length() - 1); }

    /**
     * Get the list of vertices in this path in order defined by the edges of the path
     * @return a non-null, non-empty list of vertices
     */
    public List<V> getVertices() {
        final List<V> result = new ArrayList<>(edgesInOrder.size()+1);
        result.add(getFirstVertex());
        result.addAll(edgesInOrder.stream().map(graph::getEdgeTarget).collect(Collectors.toList()));
        return result;
    }

    /**
     * Get the final vertex of the path
     * @return a non-null vertex
     */
    public V getLastVertex() { return lastVertex; }

    /**
     * Get the first vertex in this path
     * @return a non-null vertex
     */
    public V getFirstVertex() {
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
        if( getEdges().isEmpty() ) { return BaseGraph.getAdditionalSequence(lastVertex, true); }

        byte[] bases = BaseGraph.getAdditionalSequence(graph.getEdgeSource(edgesInOrder.get(0)), true);
        for( final E e : edgesInOrder ) {
            bases = ArrayUtils.addAll(bases, BaseGraph.getAdditionalSequence(graph.getEdgeTarget(e), false));
        }
        return bases;
    }


    /**
     * Calculate the cigar elements for this path against the reference sequence
     *
     * @param refSeq the reference sequence that all of the bases in this path should align to
     * @return a Cigar mapping this path to refSeq, or null if no reasonable alignment could be found
     */
    public  Cigar calculateCigar(final byte[] refSeq, final SmithWatermanAligner aligner, final SWParameters pathToReferenceSWParameters) {
        //Note: CigarUtils.calculateCigar already checks for null
        return CigarUtils.calculateCigar(refSeq, getBases(), aligner, pathToReferenceSWParameters, SWOverhangStrategy.SOFTCLIP);
    }

}
