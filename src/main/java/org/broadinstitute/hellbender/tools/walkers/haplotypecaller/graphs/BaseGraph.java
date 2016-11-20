package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import com.google.common.collect.Sets;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.jgrapht.EdgeFactory;
import org.jgrapht.alg.CycleDetector;
import org.jgrapht.graph.DefaultDirectedGraph;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;
import java.util.function.Supplier;
import java.util.stream.Collectors;

/**
 * Common code for graphs used for local assembly.
 */
public abstract class BaseGraph<V extends BaseVertex, E extends BaseEdge> extends DefaultDirectedGraph<V, E> {
    private static final long serialVersionUID = 1l;
    protected final int kmerSize;

    /**
     * Construct a TestGraph with kmerSize
     * @param kmerSize
     */
    protected BaseGraph(final int kmerSize, final EdgeFactory<V,E> edgeFactory) {
        super(edgeFactory);
        Utils.validateArg(kmerSize > 0, () -> "kmerSize must be > 0 but got " + kmerSize);
        this.kmerSize = kmerSize;
    }

    /**
     * How big of a kmer did we use to create this graph?
     * @return
     */
    public final int getKmerSize() {
        return kmerSize;
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a reference node (meaning that it appears on the reference path in the graph)
     */
    public final boolean isReferenceNode( final V v ) {
        Utils.nonNull(v, "Attempting to test a null vertex.");

        if (edgesOf(v).stream().anyMatch(e -> e.isRef())){
            return true;
        }

        // edge case: if the graph only has one node then it's a ref node, otherwise it's not
        return vertexSet().size() == 1;
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a source node (in degree == 0)
     */
    public final boolean isSource( final V v ) {
        Utils.nonNull(v, "Attempting to test a null vertex.");
        return inDegreeOf(v) == 0;
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a sink node (out degree == 0)
     */
    public final boolean isSink( final V v ) {
        Utils.nonNull(v, "Attempting to test a null vertex.");
        return outDegreeOf(v) == 0;
    }

    /**
     * Get the set of source vertices of this graph
     * @return a non-null set
     */
    public final Set<V> getSources() {
        return vertexSet().stream().filter(v -> isSource(v)).collect(Collectors.toSet());
    }

    /**
     * Get the set of sink vertices of this graph
     * @return a non-null set
     */
    public final Set<V> getSinks() {
        return vertexSet().stream().filter(v -> isSink(v)).collect(Collectors.toSet());
    }

    /**
     * Convert this kmer graph to a simple sequence graph.
     *
     * Each kmer suffix shows up as a distinct SeqVertex, attached in the same structure as in the kmer
     * graph.  Nodes that are sources are mapped to SeqVertex nodes that contain all of their sequence
     *
     * @return a newly allocated SequenceGraph
     */
    public SeqGraph toSequenceGraph() {
        final SeqGraph seqGraph = new SeqGraph(kmerSize);
        final Map<V, SeqVertex> vertexMap = new HashMap<>();

        // create all of the equivalent seq graph vertices
        for ( final V dv : vertexSet() ) {
            final SeqVertex sv = new SeqVertex(dv.getAdditionalSequence(isSource(dv)));
            sv.setAdditionalInfo(dv.getAdditionalInfo());
            vertexMap.put(dv, sv);
            seqGraph.addVertex(sv);
        }

        // walk through the nodes and connect them to their equivalent seq vertices
        for( final E e : edgeSet() ) {
            final SeqVertex seqInV = vertexMap.get(getEdgeSource(e));
            final SeqVertex seqOutV = vertexMap.get(getEdgeTarget(e));
            seqGraph.addEdge(seqInV, seqOutV, new BaseEdge(e.isRef(), e.getMultiplicity()));
        }

        return seqGraph;
    }

    /**
     * Pull out the additional sequence implied by traversing this node in the graph
     * @param v the vertex from which to pull out the additional base sequence
     * @return  non-null byte array
     */
    public final byte[] getAdditionalSequence( final V v ) {
        Utils.nonNull(v, "Attempting to pull sequence from a null vertex.");
        return v.getAdditionalSequence(isSource(v));
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a reference source
     */
    public final boolean isRefSource( final V v ) {
        Utils.nonNull(v, "Attempting to pull sequence from a null vertex.");

        // confirm that no incoming edges are reference edges
        if (incomingEdgesOf(v).stream().anyMatch(e -> e.isRef())) {
            return false;
        }

        // confirm that there is an outgoing reference edge
        if (outgoingEdgesOf(v).stream().anyMatch(e -> e.isRef())) {
            return true;
        }

        // edge case: if the graph only has one node then it's a ref source, otherwise it's not
        return vertexSet().size() == 1;
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a reference sink
     */
    public final boolean isRefSink( final V v ) {
        Utils.nonNull(v, "Attempting to pull sequence from a null vertex.");

        // confirm that no outgoing edges are reference edges
        if (outgoingEdgesOf(v).stream().anyMatch(e -> e.isRef())) {
            return false;
        }

        // confirm that there is an incoming reference edge
        if (incomingEdgesOf(v).stream().anyMatch(e -> e.isRef())) {
            return true;
        }

        // edge case: if the graph only has one node then it's a ref sink, otherwise it's not
        return vertexSet().size() == 1;
    }

    /**
     * @return the reference source vertex pulled from the graph, can be null if it doesn't exist in the graph
     */
    public final V getReferenceSourceVertex( ) {
        return vertexSet().stream().filter(v -> isRefSource(v)).findFirst().orElse(null);
    }

    /**
     * @return the reference sink vertex pulled from the graph, can be null if it doesn't exist in the graph
     */
    public final V getReferenceSinkVertex( ) {
        return vertexSet().stream().filter(v -> isRefSink(v)).findFirst().orElse(null);
    }

    /**
     * Traverse the graph and get the next reference vertex if it exists
     * @param v the current vertex, can be null
     * @return  the next reference vertex if it exists, otherwise null
     */
    public final V getNextReferenceVertex( final V v ) {
        return getNextReferenceVertex(v, false, Optional.<E>empty());
    }

    /**
     * Traverse the graph and get the next reference vertex if it exists
     * @param v the current vertex, can be null
     * @param allowNonRefPaths if true, allow sub-paths that are non-reference if there is only a single outgoing edge
     * @param blacklistedEdge optional edge to ignore in the traversal down; useful to exclude the non-reference dangling paths
     * @return the next vertex (but not necessarily on the reference path if allowNonRefPaths is true) if it exists, otherwise null
     */
    public final V getNextReferenceVertex( final V v, final boolean allowNonRefPaths, final Optional<E> blacklistedEdge ) {
        if( v == null ) { return null; }

        final Set<E> outgoingEdges = outgoingEdgesOf(v);

        if (outgoingEdges.isEmpty()){
            return null;
        }

        for( final E edgeToTest : outgoingEdges ) {
            if( edgeToTest.isRef() ) {
                return getEdgeTarget(edgeToTest);
            }
        }

        if (!allowNonRefPaths){
            return null;
        }

        //singleton or empty set
        final Set<E> blacklistedEdgeSet = blacklistedEdge.isPresent() ? Collections.singleton(blacklistedEdge.get()) : Collections.emptySet();

        // if we got here, then we aren't on a reference path
        final Optional<E> edge = outgoingEdges.stream().filter(e -> !blacklistedEdgeSet.contains(e)).findAny();
        return edge.isPresent() ? getEdgeTarget(edge.get()) : null;
    }

    /**
     * Traverse the graph and get the previous reference vertex if it exists
     * @param v the current vertex, can be null
     * @return  the previous reference vertex if it exists or null otherwise.
     */
    public final V getPrevReferenceVertex( final V v ) {
        if( v == null ) { return null; }
        return incomingEdgesOf(v).stream().map(e -> getEdgeSource(e)).filter(vrtx -> isReferenceNode(vrtx)).findFirst().orElse(null);
    }

    /**
     * Walk along the reference path in the graph and pull out the corresponding bases
     * @param fromVertex    starting vertex
     * @param toVertex      ending vertex
     * @param includeStart  should the starting vertex be included in the path
     * @param includeStop   should the ending vertex be included in the path
     * @return              byte[] array holding the reference bases, this can be null if there are no nodes between the starting and ending vertex (insertions for example)
     */
    public final byte[] getReferenceBytes( final V fromVertex, final V toVertex, final boolean includeStart, final boolean includeStop ) {
        Utils.nonNull(fromVertex, "Starting vertex in requested path cannot be null.");
        Utils.nonNull(toVertex, "From vertex in requested path cannot be null.");

        byte[] bytes = null;
        V v = fromVertex;
        if( includeStart ) {
            bytes = ArrayUtils.addAll(bytes, getAdditionalSequence(v));
        }
        v = getNextReferenceVertex(v); // advance along the reference path
        while( v != null && !v.equals(toVertex) ) {
            bytes = ArrayUtils.addAll(bytes, getAdditionalSequence(v));
            v = getNextReferenceVertex(v); // advance along the reference path
        }
        if( includeStop && v != null && v.equals(toVertex)) {
            bytes = ArrayUtils.addAll(bytes, getAdditionalSequence(v));
        }
        return bytes;
    }

    /**
     * Convenience function to add multiple vertices to the graph at once
     * @param vertices one or more vertices to add
     */
    @SafeVarargs
    @SuppressWarnings("varargs")
    public final void addVertices(final V... vertices) {
        Utils.nonNull(vertices);
        addVertices(Arrays.asList(vertices));
    }

    /**
     * Convenience function to add multiple vertices to the graph at once
     * @param vertices one or more vertices to add
     */
    public final void addVertices(final Collection<V> vertices) {
        Utils.nonNull(vertices);
        vertices.forEach(v -> addVertex(v));
    }

    /**
     * Convenience function to add multiple edges to the graph
     * @param start the first vertex to connect
     * @param remaining all additional vertices to connect
     */
    @SafeVarargs
    public final void addEdges(final V start, final V... remaining) {
        Utils.nonNull(start, "start vertex");
        if (remaining == null || remaining.length == 0){
            return;
        }
        V prev = start;
        for ( final V next : remaining ) {
            Utils.nonNull(next, "null vertex");
            addEdge(prev, next);
            prev = next;
        }
    }

    /**
     * Convenience function to add multiple edges to the graph
     * @param start the first vertex to connect
     * @param remaining all additional vertices to connect
     */
    @SafeVarargs
    public final void addEdges(final Supplier<E> template, final V start, final V... remaining) {
        Utils.nonNull(template, "template edge");
        Utils.nonNull(start, "start vertex");

        V prev = start;
        for ( final V next : remaining ) {
            Utils.nonNull(next, "null vertex");
            addEdge(prev, next, template.get());
            prev = next;
        }
    }

    /**
     * Get the set of vertices connected by outgoing edges of V
     * @param v a non-null vertex
     * @return a set of vertices connected by outgoing edges from v
     */
    public final Set<V> outgoingVerticesOf(final V v) {
        Utils.nonNull(v);
        return outgoingEdgesOf(v).stream().map(e -> getEdgeTarget(e)).collect(Collectors.toSet());
    }

    /**
     * Get the set of vertices connected to v by incoming edges
     * @param v a non-null vertex
     * @return a set of vertices {X} connected X -> v
     */
    public final Set<V> incomingVerticesOf(final V v) {
        Utils.nonNull(v);
        return incomingEdgesOf(v).stream().map(e -> getEdgeSource(e)).collect(Collectors.toSet());
    }

    /**
     * Get the set of vertices connected to v by incoming or outgoing edges
     * @param v a non-null vertex
     * @return a set of vertices {X} connected X -> v or v -> Y
     */
    public final Set<V> neighboringVerticesOf(final V v) {
        Utils.nonNull(v);
        return Sets.union(incomingVerticesOf(v), outgoingVerticesOf(v));
    }

    /**
     * Print out the graph in the dot language for visualization
     * @param destination File to write to
     */
    public final void printGraph(final File destination, final int pruneFactor) {
        try (PrintStream stream = new PrintStream(new FileOutputStream(destination))) {
            printGraph(stream, true, pruneFactor);
        } catch ( final FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(destination, e);
        }
    }

    public final void printGraph(final PrintStream graphWriter, final boolean writeHeader, final int pruneFactor) {
        if ( writeHeader ) {
            graphWriter.println("digraph assemblyGraphs {");
        }

        for( final E edge : edgeSet() ) {
            final String edgeString =  String.format("\t%s -> %s ", getEdgeSource(edge).toString(), getEdgeTarget(edge).toString());
            final String edgeLabelString;
            if (edge.getMultiplicity() > 0 && edge.getMultiplicity() < pruneFactor){
                edgeLabelString = String.format("[style=dotted,color=grey,label=\"%s\"];", edge.getDotLabel());
            } else {
                edgeLabelString = String.format("[label=\"%s\"];", edge.getDotLabel());
            }
            graphWriter.print(edgeString);
            graphWriter.print(edgeLabelString);
            if( edge.isRef() ) {
                graphWriter.println(edgeString + " [color=red];");
            }
        }

        for( final V v : vertexSet() ) {
            graphWriter.println(String.format("\t%s [label=\"%s\",shape=box]", v.toString(), new String(getAdditionalSequence(v)) + v.getAdditionalInfo()));
        }

        if ( writeHeader ) {
            graphWriter.println("}");
        }
    }

    /**
     * Remove edges that are connected before the reference source and after the reference sink
     *
     * Also removes all vertices that are orphaned by this process
     */
    public final void cleanNonRefPaths() {
        if( getReferenceSourceVertex() == null || getReferenceSinkVertex() == null ) {
            return;
        }

        // Remove non-ref edges connected before and after the reference path
        final Collection<E> edgesToCheck = new HashSet<>();
        edgesToCheck.addAll(incomingEdgesOf(getReferenceSourceVertex()));
        while( !edgesToCheck.isEmpty() ) {
            final E e = edgesToCheck.iterator().next();
            if( !e.isRef() ) {
                edgesToCheck.addAll( incomingEdgesOf(getEdgeSource(e)) );
                removeEdge(e);
            }
            edgesToCheck.remove(e);
        }

        edgesToCheck.addAll(outgoingEdgesOf(getReferenceSinkVertex()));
        while( !edgesToCheck.isEmpty() ) {
            final E e = edgesToCheck.iterator().next();
            if( !e.isRef() ) {
                edgesToCheck.addAll( outgoingEdgesOf(getEdgeTarget(e)) );
                removeEdge(e);
            }
            edgesToCheck.remove(e);
        }

        removeSingletonOrphanVertices();
    }

    /**
     * Prune all chains from this graph where all edges in the path have multiplicity < pruneFactor
     *
     * @see LowWeightChainPruner for more information
     *
     * @param pruneFactor all edges with multiplicity < this factor that aren't ref edges will be removed
     */
    public final void pruneLowWeightChains( final int pruneFactor ) {
        new LowWeightChainPruner<V,E>(pruneFactor).pruneLowWeightChains(this);
    }

    /**
     * Remove all vertices in the graph that have in and out degree of 0
     */
    public void removeSingletonOrphanVertices() {
        // Run through the graph and clean up singular orphaned nodes
        //Note: need to collect nodes to remove first because we can't directly modify the list we're iterating over
        final List<V> toRemove = vertexSet().stream().filter(v -> isSingletonOrphan(v)).collect(Collectors.toList());
        removeAllVertices(toRemove);
    }

    private boolean isSingletonOrphan(final V v) {
        Utils.nonNull(v);
        return inDegreeOf(v) == 0 && outDegreeOf(v) == 0 && !isRefSource(v);
    }

    /**
     * Remove all vertices on the graph that cannot be accessed by following any edge,
     * regardless of its direction, from the reference source vertex
     */
    public final void removeVerticesNotConnectedToRefRegardlessOfEdgeDirection() {
        final Collection<V> toRemove = new HashSet<>(vertexSet());

        final V refV = getReferenceSourceVertex();
        if ( refV != null ) {
            for ( final V v : new BaseGraphIterator<>(this, refV, true, true) ) {
                toRemove.remove(v);
            }
        }

        removeAllVertices(toRemove);
    }

    /**
     * Remove all vertices in the graph that aren't on a path from the reference source vertex to the reference sink vertex
     *
     * More aggressive reference pruning algorithm than removeVerticesNotConnectedToRefRegardlessOfEdgeDirection,
     * as it requires vertices to not only be connected by a series of directed edges but also prunes away
     * paths that do not also meet eventually with the reference sink vertex
     */
    public final void removePathsNotConnectedToRef() {
        if ( getReferenceSourceVertex() == null || getReferenceSinkVertex() == null ) {
            throw new IllegalStateException("Graph must have ref source and sink vertices");
        }

        // get the set of vertices we can reach by going forward from the ref source
        final Collection<V> onPathFromRefSource = new HashSet<>(vertexSet().size());
        for ( final V v : new BaseGraphIterator<>(this, getReferenceSourceVertex(), false, true) ) {
            onPathFromRefSource.add(v);
        }

        // get the set of vertices we can reach by going backward from the ref sink
        final Collection<V> onPathFromRefSink = new HashSet<>(vertexSet().size());
        for ( final V v : new BaseGraphIterator<>(this, getReferenceSinkVertex(), true, false) ) {
            onPathFromRefSink.add(v);
        }

        // we want to remove anything that's not in both the sink and source sets
        final Collection<V> verticesToRemove = new HashSet<>(vertexSet());
        onPathFromRefSource.retainAll(onPathFromRefSink);
        verticesToRemove.removeAll(onPathFromRefSource);
        removeAllVertices(verticesToRemove);

        // simple sanity checks that this algorithm is working.
        if ( getSinks().size() > 1 ) {
            throw new IllegalStateException("Should have eliminated all but the reference sink, but found " + getSinks());
        }

        if ( getSources().size() > 1 ) {
            throw new IllegalStateException("Should have eliminated all but the reference source, but found " + getSources());
        }
    }

    /**
     * Semi-lenient comparison of two graphs, truing true if g1 and g2 have similar structure
     *
     * By similar this means that both graphs have the same number of vertices, where each vertex can find
     * a vertex in the other graph that's seqEqual to it.  A similar constraint applies to the edges,
     * where all edges in g1 must have a corresponding edge in g2 where both source and target vertices are
     * seqEqual
     *
     * @param g1 the first graph to compare
     * @param g2 the second graph to compare
     * @param <T> the type of the nodes in those graphs
     * @return true if g1 and g2 are equals
     */
    public static <T extends BaseVertex, E extends BaseEdge> boolean graphEquals(final BaseGraph<T,E> g1, final BaseGraph<T,E> g2) {
        Utils.nonNull(g1, "g1");
        Utils.nonNull(g2, "g2");
        final Set<T> vertices1 = g1.vertexSet();
        final Set<T> vertices2 = g2.vertexSet();
        final Set<E> edges1 = g1.edgeSet();
        final Set<E> edges2 = g2.edgeSet();

        if ( vertices1.size() != vertices2.size() || edges1.size() != edges2.size() ) {
            return false;
        }

        //for every vertex in g1 there is a vertex in g2 with an equal getSequenceString
        final boolean ok= vertices1.stream().map(v1 -> v1.getSequenceString()).allMatch(v1seqString -> vertices2.stream().anyMatch(v2 -> v1seqString.equals(v2.getSequenceString())));
        if (! ok){
            return false;
        }

        //for every edge in g1 there is an equal edge in g2
        final boolean okG1 = edges1.stream().allMatch(e1 -> edges2.stream().anyMatch(e2 -> g1.seqEquals(e1, e2, g2)));
        if (! okG1){
            return false;
        }
        //for every edge in g2 there is an equal edge in g1
        final boolean okG2 = edges2.stream().allMatch(e2 -> edges1.stream().anyMatch(e1 -> g2.seqEquals(e2, e1, g1)));
        return okG2;
    }

    // For use when comparing edges across graphs!
    private boolean seqEquals( final E edge1, final E edge2, final BaseGraph<V,E> graph2 ) {
        return (this.getEdgeSource(edge1).seqEquals(graph2.getEdgeSource(edge2))) && (this.getEdgeTarget(edge1).seqEquals(graph2.getEdgeTarget(edge2)));
    }


    /**
     * Get the incoming edge of v.  Requires that there be only one such edge or throws an error
     * @param v our vertex
     * @return the single incoming edge to v, or null if none exists
     */
    public final E incomingEdgeOf(final V v) {
        Utils.nonNull(v);
        return getSingletonEdge(incomingEdgesOf(v));
    }

    /**
     * Get the outgoing edge of v.  Requires that there be only one such edge or throws an error
     * @param v our vertex
     * @return the single outgoing edge from v, or null if none exists
     */
    public final E outgoingEdgeOf(final V v) {
        Utils.nonNull(v);
        return getSingletonEdge(outgoingEdgesOf(v));
    }

    /**
     * Helper function that gets the a single edge from edges, null if edges is empty, or
     * throws an error is edges has more than 1 element
     * @param edges a set of edges
     * @return a edge
     */
    private E getSingletonEdge(final Collection<E> edges) {
        Utils.validateArg(edges.size() <= 1, "Cannot get a single incoming edge for a vertex with multiple incoming edges " + edges);
        return edges.isEmpty() ? null : edges.iterator().next();
    }

    /**
     * Add edge between source -> target if none exists, or add e to an already existing one if present
     *
     * @param source source vertex
     * @param target vertex
     * @param e edge to add
     */
    public final void addOrUpdateEdge(final V source, final V target, final E e) {
        Utils.nonNull(source, "source");
        Utils.nonNull(target, "target");
        Utils.nonNull(e, "edge");

        final E prev = getEdge(source, target);
        if ( prev != null ) {
            prev.add(e);
        } else {
            addEdge(source, target, e);
        }
    }

    @Override
    public String toString() {
        return "BaseGraph{" +
                "kmerSize=" + kmerSize +
                '}';
    }

    /**
     * Get the set of vertices within distance edges of source, regardless of edge direction
     *
     * @param source the source vertex to consider
     * @param distance the distance
     * @return a set of vertices within distance of source
     */
    private Set<V> verticesWithinDistance(final V source, final int distance) {
        if ( distance == 0 ) {
            return Collections.singleton(source);
        }

        final Set<V> found = new HashSet<>();
        found.add(source);
        for ( final V v : neighboringVerticesOf(source) ) {
            found.addAll(verticesWithinDistance(v, distance - 1));
        }

        return found;
    }

    /**
     * Get a graph containing only the vertices within distance edges of target
     * @param target a vertex in graph
     * @param distance the max distance
     * @return a non-null graph
     */
    public final BaseGraph<V,E> subsetToNeighbors(final V target, final int distance) {
        Utils.nonNull(target, "Target cannot be null");
        Utils.validateArg(containsVertex(target), () -> "Graph doesn't contain vertex " + target);
        Utils.validateArg(distance >= 0, () -> "Distance must be >= 0 but got " + distance);

        final Set<V> toKeep = verticesWithinDistance(target, distance);
        final Collection<V> toRemove = new HashSet<>(vertexSet());
        toRemove.removeAll(toKeep);

        final BaseGraph<V,E> result = clone();
        result.removeAllVertices(toRemove);

        return result;
    }

    /**
     * Get a subgraph of graph that contains only vertices within a given number of edges of the ref source vertex
     * @return a non-null subgraph of this graph
     */
    public final BaseGraph<V,E> subsetToRefSource(final int refSourceNeighborhood) {
        Utils.validateArg(refSourceNeighborhood > 0, () -> "refSourceNeighborhood needs to be positive but was " + refSourceNeighborhood);
        return subsetToNeighbors(getReferenceSourceVertex(), refSourceNeighborhood);
    }

    /**
     * Checks whether the graph contains all the vertices in a collection.
     *
     * @param vertices the vertices to check. Must not be null and must not contain a null.
     *
     * @throws IllegalArgumentException if {@code vertices} is {@code null}.
     *
     * @return {@code true} if all the vertices in the input collection are present in this graph.
     * Also if the input collection is empty. Otherwise it returns {@code false}.
     */
    public final boolean containsAllVertices(final Collection<? extends V> vertices) {
        Utils.nonNull(vertices, "the input vertices collection cannot be null");
        Utils.containsNoNull(vertices, "null vertex");
        return vertices.stream().allMatch(v -> containsVertex(v));
    }

    /**
     * Checks for the presence of directed cycles in the graph.
     *
     * @return {@code true} if the graph has cycles, {@code false} otherwise.
     */
    public final boolean hasCycles() {
        return new CycleDetector<>(this).detectCycles();
    }

    @Override
    @SuppressWarnings("unchecked")
    public BaseGraph<V,E> clone()  {
        return (BaseGraph<V,E>) super.clone();
    }

    /**
     * General iterator that can iterate over all vertices in a BaseGraph, following either
     * incoming, outgoing edge (as well as both or none) edges.  Supports traversal of graphs
     * with cycles and other crazy structures.  Will only ever visit each vertex once.  The
     * order in which the vertices are visited is undefined.
     */
    private static final class BaseGraphIterator<T extends BaseVertex, E extends BaseEdge> implements Iterator<T>, Iterable<T> {
        final Collection<T> visited = new HashSet<>();
        final Deque<T> toVisit = new LinkedList<>();
        final BaseGraph<T,E> graph;
        final boolean followIncomingEdges;
        final boolean followOutgoingEdges;

        /**
         * Create a new BaseGraphIterator starting its traversal at start
         *
         * Note that if both followIncomingEdges and followOutgoingEdges are false, we simply return the
         * start vertex
         *
         * @param graph the graph to iterator over.  Cannot be null
         * @param start the vertex to start at.  Cannot be null
         * @param followIncomingEdges should we follow incoming edges during our
         *                            traversal? (goes backward through the graph)
         * @param followOutgoingEdges should we follow outgoing edges during out traversal?
         */
        private BaseGraphIterator(final BaseGraph<T,E> graph, final T start,
                                 final boolean followIncomingEdges, final boolean followOutgoingEdges) {
            Utils.nonNull(graph, "graph cannot be null");
            Utils.nonNull(start, "start cannot be null");
            Utils.validateArg(graph.containsVertex(start), () -> "start " + start + " must be in graph but it isn't");
            this.graph = graph;
            this.followIncomingEdges = followIncomingEdges;
            this.followOutgoingEdges = followOutgoingEdges;

            toVisit.add(start);
        }

        @Override
        public Iterator<T> iterator() {
            return this;
        }

        @Override
        public boolean hasNext() {
            return ! toVisit.isEmpty();
        }

        @Override
        public T next() {
            final T v = toVisit.pop();

            if ( ! visited.contains(v) ) {
                visited.add(v);
                if ( followIncomingEdges ) {
                    toVisit.addAll(graph.incomingVerticesOf(v));
                }
                if ( followOutgoingEdges ) {
                    toVisit.addAll(graph.outgoingVerticesOf(v));
                }
            }

            return v;
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException("Doesn't implement remove");
        }
    }
}
