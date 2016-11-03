package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
 * Split a collection of middle nodes in a graph into their shared prefix and suffix values
 *
 * This code performs the following transformation.  Suppose I have a set of vertices V, such
 * that each vertex is composed of sequence such that
 *
 * Vi = prefix + seq_i + suffix
 *
 * where prefix and suffix are shared sequences across all vertices V
 *
 * This algorithm creates a new SeqGraph with the following configuration
 *
 * prefix -> has outgoing edges to all seq_i
 * suffix -> has incoming edges for all seq_i
 *
 * There are a few special cases that must be handled.  First, Vi could be simply
 * == to the prefix or the suffix.  These generate direct connections between
 * the prefix and suffix nodes, and they are handled internally by the algorithm.
 *
 * Note that for convenience, we will always create newTop and newBottom nodes, but
 * these may be empty node (i.e., they contain no sequence).  That allows them to be
 * trivially merged, if desired, when the graph is incorporated into an overall
 * graph.
 *
 * The product of this operation is a SeqGraph that contains the split.  There's a
 * function to merge reconnect this graph into the graph that contains the middle nodes
 *
 * The process guarentees a few things about the output:
 *
 * -- Preserves the paths and weights among all vertices
 *
 * It produces a graph that has some unusual properties
 *
 * -- May add nodes with no sequence (isEmpty() == true) to preserve connectivity among the graph
 * -- May introduce edges with no multiplicity to preserve paths through the graph
 *
 * The overall workflow of using this class is simple:
 *
 * find vertices V in graph that you want to split out
 * s = new SharedVertexSequenceSplitter(graph, V)
 * s.updateGraph(graph)
 *
 * to update the graph with the modifications created by this splitter
 */
public final class SharedVertexSequenceSplitter {
    private final SeqGraph outer;
    private final SeqVertex prefixV;
    private final SeqVertex suffixV;
    private final Collection<SeqVertex> toSplits;

    // updated in split routine
    private SeqGraph splitGraph = null;
    private Collection<SeqVertex> newMiddles = null;
    private List<BaseEdge> edgesToRemove = null;

    /**
     * Create a utility that will change the given graph so that the vertices in toSplitsArg (with their shared suffix and prefix
     * sequences) are extracted out.
     *
     * @param graph the graph containing the vertices in toSplitsArg
     * @param toSplitsArg a collection of vertices to split.  Must be contained within graph, and have only connections
     *                    from a single shared top and/or bottom node
     */
    public SharedVertexSequenceSplitter(final SeqGraph graph, final Collection<SeqVertex> toSplitsArg) {
        Utils.nonNull(graph, "graph cannot be null");
        Utils.nonNull(toSplitsArg, "toSplitsArg cannot be null");
        Utils.validateArg( toSplitsArg.size() > 1, () -> "Can only split at least 2 vertices but only got " + toSplitsArg);
        Utils.validateArg(graph.vertexSet().containsAll(toSplitsArg), "graph doesn't contain all of the vertices to split");

        this.outer = graph;
        this.toSplits = toSplitsArg;

        // all of the edges point to the same sink, so it's time to merge
        final Pair<SeqVertex, SeqVertex> prefixAndSuffix = commonPrefixAndSuffixOfVertices(toSplits);
        prefixV = prefixAndSuffix.getLeft();
        suffixV = prefixAndSuffix.getRight();
    }

    /**
     * Simple single-function interface to split and then update a graph
     *
     * @see #updateGraph(SeqVertex, SeqVertex) for a full description of top and bottom
     *
     * @param top the top vertex, may be null
     * @param bottom the bottom vertex, may be null
     * @return true if some useful splitting was done, false otherwise
     */
    public boolean splitAndUpdate(final SeqVertex top, final SeqVertex bottom) {
        split();
        updateGraph(top, bottom);
        return true;
    }

    /**
     * Does either the common suffix or prefix have at least minCommonSequence bases in it?
     * @param minCommonSequence a minimum length of the common sequence, must be >= 0
     * @return true if either suffix or prefix length >= minCommonSequence
     */
    public boolean meetsMinMergableSequenceForEitherPrefixOrSuffix(final int minCommonSequence) {
        return meetsMinMergableSequenceForPrefix(minCommonSequence) || meetsMinMergableSequenceForSuffix(minCommonSequence);
    }

    /**
     * Does the common prefix have at least minCommonSequence bases in it?
     * @param minCommonSequence a minimum length of the common sequence, must be >= 0
     * @return true if prefix length >= minCommonSequence
     */
    public boolean meetsMinMergableSequenceForPrefix(final int minCommonSequence) {
        return getPrefixV().length() >= minCommonSequence;
    }

    /**
     * Does the common suffix have at least minCommonSequence bases in it?
     * @param minCommonSequence a minimum length of the common sequence, must be >= 0
     * @return true if suffix length >= minCommonSequence
     */
    public boolean meetsMinMergableSequenceForSuffix(final int minCommonSequence) {
        return getSuffixV().length() >= minCommonSequence;
    }

    /**
     * Actually do the splitting up of the vertices
     *
     * Must be called before calling updateGraph
     */
    public void split() {
        splitGraph = new SeqGraph(outer.getKmerSize());
        newMiddles = new LinkedList<>();
        edgesToRemove = new LinkedList<>();

        splitGraph.addVertices(getPrefixV(), getSuffixV());

        for ( final SeqVertex mid : toSplits ) {
            final BaseEdge toMid = processEdgeToRemove(mid, outer.incomingEdgeOf(mid));
            final BaseEdge fromMid = processEdgeToRemove(mid, outer.outgoingEdgeOf(mid));

            final SeqVertex remaining = mid.withoutPrefixAndSuffix(getPrefixV().getSequence(), getSuffixV().getSequence());
            if ( remaining != null ) {
                // there's some sequence prefix + seq + suffix, so add the node and make edges
                splitGraph.addVertex(remaining);
                getNewMiddles().add(remaining);
                // update edge from top -> middle to be top -> without suffix
                splitGraph.addEdge(getPrefixV(), remaining, toMid);
                splitGraph.addEdge(remaining, getSuffixV(), fromMid);
            } else {
                // prefix + suffix completely explain this node
                splitGraph.addOrUpdateEdge(getPrefixV(), getSuffixV(), toMid.copy().add(fromMid));
            }
        }
    }

    /**
     * Update graph outer, replacing the previous middle vertices that were split out with the new
     * graph structure of the split, linking this subgraph into the graph at top and bot (the
     * vertex connecting the middle nodes and the vertex outgoing of all middle node)
     *
     * @param top an optional top node that must have outgoing edges to all split vertices.  If null, this subgraph
     *            will be added without any incoming edges
     * @param bot an optional bottom node that must have incoming edges to all split vertices.  If null, this subgraph
     *            will be added without any outgoing edges to the rest of the graph
     */
    public void updateGraph(final SeqVertex top, final SeqVertex bot) {
        Utils.validateArg(outer.vertexSet().containsAll(toSplits), "graph doesn't contain all of the original vertices to split");
        Utils.validateArg( top != null || bot != null, "Cannot update graph without at least one top or bot vertex, but both were null");
        Utils.validateArg( top == null || outer.containsVertex(top), () -> "top " + top + " not found in graph " + outer);
        Utils.validateArg( bot == null || outer.containsVertex(bot), () -> "bot " + bot + " not found in graph " + outer);
        if ( splitGraph == null ) {
            throw new IllegalStateException("Cannot call updateGraph until split() has been called");
        }

        outer.removeAllVertices(toSplits);
        outer.removeAllEdges(edgesToRemove);

        outer.addVertices(getNewMiddles());

        final boolean hasPrefixSuffixEdge = splitGraph.getEdge(getPrefixV(), getSuffixV()) != null;
        final boolean hasOnlyPrefixSuffixEdges = hasPrefixSuffixEdge && splitGraph.outDegreeOf(getPrefixV()) == 1;
        final boolean needPrefixNode = ! getPrefixV().isEmpty() || (top == null && ! hasOnlyPrefixSuffixEdges);
        final boolean needSuffixNode = ! getSuffixV().isEmpty() || (bot == null && ! hasOnlyPrefixSuffixEdges);

        // if prefix / suffix are needed, keep them
        final SeqVertex topForConnect = needPrefixNode ? getPrefixV() : top;
        final SeqVertex botForConnect = needSuffixNode ? getSuffixV() : bot;

        if ( needPrefixNode ) {
            addPrefixNodeAndEdges(top);
        }

        if ( needSuffixNode ) {
            addSuffixNodeAndEdges(bot);
        }

        if ( topForConnect != null ) {
            addEdgesFromTopNode(topForConnect, botForConnect);
        }

        if ( botForConnect != null ) {
            addEdgesToBottomNode(botForConnect);
        }
    }

    private void addEdgesToBottomNode(final SeqVertex botForConnect) {
        for ( final BaseEdge e : splitGraph.incomingEdgesOf(getSuffixV()) ) {
            outer.addEdge(splitGraph.getEdgeSource(e), botForConnect, e);
        }
    }

    private void addEdgesFromTopNode(final SeqVertex topForConnect, final SeqVertex botForConnect) {
        for ( final BaseEdge e : splitGraph.outgoingEdgesOf(getPrefixV()) ) {
            final SeqVertex target = splitGraph.getEdgeTarget(e);

            if ( target == getSuffixV()) { // going straight from prefix -> suffix
                if ( botForConnect != null ) {
                    outer.addEdge(topForConnect, botForConnect, e);
                }
            } else {
                outer.addEdge(topForConnect, target, e);
            }
        }
    }

    private void addSuffixNodeAndEdges(final SeqVertex bot) {
        outer.addVertex(getSuffixV());
        if ( bot != null ) {
            outer.addEdge(getSuffixV(), bot, BaseEdge.makeOREdge(splitGraph.incomingEdgesOf(getSuffixV()), 1));
        }
    }

    private void addPrefixNodeAndEdges(final SeqVertex top) {
        outer.addVertex(getPrefixV());
        if ( top != null ) {
            outer.addEdge(top, getPrefixV(), BaseEdge.makeOREdge(splitGraph.outgoingEdgesOf(getPrefixV()), 1));
        }
    }

    /**
     * Return the longest suffix of bases shared among all provided vertices
     *
     * For example, if the vertices have sequences AC, CC, and ATC, this would return
     * a single C.  However, for ACC and TCC this would return CC.  And for AC and TG this
     * would return null;
     *
     * @param middleVertices a non-empty set of vertices
     * @return
     */
    @VisibleForTesting
    static Pair<SeqVertex, SeqVertex> commonPrefixAndSuffixOfVertices(final Collection<SeqVertex> middleVertices) {
        final List<byte[]> kmers = new ArrayList<>(middleVertices.size());

        int min = Integer.MAX_VALUE;
        for ( final SeqVertex v : middleVertices ) {
            kmers.add(v.getSequence());
            min = Math.min(min, v.getSequence().length);
        }

        final int prefixLen = GraphUtils.commonMaximumPrefixLength(kmers);
        final int suffixLen = GraphUtils.commonMaximumSuffixLength(kmers, min - prefixLen);

        final byte[] kmer = kmers.get(0);
        final byte[] prefix = Arrays.copyOfRange(kmer, 0, prefixLen);
        final byte[] suffix = Arrays.copyOfRange(kmer, kmer.length - suffixLen, kmer.length);
        return new MutablePair<>(new SeqVertex(prefix), new SeqVertex(suffix));
    }

    /**
     * Helper function that returns an edge that we should use for splitting
     *
     * If e is null, creates a new 0 multiplicity edge, set to ref is any edges to V are ref
     * If e is not null, returns a new copy of e, and schedules e for removal
     *
     * @param e a non-null edge
     * @return a non-null edge
     */
    private BaseEdge processEdgeToRemove(final SeqVertex v, final BaseEdge e) {
        if ( e == null ) {
            // there's no edge, so we return a newly allocated one and don't schedule e for removal
            // the weight must be 0 to preserve sum through the diamond
            return new BaseEdge(outer.isReferenceNode(v), 0);
        } else {
            // schedule edge for removal, and return a freshly allocated one for our graph to use
            edgesToRemove.add(e);
            return e.copy();
        }
    }

    @VisibleForTesting
    SeqVertex getPrefixV() {
        return prefixV;
    }

    @VisibleForTesting
    SeqVertex getSuffixV() {
        return suffixV;
    }

    @VisibleForTesting
    SeqGraph getSplitGraph() {
        return splitGraph;
    }

    @VisibleForTesting
    Collection<SeqVertex> getNewMiddles() {
        return newMiddles;
    }
}