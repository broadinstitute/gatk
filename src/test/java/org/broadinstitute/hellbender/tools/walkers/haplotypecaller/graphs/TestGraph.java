package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.MultiDeBruijnVertex;
import org.broadinstitute.hellbender.utils.Utils;
import org.jgrapht.EdgeFactory;

import java.util.HashMap;
import java.util.Map;

/**
 * A Test kmer graph
 */
public final class TestGraph extends BaseGraph<MultiDeBruijnVertex, BaseEdge> {
    private static final long serialVersionUID = 1L;

    @Override
    public TestGraph clone() {
        return (TestGraph) super.clone();
    }

    /**
     * Edge factory that creates non-reference multiplicity 1 edges
     */
    private static class MyEdgeFactory implements EdgeFactory<MultiDeBruijnVertex, BaseEdge> {
        @Override
        public BaseEdge createEdge(MultiDeBruijnVertex sourceVertex, MultiDeBruijnVertex targetVertex) {
            return new BaseEdge(false, 1);
        }
    }

    /**
     * Create an empty TestGraph with default kmer size
     */
    public TestGraph() {
        this(11);
    }

    /**
     * Create an empty TestGraph with kmer size
     * @param kmerSize kmer size, must be >= 1
     */
    public TestGraph(int kmerSize) {
        super(kmerSize, new MyEdgeFactory());
    }


    /**
     * Add edge to assembly graph connecting the two kmers
     * @param kmer1 the source kmer for the edge
     * @param kmer2 the target kmer for the edge
     * @param isRef true if the added edge is a reference edge
     */
    public void addKmersToGraph( final byte[] kmer1, final byte[] kmer2, final boolean isRef, final int multiplicity ) {
        Utils.nonNull( kmer1, "Attempting to add a null kmer to the graph.");
        Utils.nonNull(kmer2,"Attempting to add a null kmer to the graph.");
        Utils.validateArg( kmer1.length == kmer2.length, "Attempting to add a kmers to the graph with different lengths.");

        final MultiDeBruijnVertex v1 = new MultiDeBruijnVertex( kmer1 , true);
        final MultiDeBruijnVertex v2 = new MultiDeBruijnVertex( kmer2 , true);
        final BaseEdge toAdd = new BaseEdge(isRef, multiplicity);

        addVertices(v1, v2);
        addOrUpdateEdge(v1, v2, toAdd);
    }

    /**
     * Convert this kmer graph to a simple sequence graph.
     *
     * Each kmer suffix shows up as a distinct SeqVertex, attached in the same structure as in the kmer
     * graph.  Nodes that are sources are mapped to SeqVertex nodes that contain all of their sequence
     *
     * @return a newly allocated SequenceGraph
     */
    @Override
    public SeqGraph toSequenceGraph() {
        final SeqGraph seqGraph = new SeqGraph(getKmerSize());
        final Map<MultiDeBruijnVertex, SeqVertex> vertexMap = new HashMap<>();

        // create all of the equivalent seq graph vertices
        for ( final MultiDeBruijnVertex dv : vertexSet() ) {
            final SeqVertex sv = new SeqVertex(dv.getAdditionalSequence(isSource(dv)));
            vertexMap.put(dv, sv);
            seqGraph.addVertex(sv);
        }

        // walk through the nodes and connect them to their equivalent seq vertices
        for( final BaseEdge e : edgeSet() ) {
            final SeqVertex seqOutV = vertexMap.get(getEdgeTarget(e));
            final SeqVertex seqInV = vertexMap.get(getEdgeSource(e));
            seqGraph.addEdge(seqInV, seqOutV, e);
        }

        return seqGraph;
    }
}
