package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.Kmer;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqGraph;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

public class ExperimentalReadThreadingGraph extends ReadThreadingGraph {
    private Kmer referenceStartKmer;

    public ExperimentalReadThreadingGraph(int kmerSize) {
        super(kmerSize);
    }

    /**
     * Create a new ReadThreadingAssembler using kmerSize for matching
     * @param kmerSize must be >= 1
     */
    ExperimentalReadThreadingGraph(final int kmerSize, final boolean debugGraphTransformations, final byte minBaseQualityToUseInAssembly, final int numPruningSamples) {
        super(kmerSize, debugGraphTransformations, minBaseQualityToUseInAssembly, numPruningSamples);
    }

    /**
     * Find vertex and its position in seqForKmers where we should start assembling seqForKmers
     *
     * @param seqForKmers the sequence we want to thread into the graph
     * @return the position of the starting vertex in seqForKmer, or -1 if it cannot find one
     */
    @Override
    protected int findStart(final SequenceForKmers seqForKmers) {
        if ( seqForKmers.isRef ) {
            referenceStartKmer = new Kmer(seqForKmers.sequence, 0, kmerSize);
            return 0;
        }

        for ( int i = seqForKmers.start; i < seqForKmers.stop - kmerSize; i++ ) {
            final Kmer kmer1 = new Kmer(seqForKmers.sequence, i, kmerSize);
            if ( isThreadingStart(kmer1) ) {
                return i;
            }
        }

        return -1;
    }


    /**
     * Compute the smallest kmer size >= minKmerSize and <= maxKmerSize that has no non-unique kmers
     * among all sequences added to the current graph.  Will always return a result for maxKmerSize if
     * all smaller kmers had non-unique kmers.
     *
     * @param kmerSize the kmer size to check for non-unique kmers of
     * @return a non-null NonUniqueResult
     */
    @Override
    protected Set<Kmer> determineNonUniques(final int kmerSize) {
        return Collections.emptySet();
    }

    /**
     * In this version of the graph, we treat all kmers as valid starts
     * @param kmer the query kmer.
     * @return
     */
    @Override
    protected boolean isThreadingStart(final Kmer kmer) {
        return !startThreadingOnlyAtExistingVertex || uniqueKmers.containsKey(kmer);
    }


    /**
     * Workhorse routine of the assembler.  Given a sequence whose last vertex is anchored in the graph, extend
     * the graph one bp according to the bases in sequence.
     *
     * @param prevVertex a non-null vertex where sequence was last anchored in the graph
     * @param sequence the sequence we're threading through the graph
     * @param kmerStart the start of the current kmer in graph we'd like to add
     * @param count the number of observations of this kmer in graph (can be > 1 for GGA)
     * @param isRef is this the reference sequence?
     * @return a non-null vertex connecting prevVertex to in the graph based on sequence
     */
    @Override
    protected MultiDeBruijnVertex extendChainByOne(final MultiDeBruijnVertex prevVertex, final byte[] sequence, final int kmerStart, final int count, final boolean isRef) {
        final Set<MultiSampleEdge> outgoingEdges = outgoingEdgesOf(prevVertex);

        final int nextPos = kmerStart + kmerSize - 1;
        for ( final MultiSampleEdge outgoingEdge : outgoingEdges ) {
            final MultiDeBruijnVertex target = getEdgeTarget(outgoingEdge);
            if ( target.getSuffix() == sequence[nextPos] ) {
                // we've got a match in the chain, so simply increase the count of the edge by 1 and continue
                outgoingEdge.incMultiplicity(count);
                return target;
            }
        }

        // none of our outgoing edges had our unique suffix base, so we check for an opportunity to merge back in
        final Kmer kmer = new Kmer(sequence, kmerStart, kmerSize);
        final MultiDeBruijnVertex mkergeVertex = uniqueKmers.get(kmer);

        // either use our unique merge vertex, or create a new one in the chain
        final MultiDeBruijnVertex nextVertex = mkergeVertex == null ? createVertex(kmer) : mkergeVertex;
        addEdge(prevVertex, nextVertex, ((MyEdgeFactory)getEdgeFactory()).createEdge(isRef, count));
        return nextVertex;
    }


    // Generate a SeqGraph that is special
    @Override
    public SeqGraph toSequenceGraph() {
        buildGraphIfNecessary();
        return super.toSequenceGraph();
    }






    
    private class ThreadingTree {
        ThreadingNode rootNode;
        Kmer nodeBase;

        public ThreadingTree(Kmer vertex) {
            nodeBase = vertex;
        }


        public void threadRead(List<MultiSampleEdge> branchedEdges) {
            rootNode.incrementNode(branchedEdges, 0);
        }

    }

    //TODO these maps probably could just be a 4 element array keyed on the base
    private class ThreadingNode {
        private Map<MultiSampleEdge, ThreadingNode> childrenNodes;
        private MultiSampleEdge edgeInQuesiton;
        private int count = 0;

        public ThreadingNode(MultiSampleEdge edge) {
            edgeInQuesiton = edge;
            childrenNodes = new HashMap<>(2);
        }

        public void incrementNode(List<MultiSampleEdge> edges, int index) {
            count++;
            if (index < edges.size()) {
                getNode(edges.get(index)).incrementNode(edges, index + 1);
            }
        }

        private ThreadingNode getNode(MultiSampleEdge edge) {
            ThreadingNode nextNode = childrenNodes.get(edge);
            if (nextNode == null) {
                nextNode = new ThreadingNode(edge);
                childrenNodes.put(edge, nextNode);
            }
            return nextNode;
        }

        public int getCount() {
            return count;
        }
    }

}
