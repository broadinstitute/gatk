package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.Kmer;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqGraph;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;

import java.util.*;
import java.util.stream.Collectors;

public class ExperimentalReadThreadingGraph extends ReadThreadingGraph {
    private Kmer referenceStartKmer;
    private Map<MultiDeBruijnVertex, ThreadingTree> readThreadingJunctionTrees;

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

    protected int findStartForJuncitonThreading(final SequenceForKmers seqForKmers) {
        for ( int i = seqForKmers.start; i < seqForKmers.stop - kmerSize; i++ ) {
            final Kmer kmer1 = new Kmer(seqForKmers.sequence, i, kmerSize);
            if ( isJTStart(kmer1) ) {
                return i;
            }
        }

        return -1;
    }

    @Override
    // We don't want to remove pending sequences as they are the data we need for read threading
    protected void removePendingSequencesIfNecessary() {
        return;
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
    protected boolean isJTStart(final Kmer kmer) {
        return uniqueKmers.containsKey(kmer);
    }

    /**
     * Finds the path in the graph from this vertex to the reference sink, including this vertex
     *
     * @param start   the reference vertex to start from
     * @param direction describes which direction to move in the graph (i.e. down to the reference sink or up to the source)
     * @param blacklistedEdge edge to ignore in the traversal down; useful to exclude the non-reference dangling paths
     * @return the path (non-null, non-empty)
     */
    @Override
    protected List<MultiDeBruijnVertex> getReferencePath(final MultiDeBruijnVertex start,
                                                         final TraversalDirection direction,
                                                         final Optional<MultiSampleEdge> blacklistedEdge) {
        return direction == TraversalDirection.downwards ? referencePath : Lists.reverse(referencePath);
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



    /**
     * Extends the current vertex chain by one, making sure to build junction tree objects at each node with out-degree > 1.
     *
     *
     * @param prevVertex
     * @param sequence
     * @param kmerStart
     * @return
     */
    protected MultiDeBruijnVertex extendJunctionThreadingByOne(final MultiDeBruijnVertex prevVertex, final byte[] sequence, final int kmerStart, List<ThreadingNode> nodesToExtend) {
        final Set<MultiSampleEdge> outgoingEdges = outgoingEdgesOf(prevVertex);
        ThreadingNode nodeAtCurrentBase = null;

        if (outgoingEdges.size() != 1) {
            nodeAtCurrentBase = readThreadingJunctionTrees.computeIfAbsent(prevVertex, k -> new ThreadingTree(prevVertex)).getAndIncrementRootNode();
        }

        final int nextPos = kmerStart + kmerSize - 1;
        for (final MultiSampleEdge outgoingEdge : outgoingEdges) {
            final MultiDeBruijnVertex target = getEdgeTarget(outgoingEdge);
            if (target.getSuffix() == sequence[nextPos]) {
                // If this vertex is a junction that we must record
                if (nodeAtCurrentBase != null) {
                    //TODO this needs to be handled in a better way...
                    nodesToExtend.add(nodeAtCurrentBase);
                    List<ThreadingNode> newNodes = nodesToExtend.stream().map(n -> n.addEdge(outgoingEdge)).collect(Collectors.toList());
                    nodesToExtend.clear();
                    nodesToExtend.addAll(newNodes);
                }
                return target;
            }
        }
        return null;
    }

    //TODO this is a placeholder for when these mechanisms should be more smart, as of right now I can live with the consequnces
    //TODO the primary issue with this approach is the problem with choosing reference paths for smith waterman...
    @Override
    public void recoverDanglingHeads(final int pruneFactor, final int minDanglingBranchLength, final boolean recoverAll, final SmithWatermanAligner aligner) {
        return;
    }
    @Override
    public void recoverDanglingTails(final int pruneFactor, final int minDanglingBranchLength, final boolean recoverAll, final SmithWatermanAligner aligner) {
        return;
    }


    // Generate a SeqGraph that is special
    @Override
    public SeqGraph toSequenceGraph() {
        buildGraphIfNecessary();
        return super.toSequenceGraph();
    }

    // Generate threading trees
    public void generateJunctionTrees() {
        buildGraphIfNecessary();
        readThreadingJunctionTrees = new HashMap<>();
        pending.values().stream().flatMap(Collection::stream).forEach(this::threadSequenceForJuncitonTree);
    }

    //TODO update this mehtod to make some sense
    public void threadSequenceForJuncitonTree(SequenceForKmers seqForKmers) {
        // Maybe handle this differently, the reference junction tree should be held seperatedly from everything else.
        if (seqForKmers.isRef) {
            return;
        }

        // List we will use to keep track of sequences
        List<ThreadingNode> nodesUpdated = new ArrayList<>(3);

        // Find the first kmer in the read that exists on the graph
        final int startPos = findStartForJuncitonThreading(seqForKmers);
        if ( startPos == -1 ) {
            return;
        }

        final MultiDeBruijnVertex startingVertex = uniqueKmers.get(new Kmer(seqForKmers.sequence, startPos, kmerSize));

        // loop over all of the bases in sequence, extending the graph by one base at each point, as appropriate
        MultiDeBruijnVertex lastVertex = startingVertex;
        int kmersPastSinceLast = 0;
        for ( int i = startPos + 1; i <= seqForKmers.stop - kmerSize; i++ ) {
            MultiDeBruijnVertex vertex;
            if (kmersPastSinceLast == 0) {
                vertex = extendJunctionThreadingByOne(lastVertex, seqForKmers.sequence, i, nodesUpdated);
            } else {
                Kmer kmer = new Kmer(seqForKmers.sequence, i, kmerSize);
                vertex = uniqueKmers.get(kmer);
                if (vertex != null) {
                   nodesUpdated.addAll(attemptToResolveThreadingBetweenVertexes(startingVertex , vertex));
                }
            }
            // If for whatever reason vertex = null, then we have fallen off the corrected graph so we don't update anything
            if (vertex != null) {
                lastVertex = vertex;
            } else {
                kmersPastSinceLast++;
            }

        }
    }

    // TODO this needs to be filled out and resolved
    // TODO as an extension, this should be made to intelligently pick nodes if there is a fork (possibly an expensive step)
    private List<ThreadingNode> attemptToResolveThreadingBetweenVertexes(MultiDeBruijnVertex startingVertex , MultiDeBruijnVertex vertex) {
        return Collections.emptyList();
    }

    @VisibleForTesting
    public Map<MultiDeBruijnVertex, ThreadingTree> getReadThreadingJunctionTrees() {
        return Collections.unmodifiableMap(readThreadingJunctionTrees);
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Classes associated with junction trees
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    public class ThreadingTree {
        private ThreadingNode rootNode;
        private MultiDeBruijnVertex graphBase;

        public ThreadingTree(MultiDeBruijnVertex vertex) {
            graphBase = vertex;
            rootNode = new ThreadingNode(null);
        }

        public ThreadingNode getAndIncrementRootNode() {
            rootNode.incrementCount();
            return rootNode;
        }

        public void threadRead(List<MultiSampleEdge> branchedEdges) {
            rootNode.incrementNode(branchedEdges, 0);
        }

        public List<List<MultiSampleEdge>> enumeratePathsPresent() {
            return rootNode.getEdgesThroughNode();
        }

        // Returns the junction choices as a list of suffixes to follow
        @VisibleForTesting
        public List<String> getPathsPresentAsBaseChoiceStrings() {
            return rootNode.getEdgesThroughNode().stream()
                    .map(path -> path.stream()
                                    .map(edge -> getEdgeTarget(edge).getSuffix())
                                    .map(b -> Character.toString((char)b.byteValue()))
                                    .collect(Collectors.joining()))
                    .collect(Collectors.toList());

        }
    }

    private class ThreadingNode {
        private Map<MultiSampleEdge, ThreadingNode> childrenNodes;
        private MultiSampleEdge prevEdge = null; // This may be null if this node corresponds to the root of the graph
        private int count = 0;

        private ThreadingNode(MultiSampleEdge edge) {
            prevEdge = edge;
            childrenNodes = new HashMap<>(2);
        }

        /**
         * Adds an edge to the tree if this node has not been traversed before, or increments the corresponding node
         * if the edge already exists. Then returns the corresponding next node
         *
         * @param edge edge to add to this current node
         * @return the node corresponding to the one that was added to this graph
         */
        public ThreadingNode addEdge(MultiSampleEdge edge) {
            ThreadingNode nextNode = childrenNodes.computeIfAbsent(edge, k -> new ThreadingNode(edge));
            nextNode.incrementCount();
            return nextNode;
        }

        void incrementCount() {
            count++;
        }

        private void incrementNode(List<MultiSampleEdge> edges, int index) {
            incrementCount();
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

        @VisibleForTesting
        public List<List<MultiSampleEdge>> getEdgesThroughNode() {
            List<List<MultiSampleEdge>> paths = new ArrayList<>();
            if (childrenNodes.isEmpty()) {
                paths.add(prevEdge == null ? new ArrayList<>() : Lists.newArrayList(prevEdge));
            }

            for (ThreadingNode childNode : childrenNodes.values()) {
                for ( List<MultiSampleEdge> subPath : childNode.getEdgesThroughNode()) {
                    List<MultiSampleEdge> pathRoot = (prevEdge == null ? new ArrayList<>() : Lists.newArrayList(prevEdge));
                    pathRoot.addAll(subPath);
                    paths.add(pathRoot);
                }
            }

            return paths;
        }

        public int getCount() {
            return count;
        }
    }

}
