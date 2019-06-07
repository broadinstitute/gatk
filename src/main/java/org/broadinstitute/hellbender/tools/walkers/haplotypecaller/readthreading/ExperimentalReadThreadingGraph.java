package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.Kmer;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqGraph;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignment;

import java.lang.annotation.Target;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Experimental version of the ReadThreadingGraph which does not support seqGraph construction but does support the ability
 * to
 */
public class ExperimentalReadThreadingGraph extends ReadThreadingGraphInterface {
    private final OneShotLogger oneShotDanglingTailLogger = new OneShotLogger(this.getClass());

    private Map<MultiDeBruijnVertex, ThreadingTree> readThreadingJunctionTrees;

    public ExperimentalReadThreadingGraph(int kmerSize) {
        this(kmerSize, false, (byte)6, 1);
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
    protected void resetToInitialState() {
        // We don't clear pending here in order to later build the read threading
        nonUniqueKmers = null;
        uniqueKmers.clear();
        refSource = null;
        alreadyBuilt = false;
    }

    @Override
    // We don't want to remove pending sequences as they are the data we need for read threading
    protected void removePendingSequencesIfNecessary() {
        return;
    }

    @Override
    //TODO come up with some huristic for when we think one of these graphs is "too difficult" to call properly
    public boolean isLowComplexity() {
        return false;
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

    @Override
    protected boolean baseIsUsableForAssembly(byte base, byte qual) {
        return base != BaseUtils.Base.N.base && qual >= minBaseQualityToUseInAssembly;
    }

    /**
     * Extends the current vertex chain by one, making sure to build junction tree objects at each node with out-degree > 1.
     *
     * If a node has an in-degree > 1, then it has a junction tree added to the previous node.
     *
     * @param prevVertex Current vertex to extend off of
     * @param sequence Sequence of kmers to extend
     * @param kmerStart index of where the current kmer starts
     * @return The next vertex in the sequence. Null if the corresponding edge doesn't exist.
     */
    protected MultiDeBruijnVertex extendJunctionThreadingByOne(final MultiDeBruijnVertex prevVertex, final byte[] sequence, final int kmerStart, List<ThreadingNode> nodesToExtend) {
        final Set<MultiSampleEdge> outgoingEdges = outgoingEdgesOf(prevVertex);

        final int nextPos = kmerStart + kmerSize - 1;
        for (final MultiSampleEdge outgoingEdge : outgoingEdges) {
            final MultiDeBruijnVertex target = getEdgeTarget(outgoingEdge);
            if (target.getSuffix() == sequence[nextPos]) {
                // If this node has an out-degree > 1, add the edge we took to existing trees
                if (outgoingEdges.size() != 1) {
                    // TODO, make an object to encapsulate this operation better
                    List<ThreadingNode> newNodes = nodesToExtend.stream().map(n -> n.addEdge(outgoingEdge)).collect(Collectors.toList());
                    nodesToExtend.clear();
                    nodesToExtend.addAll(newNodes);
                }

                // Only want to create a new tree if we walk into a node with
                // NOTE: we do this after extending the previous node
                if (incomingEdgesOf(target).size() > 1) {
                    nodesToExtend.add(readThreadingJunctionTrees.computeIfAbsent(prevVertex, k -> new ThreadingTree(prevVertex)).getAndIncrementRootNode());
                }

                return target;
            }
        }
        return null;
    }

    //TODO this is a placeholder for when these mechanisms should be more smart, as of right now I can live with the consequnces
    //TODO the primary issue with this approach is the problem with choosing reference paths for smith waterman...

    /**
     * Generates the CIGAR string from the Smith-Waterman alignment of the dangling path (where the
     * provided vertex is the source) and the reference path.
     *
     * @param aligner
     * @param vertex   the source of the dangling head
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @param recoverAll recover even branches with forks
     * @return a SmithWaterman object which can be null if no proper alignment could be generated
     */
    @VisibleForTesting
    @Override
    DanglingChainMergeHelper generateCigarAgainstUpwardsReferencePath(final MultiDeBruijnVertex vertex, final int pruneFactor, final int minDanglingBranchLength, final boolean recoverAll, SmithWatermanAligner aligner) {

        // find the highest common descendant path between vertex and the reference source if available
        final List<MultiDeBruijnVertex> altPath = findPathDownwardsToHighestCommonDescendantOfReference(vertex, pruneFactor, !recoverAll);
        if ( altPath == null || isRefSink(altPath.get(0)) || altPath.size() < minDanglingBranchLength + 1 ) // add 1 to include the LCA
        {
            return null;
        }

        // now get the reference path from the LCA
        final List<MultiDeBruijnVertex> refPath = getReferencePathBackwardsForKmer(altPath.get(0));

        // create the Smith-Waterman strings to use
        final byte[] refBases = getBasesForPath(refPath, true);
        final byte[] altBases = getBasesForPath(altPath, true);

        // run Smith-Waterman to determine the best alignment (and remove trailing deletions since they aren't interesting)
        final SmithWatermanAlignment alignment = aligner.align(refBases, altBases, SmithWatermanAligner.STANDARD_NGS, SWOverhangStrategy.LEADING_INDEL);
        return new DanglingChainMergeHelper(altPath, refPath, altBases, refBases, AlignmentUtils.removeTrailingDeletions(alignment.getCigar()));
    }


    /**
     * NOTE: we override the old behavior here to allow for looping/nonunique references to be handled properly
     *
     * @return the reference source vertex pulled from the graph, can be null if it doesn't exist in the graph
     */
    @Override
    public final MultiDeBruijnVertex getReferenceSourceVertex( ) {
        return referencePath != null && !referencePath.isEmpty()? referencePath.get(0) : null;
    }

    /**
     * NOTE: we override the old behavior here to allow for looping/nonunique references to be handled properly
     *
     * @return the reference sink vertex pulled from the graph, can be null if it doesn't exist in the graph
     */
    @Override
    public final MultiDeBruijnVertex getReferenceSinkVertex( ) {
        return referencePath != null && !referencePath.isEmpty()? referencePath.get(referencePath.size() - 1) : null;
    }

    /**
     * Generates the CIGAR string from the Smith-Waterman alignment of the dangling path (where the
     * provided vertex is the sink) and the reference path.
     *
     * @param aligner
     * @param vertex   the sink of the dangling chain
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @param recoverAll recover even branches with forks
     * @return a SmithWaterman object which can be null if no proper alignment could be generated
     */
    @Override
    DanglingChainMergeHelper generateCigarAgainstDownwardsReferencePath(final MultiDeBruijnVertex vertex, final int pruneFactor, final int minDanglingBranchLength, final boolean recoverAll, SmithWatermanAligner aligner) {
        final int minTailPathLength = Math.max(1, minDanglingBranchLength); // while heads can be 0, tails absolutely cannot

        // find the lowest common ancestor path between this vertex and the diverging master path if available
        final List<MultiDeBruijnVertex> altPath = findPathUpwardsToLowestCommonAncestor(vertex, pruneFactor, !recoverAll);
        if ( altPath == null || isRefSource(altPath.get(0)) || altPath.size() < minTailPathLength + 1 ) // add 1 to include the LCA
        {
            return null;
        }

        // now get the reference path from the LCA
        final List<MultiDeBruijnVertex> refPath = getReferencePathForwardFromKmer(altPath.get(0));

        // create the Smith-Waterman strings to use
        final byte[] refBases = getBasesForPath(refPath, false);
        final byte[] altBases = getBasesForPath(altPath, false);

        // run Smith-Waterman to determine the best alignment (and remove trailing deletions since they aren't interesting)
        final SmithWatermanAlignment alignment = aligner.align(refBases, altBases, SmithWatermanAligner.STANDARD_NGS, SWOverhangStrategy.LEADING_INDEL);
        return new DanglingChainMergeHelper(altPath, refPath, altBases, refBases, AlignmentUtils.removeTrailingDeletions(alignment.getCigar()));
    }

    /**
     * This is a heruistic for deciding what the proper reference path is to align the dangling ends to. The last possible
     * reference path emanating from the kmer (in the case of the root kmer for a dangling end being a repeated kmer
     * from the reference) is chose (or the first in the case for dangling heads). This is based on the assumption that
     * dangling tails worthy of recovering are often a result of the assembly window and thus we choose the last possible
     * kmer as the option.
     *
     *
     * @param targetKmer
     * @return
     */
    private List<MultiDeBruijnVertex> getReferencePathForwardFromKmer(final MultiDeBruijnVertex targetKmer) {
        int finalIndex = referencePath.lastIndexOf(targetKmer);
        return referencePath.subList(finalIndex, referencePath.size());
    }

    private List<MultiDeBruijnVertex> getReferencePathBackwardsForKmer(final MultiDeBruijnVertex targetKmer) {
        int firstIndex = referencePath.indexOf(targetKmer);
        return Lists.reverse(referencePath.subList(0, firstIndex + 1));
    }


    // Generate a SeqGraph that is special
    @Override
    public SeqGraph toSequenceGraph() {
        throw new UnsupportedOperationException("Cannot construct a sequence graph using ExperimentalReadThreadingGraph");
    }

    // Generate threading trees
    public void generateJunctionTrees() {
        buildGraphIfNecessary();
        readThreadingJunctionTrees = new HashMap<>();
        pending.values().stream().flatMap(Collection::stream).forEach(this::threadSequenceForJuncitonTree);
    }

    /**
     * Takes a contiguous sequence of kmers and threads it through the graph, generating and subsequently adding to
     * juncition trees at each relevant node.
     *
     * @param seqForKmers
     */
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
                   nodesUpdated.addAll(attemptToResolveThreadingBetweenVertexes(lastVertex , vertex));
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
    // Test method for returning all existing junction trees.
    public Map<MultiDeBruijnVertex, ThreadingTree> getReadThreadingJunctionTrees(boolean pruned) {
        return pruned ? Maps.filterValues( Collections.unmodifiableMap(readThreadingJunctionTrees), n -> n.isEmptyTree())
                : Collections.unmodifiableMap(readThreadingJunctionTrees);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Classes associated with junction trees
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Class for storing a junction tree. A Juntion tree is a tree whose nodes correspond to edges in the graph where a
     * traversal decision was made (i.e. the path forks and one edge was taken). Each node keeps track of a list of
     * its children nodes, each of which stores a count of the evidence seen supporting the given edge.
     */
    public class ThreadingTree {
        private ThreadingNode rootNode;
        private MultiDeBruijnVertex graphBase;

        public ThreadingTree(MultiDeBruijnVertex vertex) {
            graphBase = vertex;
            rootNode = new ThreadingNode(null);
        }

        /**
         * Returns the root node for this edge (ensuring that it has had its count intcremented to keep track of the total evidence spanning this tree.
         *
         * @return
         */
        public ThreadingNode getAndIncrementRootNode() {
            rootNode.incrementCount();
            return rootNode;
        }

        // Determines if this tree actually has any data associated with it
        public boolean isEmptyTree() {
            return !rootNode.childrenNodes.isEmpty();
        }

        // Returns a list of all the paths (illustrated as sequential edges) observed from this point throug hthe graph
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

    // Linked node object for storing tree topography
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

        /**
         * Helper method that returns a list of all the possible branching paths through the tree originiating from this node.
         *
         * NOTE: This method is intended for debugging and testing and thus may not be efficient for returning the list
         *
         * @return A list of potential paths originating from this node as observed
         */
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

        // Return the count of total evidence supporting this node in the tree
        public int getCount() {
            return count;
        }
    }

}
