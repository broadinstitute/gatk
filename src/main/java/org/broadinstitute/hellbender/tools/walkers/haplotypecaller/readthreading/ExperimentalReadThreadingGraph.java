package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import org.apache.commons.collections.ListUtils;
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
    private static final MultiDeBruijnVertex SYMBOLIC_END_VETEX = new MultiDeBruijnVertex(new byte[]{'_'});
    private MultiSampleEdge SYMBOLIC_END_EDGE;

    private final OneShotLogger oneShotDanglingTailLogger = new OneShotLogger(this.getClass());

    private Map<MultiDeBruijnVertex, ThreadingTree> readThreadingJunctionTrees = new HashMap<>();

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
        final int countToUse = isRef ? 0 : count; //NOTE we do not count reference multiplicity in the graph now

        final int nextPos = kmerStart + kmerSize - 1;
        for ( final MultiSampleEdge outgoingEdge : outgoingEdges ) {
            final MultiDeBruijnVertex target = getEdgeTarget(outgoingEdge);
            if ( target.getSuffix() == sequence[nextPos] ) {
                // we've got a match in the chain, so simply increase the count of the edge by 1 and continue
                outgoingEdge.incMultiplicity(countToUse);
                return target;
            }
        }

        // none of our outgoing edges had our unique suffix base, so we check for an opportunity to merge back in
        final Kmer kmer = new Kmer(sequence, kmerStart, kmerSize);
        final MultiDeBruijnVertex mkergeVertex = uniqueKmers.get(kmer);

        // either use our unique merge vertex, or create a new one in the chain
        final MultiDeBruijnVertex nextVertex = mkergeVertex == null ? createVertex(kmer) : mkergeVertex;
        addEdge(prevVertex, nextVertex, ((MyEdgeFactory)getEdgeFactory()).createEdge(isRef, countToUse));
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
                // Only want to create a new tree if we walk into a node with
                // NOTE: we do this deciding on our path
                if (vertexWarrantsJunctionTree(prevVertex)) {
                    nodesToExtend.add(readThreadingJunctionTrees.computeIfAbsent(prevVertex, k -> new ThreadingTree(prevVertex)).getAndIncrementRootNode());
                }

                // If this node has an out-degree > 1, add the edge we took to existing trees
                if (outgoingEdges.size() != 1) {
                    // TODO, make an object to encapsulate this operation better
                    addEdgeToJunctionTreeNodes(nodesToExtend, outgoingEdge);
                }

                return target;
            }
        }
        return null;
    }

    // Helper method that adds a single edge to all of the nodes in nodesToExtend.
    private static void addEdgeToJunctionTreeNodes(List<ThreadingNode> nodesToExtend, MultiSampleEdge outgoingEdge) {
        List<ThreadingNode> newNodes = nodesToExtend.stream().map(n -> n.addEdge(outgoingEdge)).collect(Collectors.toList());
        nodesToExtend.clear();
        nodesToExtend.addAll(newNodes);
    }

    // Helper method used to determine whether a vertex meets the criteria to hold a junction tree
    // The current criteria, if any outgoing edge from a particular vertex leads to a vertex with inDegree > 1, then it warrants a tree. Or if it is the reference start vertex.
    // NOTE: this check is necessary to handle the edge cases that may arise when a vertex has multiple exit paths but happens to lead to a vetex that needs a junction tree
    private boolean vertexWarrantsJunctionTree(final MultiDeBruijnVertex vertex) {
        // The reference source vertex warrants a junction tree
        if (getReferenceSourceVertex() == vertex) {
            return true;
        }

        for (MultiSampleEdge edge : outgoingEdgesOf(vertex)) {
            if (inDegreeOf(getEdgeTarget(edge)) > 1) {
                return true;
            }
        }
        return false;
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
        final List<MultiDeBruijnVertex> refPath = getReferencePathForwardFromKmer(altPath.get(0), Optional.ofNullable(getHeaviestIncomingEdge(altPath.get(1))));

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
     * //TODO both of these are a hack to emulate the current behavior when targetkmer isn't actually a reference kmer. The old behaior
     * //TODO was to return a singleton list of targetKmer. The real solution is to walk back target kmer to the actual ref base
     *
     * @param targetKmer vertex corresponding to the root
     * @return
     */
    private List<MultiDeBruijnVertex> getReferencePathForwardFromKmer(final MultiDeBruijnVertex targetKmer,
                                                                      final Optional<MultiSampleEdge> blacklistedEdge) {
        List<MultiDeBruijnVertex> extraSequence = new ArrayList<>(2);
        MultiDeBruijnVertex vert = targetKmer;
        int finalIndex = referencePath.lastIndexOf(vert);

        while (finalIndex == -1 &&  vert != null) {
            // If the current verex is not a reference vertex but exists, add it to our list of extra bases to append to the reference
            extraSequence.add(vert);

            final Set<MultiSampleEdge> outgoingEdges = outgoingEdgesOf(vert);

            // singleton or empty set
            final Set<MultiSampleEdge> blacklistedEdgeSet = blacklistedEdge.isPresent() ? Collections.singleton(blacklistedEdge.get()) : Collections.emptySet();

            // walk forward while the path is unambiguous
            final List<MultiSampleEdge> edges = outgoingEdges.stream().filter(e -> !blacklistedEdgeSet.contains(e)).limit(2).collect(Collectors.toList());

            vert = edges.size() == 1 ? getEdgeTarget(edges.get(0)) : null;
            finalIndex = vert == null ? -1 : referencePath.lastIndexOf(vert);
        }

        // if we found extra sequence append it to the front of the
        extraSequence.addAll(finalIndex != -1 ? referencePath.subList(finalIndex, referencePath.size()) : Collections.emptyList());
        return extraSequence;
    }

    // TODO this behavior is frankly silly and needs to be fixed, there is no way upwards paths should be dangingling head recovered differently
    private List<MultiDeBruijnVertex> getReferencePathBackwardsForKmer(final MultiDeBruijnVertex targetKmer) {
        int firstIndex = referencePath.indexOf(targetKmer);
        if (firstIndex == -1) return Collections.singletonList(targetKmer);
        return Lists.reverse(referencePath.subList(0, firstIndex + 1));
    }


    /**
     * Filters empty or uninformative junction trees from the graph.
     *
     * @return
     */
    public void pruneJunctionTrees(final int pruneFactor) {
        if (pruneFactor > 0) {
            throw new UnsupportedOperationException("Currently pruning based JT evidence is not supported.");
        }
        readThreadingJunctionTrees = Maps.filterValues( Collections.unmodifiableMap(readThreadingJunctionTrees), ThreadingTree::isEmptyTree);
    }

    // Generate a SeqGraph that is special
    @Override
    public SeqGraph toSequenceGraph() {
        throw new UnsupportedOperationException("Cannot construct a sequence graph using ExperimentalReadThreadingGraph");
    }

    // Generate threading trees
    public void generateJunctionTrees() {
        buildGraphIfNecessary();
        // Adding handle vertex to support symbolic end alleles
        addVertex(SYMBOLIC_END_VETEX);
        SYMBOLIC_END_EDGE = addEdge(getReferenceSinkVertex(), SYMBOLIC_END_VETEX);

        readThreadingJunctionTrees = new HashMap<>();
        pending.values().stream().flatMap(Collection::stream).forEach(this::threadSequenceForJuncitonTree);
    }

    /**
     * Takes a contiguous sequence of kmers and threads it through the graph, generating and subsequently adding to
     * juncition trees at each relevant node.
     *
     * @param seqForKmers
     */
    public void threadSequenceForJuncitonTree(final SequenceForKmers seqForKmers) {
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
            final MultiDeBruijnVertex vertex;
            if (kmersPastSinceLast == 0) {
                vertex = extendJunctionThreadingByOne(lastVertex, seqForKmers.sequence, i, nodesUpdated);
            } else {
                Kmer kmer = new Kmer(seqForKmers.sequence, i, kmerSize);
                vertex = uniqueKmers.get(kmer);
                // TODO this might cause problems
                if (vertex != null) {
                   nodesUpdated.addAll(attemptToResolveThreadingBetweenVertexes(lastVertex , vertex));
                }
            }
            // If for whatever reason vertex = null, then we have fallen off the corrected graph so we don't update anything
            if (vertex != null) {
                lastVertex = vertex;
                kmersPastSinceLast = 0;
            } else {
                kmersPastSinceLast++;
            }
        }

        // As a final step, if the last vetex happens to be the ref-stop vertex then we want to append a symbolic node to the junciton trees
        if (lastVertex == getReferenceSinkVertex()) {
            addEdgeToJunctionTreeNodes(nodesUpdated, SYMBOLIC_END_EDGE);
        }
    }

    // Extendable method intended to allow for adding extra material to the graph
    public List<String> getExtraGraphFileLines() {
        List<String> output = new ArrayList<>();
        for( Map.Entry<MultiDeBruijnVertex, ThreadingTree> entry : readThreadingJunctionTrees.entrySet()) {
            // adding the root node to the graph
            output.add(String.format("\t%s -> %s ", entry.getKey().toString(), entry.getValue().rootNode.getDotName()) +
            String.format("[color=blue];"));
            output.add(String.format("\t%s [shape=point];", entry.getValue().rootNode.getDotName()));

            output.addAll(edgesForNodeRecursive(entry.getValue().rootNode));
        }
        return output;
    }

    // Recursive search through a threading tree for nodes
    private List<String> edgesForNodeRecursive(ThreadingNode node) {
        List<String> output = new ArrayList<>();

        for ( Map.Entry<MultiSampleEdge, ThreadingNode> childrenNode : node.childrenNodes.entrySet() ) {
            output.add(String.format("\t%s -> %s ", node.getDotName(), childrenNode.getValue().getDotName() + String.format("[color=blue,label=\"%d\"];",childrenNode.getValue().count)));
            output.add(String.format("\t%s [label=\"%s\",shape=plaintext]", childrenNode.getValue().getDotName(),
                    new String(getEdgeTarget(childrenNode.getKey()).getAdditionalSequence(false))));
            output.addAll(edgesForNodeRecursive(childrenNode.getValue()));
        }
        return output;
    }

    public ThreadingTree getJunctionTreeForNode(MultiDeBruijnVertex vertex) {
        return readThreadingJunctionTrees.get(vertex);
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

        private ThreadingTree(MultiDeBruijnVertex vertex) {
            graphBase = vertex;
            rootNode = new ThreadingNode(null);
        }

        /**
         * Returns the root node for this edge (ensuring that it has had its count intcremented to keep track of the total evidence spanning this tree.
         *
         * @return
         */
        private ThreadingNode getAndIncrementRootNode() {
            rootNode.incrementCount();
            return rootNode;
        }

        // Determines if this tree actually has any data associated with it
        private boolean isEmptyTree() {
            return !rootNode.childrenNodes.isEmpty();
        }

        // Returns a list of all the paths (illustrated as sequential edges) observed from this point throug hthe graph
        private List<List<MultiSampleEdge>> enumeratePathsPresent() {
            return rootNode.getEdgesThroughNode();
        }

        // Returns the junction choices as a list of suffixes to follow
        @VisibleForTesting
        List<String> getPathsPresentAsBaseChoiceStrings() {
            return rootNode.getEdgesThroughNode().stream()
                    .map(path -> path.stream()
                                    .map(edge -> getEdgeTarget(edge).getSuffix())
                                    .map(b -> Character.toString((char)b.byteValue()))
                                    .collect(Collectors.joining()))
                    .collect(Collectors.toList());

        }

        // getter for the root node, TODO to possibly be replaced when the graph gets hidden from prying eyes.
        public ThreadingNode getRootNode() {
            return rootNode;
        }
    }

    // Linked node object for storing tree topography
    public class ThreadingNode {
        private Map<MultiSampleEdge, ThreadingNode> childrenNodes;
        private MultiSampleEdge prevEdge; // This may be null if this node corresponds to the root of the graph
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
        private ThreadingNode addEdge(MultiSampleEdge edge) {
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
        private List<List<MultiSampleEdge>> getEdgesThroughNode() {
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

        // Returns a unique name based on the memory id that conforms to the restrictions placed on .dot file nodes
        public String getDotName() {
            return "TreadingNode_" + Integer.toHexString(hashCode());
        }

        // Getter for external tools to access a node
        public Map<MultiSampleEdge, ThreadingNode> getChildrenNodes() {
            return Collections.unmodifiableMap(childrenNodes);
        }

        // Return the count of total evidence supporting this node in the tree
        public int getCount() {
            return count;
        }

        // Checks if this node is the symbolic end by determining if its previous edge ends on the symbolic edge
        public boolean isSymbolicEnd() {
            return prevEdge == SYMBOLIC_END_EDGE;
        }
    }

}
