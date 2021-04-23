package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.Kmer;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.*;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Experimental version of the ReadThreadingGraph with support for threading reads to generate JunctionTrees for resolving
 * connectivity information at longer ranges.
 *
 * Note that many of the non-DeBruijn graph alterations that are made to ReadThreadingGraph are not made here:
 * - Non-Unique kmers are not duplicated by this graph
 * - Kmers are not Zipped together to form a SeqGraph
 * - The reference path is stored in its entirety rather than being calculated on the fly
 *
 * For ease of debugging, this graph supports the method {@link #printSimplifiedGraph(File, int)}} which generates a SequenceGraph and
 * adds the junction trees to the output .dot file.
 */
public class JunctionTreeLinkedDeBruijnGraph extends AbstractReadThreadingGraph {
    private static final long serialVersionUID = 1l;
    private static final MultiDeBruijnVertex SYMBOLIC_END_VETEX = new MultiDeBruijnVertex(new byte[]{'_'});
    private MultiSampleEdge SYMBOLIC_END_EDGE;
    
    private Map<MultiDeBruijnVertex, ThreadingTree> readThreadingJunctionTrees = new HashMap<>();

    // TODO should this be constructed here or elsewhere
    private final Set<Kmer> kmers = new HashSet<>();

    @VisibleForTesting
    public JunctionTreeLinkedDeBruijnGraph(int kmerSize) {
        this(kmerSize, false, (byte)6, 1, -1);
    }

    /**
     * Create a new ReadThreadingAssembler using kmerSize for matching
     * @param kmerSize must be >= 1
     */
    JunctionTreeLinkedDeBruijnGraph(final int kmerSize, final boolean debugGraphTransformations, final byte minBaseQualityToUseInAssembly, final int numPruningSamples, final int numDanglingMatchingPrefixBases) {
        super(kmerSize, debugGraphTransformations, minBaseQualityToUseInAssembly, numPruningSamples, numDanglingMatchingPrefixBases);
    }

    /**
     * Find vertex and its position in seqForKmers where we should start assembling seqForKmers
     *
     * @param seqForKmers the sequence we want to thread into the graph
     * @return the position of the starting vertex in seqForKmer, or -1 if it cannot find one
     */

    protected int findStartForJunctionThreading(final SequenceForKmers seqForKmers) {
        for ( int i = seqForKmers.start; i < seqForKmers.stop - kmerSize; i++ ) {
            final Kmer kmer1 = new Kmer(seqForKmers.sequence, i, kmerSize);
            if ( kmerToVertexMap.containsKey(kmer1) ) {
                return i;
            }
        }

        return -1;
    }

    @Override
    // We don't need to track non-uniques here so this is a no-op
    protected void preprocessReads() {
        return;
    }

    @Override
    // We don't want to remove pending sequences as they are the data we need for read threading
    protected boolean shouldRemoveReadsAfterGraphConstruction() {
        return false;
    }

    @Override
    //TODO come up with some heuristic for when we think one of these graphs is "too difficult" to call properly
    public boolean isLowQualityGraph() {
        return false;
    }

    /**
     * Checks whether a kmer can be the threading start based on the current threading start location policy.
     *
     * @see #setThreadingStartOnlyAtExistingVertex(boolean)
     *
     * @param kmer the query kmer.
     * @return {@code true} if we can start thread the sequence at this kmer, {@code false} otherwise.
     */
    protected boolean isThreadingStart(final Kmer kmer, final boolean startThreadingOnlyAtExistingVertex) {
        Utils.nonNull(kmer);
        return !startThreadingOnlyAtExistingVertex || kmers.contains(kmer);
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
        if ( direction == TraversalDirection.downwards ) {
            return getReferencePathForwardFromKmer(start, blacklistedEdge);
        } else {
            return getReferencePathBackwardsForKmer(start);
        }
    }

    // Since there are no non-unique kmers to worry about we just add it to our map
    @Override
    protected void trackKmer(Kmer kmer, MultiDeBruijnVertex newVertex) {
        kmerToVertexMap.putIfAbsent(kmer, newVertex);
    }

    @VisibleForTesting
    List<MultiDeBruijnVertex> getReferencePath(final TraversalDirection direction) {
        return Collections.unmodifiableList(direction == TraversalDirection.downwards ? referencePath : Lists.reverse(referencePath));
    }

    /**
     * Extends the current vertex chain by one, making sure to build junction tree objects at each node with out-degree > 1.
     *
     * If a node has an in-degree > 1, then it has a junction tree added to the previous node.
     *
     * @param prevVertex Current vertex to extend off of
     * @param sequence Sequence of kmers to extend
     * @param kmerStart index of where the current kmer starts
     * @param alterTrees flag to control if we actually put anything into junction trees as a result of this path
     * @return The next vertex in the sequence. Null if the corresponding edge doesn't exist.
     */
    protected MultiDeBruijnVertex extendJunctionThreadingByOne(final MultiDeBruijnVertex prevVertex, final byte[] sequence, final int kmerStart, final JunctionTreeThreadingHelper nodesHelper, final boolean alterTrees) {
        final Set<MultiSampleEdge> outgoingEdges = outgoingEdgesOf(prevVertex);

        final int nextPos = kmerStart + kmerSize - 1;
        for (final MultiSampleEdge outgoingEdge : outgoingEdges) {
            final MultiDeBruijnVertex target = getEdgeTarget(outgoingEdge);
            if (target.getSuffix() == sequence[nextPos]) {
                if (alterTrees) {
                    // Only want to create a new tree if we walk into a node that meets our junction tree criteria
                    // NOTE: we do this deciding on our path
                    nodesHelper.addTreeIfNecessary(prevVertex);

                    // If this node has an out-degree > 1, add the edge we took to existing trees
                    if (outgoingEdges.size() != 1) {
                        nodesHelper.addEdgeToJunctionTreeNodes(outgoingEdge);
                    }
                }
                return target;
            }
        }
        return null;
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

            // Defense against loops, kill the attempt if we are stuck in a loop
            if (extraSequence.contains(vert)) {
                System.err.println("Dangling End recovery killed because of a loop (getReferencePathForwardFromKmer)");
                vert = null;
            }

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

    // Generate a SeqGraph that is special
    @Override
    public SeqGraph toSequenceGraph() {
        throw new UnsupportedOperationException("Cannot construct a sequence graph using JunctionTreeLinkedDeBruijnGraph");
    }

    /**
     * Print out the graph in the dot language for visualization for the graph after merging the graph into a seq graph.
     * Junction trees will naively add junction trees to the corresponding zipped SeqGraph node for its root vertex.
     *
     * NOTE: this is intended for debugging complex assembly graphs while preserving junction tree data.
     * @param destination File to write to
     */
    @VisibleForTesting
    public final void printSimplifiedGraph(final File destination, final int pruneFactor) {
        try (PrintStream stream = new PrintStream(new FileOutputStream(destination))) {
            printSimplifiedGraph(stream, true, pruneFactor);
        } catch ( final FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(destination.getAbsolutePath(), e);
        }
    }
    private final void printSimplifiedGraph(final PrintStream graphWriter, final boolean writeHeader, final int pruneFactor) {

        /// LOCAL CLASS TO MANAGE PRINTING THE JUNCTION TREES ONTO THE SEQ GRAPH
        class PrintingSeqGraph extends SeqGraph {
            private static final long serialVersionUID = 1l;

            private PrintingSeqGraph(int kmerSize) {
                super(kmerSize);
            }
            private Map<MultiDeBruijnVertex, SeqVertex> vertexToOrigionalSeqVertex;
            private Map<SeqVertex, SeqVertex> originalSeqVertexToMergedVerex;
            @Override
            // Extendable method intended to allow for adding extra material to the graph
            public List<String> getExtraGraphFileLines() {
                List<String> output = new ArrayList<>();
                for( Map.Entry<MultiDeBruijnVertex, ThreadingTree> entry : readThreadingJunctionTrees.entrySet()) {
                    // Resolving the chain of altered vertexes
                    SeqVertex mergedSeqVertex = originalSeqVertexToMergedVerex.get(vertexToOrigionalSeqVertex.get(entry.getKey()));

                    // adding the root node to the graph
                    output.add(String.format("\t%s -> %s ", mergedSeqVertex.toString(), entry.getValue().rootNode.getDotName()) +
                            String.format("[color=blue];"));
                    output.add(String.format("\t%s [shape=point];", entry.getValue().rootNode.getDotName()));

                    output.addAll(edgesForNodeRecursive(entry.getValue().rootNode));
                }
                return output;
            }
            // Helper class to be used for zipping linear chains for easy print outputting
            private Map<SeqVertex, SeqVertex> zipLinearChainsWithVertexMapping() {
                Map<SeqVertex, SeqVertex> vertexMapping = new HashMap<>();
                // create the list of start sites [doesn't modify graph yet]
                final Collection<SeqVertex> zipStarts = vertexSet().stream().filter(this::isLinearChainStart).collect(Collectors.toList());

                if ( zipStarts.isEmpty() ) // nothing to do, as nothing could start a chain
                {
                    return vertexMapping;
                }

                // At this point, zipStarts contains all of the vertices in this graph that might start some linear
                // chain of vertices.  We walk through each start, building up the linear chain of vertices and then
                // zipping them up with mergeLinearChain, if possible
                for ( final SeqVertex zipStart : zipStarts ) {
                    final LinkedList<SeqVertex> linearChain = traceLinearChain(zipStart);

                    SeqVertex mergedOne = mergeLinearChainVertex(linearChain);

                    // Add all of the mapped vertexes to the map.
                    if (mergedOne != null) {
                        for (SeqVertex v : linearChain) {
                            vertexMapping.put(v, mergedOne);
                        }
                    } else {
                        // otherwise we indicate that nothing has changed
                        for (SeqVertex v : linearChain) {
                            vertexMapping.put(v, v);
                        }
                    }
                }
                return vertexMapping;
            }
        };

        PrintingSeqGraph printingSeqGraph = new PrintingSeqGraph(kmerSize);

        final Map<MultiDeBruijnVertex, SeqVertex> vertexToOrigionalSeqVertex = new HashMap<>();

        // create all of the equivalent seq graph vertices
        for ( final MultiDeBruijnVertex dv : vertexSet() ) {
            final SeqVertex sv = new SeqVertex(dv.getAdditionalSequence(isSource(dv)));
            sv.setAdditionalInfo(dv.getAdditionalInfo());
            vertexToOrigionalSeqVertex.put(dv, sv);
            printingSeqGraph.addVertex(sv);
        }

        // walk through the nodes and connect them to their equivalent seq vertices
        for( final MultiSampleEdge e : edgeSet() ) {
            final SeqVertex seqInV = vertexToOrigionalSeqVertex.get(getEdgeSource(e));
            final SeqVertex seqOutV = vertexToOrigionalSeqVertex.get(getEdgeTarget(e));
            printingSeqGraph.addEdge(seqInV, seqOutV, new BaseEdge(e.isRef(), e.getMultiplicity()));
        }

        final Map<SeqVertex, SeqVertex> originalSeqVertexToMergedVerex = printingSeqGraph.zipLinearChainsWithVertexMapping();
        printingSeqGraph.originalSeqVertexToMergedVerex = originalSeqVertexToMergedVerex;
        printingSeqGraph.vertexToOrigionalSeqVertex = vertexToOrigionalSeqVertex;


        printingSeqGraph.printGraph(graphWriter, writeHeader, pruneFactor);
    }

    /**
     * Walk along the reference path in the graph and pull out the corresponding bases
     *
     * NOTE: this attempts to generate the longest sequence of refernce bases in the event that fromVertex or toVertex are non-unique
     *
     * @param fromVertex    starting vertex
     * @param toVertex      ending vertex
     * @param includeStart  should the starting vertex be included in the path
     * @param includeStop   should the ending vertex be included in the path
     * @return              byte[] array holding the reference bases, this can be null if there are no nodes between the starting and ending vertex (insertions for example)
     */
    @Override
    public byte[] getReferenceBytes( final MultiDeBruijnVertex fromVertex, final MultiDeBruijnVertex toVertex, final boolean includeStart, final boolean includeStop ) {
        Utils.nonNull(fromVertex, "Starting vertex in requested path cannot be null.");
        Utils.nonNull(toVertex, "From vertex in requested path cannot be null.");

        byte[] bytes = null;
        int fromIndex = referencePath.indexOf(fromVertex);
        int toIndex = referencePath.lastIndexOf(toVertex);

        if( includeStart ) {
            bytes = ArrayUtils.addAll(bytes, getAdditionalSequence(fromVertex, true));
        }
        for (int i = fromIndex + 1; i < toIndex; i++) {
            bytes = ArrayUtils.addAll(bytes, getAdditionalSequence(referencePath.get(i)));
        }

        if( includeStop ) {
            bytes = ArrayUtils.addAll(bytes, getAdditionalSequence(toVertex));
        }
        return bytes;
    }

    @Override
    // since we don't have to validate unique vertex merging we just find the vertex and pass
    protected MultiDeBruijnVertex getNextKmerVertexForChainExtension(final Kmer kmer, final boolean isRef, final MultiDeBruijnVertex prevVertex) {
        return kmerToVertexMap.get(kmer);
    }

    /**
     * Generate the junction trees and prune them (optionally printing the graph stages as output)
     *
     * @param debugGraphOutputPath path for graph output files
     * @param refHaplotype ref haplotype location
     */
    public void postProcessForHaplotypeFinding(final File debugGraphOutputPath, final Locatable refHaplotype) {
        annotateEdgesWithReferenceIndices();
        generateJunctionTrees();
        if (debugGraphTransformations) {
            printGraph(new File(debugGraphOutputPath, refHaplotype + "-sequenceGraph." + kmerSize + ".0.4.JT_unpruned.dot"), 10000);
        }
        pruneJunctionTrees(JunctionTreeKBestHaplotypeFinder.DEFAULT_MINIMUM_WEIGHT_FOR_JT_BRANCH_TO_NOT_BE_PRUNED);
        if (debugGraphTransformations) {
            printGraph(new File(debugGraphOutputPath, refHaplotype + "-sequenceGraph." + kmerSize + ".0.5.JT_pruned.dot"), 10000);
        }
    }

    private void annotateEdgesWithReferenceIndices() {
        final List<MultiDeBruijnVertex> referencePath = getReferencePath(TraversalDirection.downwards);
        MultiDeBruijnVertex lastVert = null;
        int refIndex = 0;
        for (MultiDeBruijnVertex nextVert : referencePath) {
            if (lastVert != null) {
                getEdge(lastVert, nextVert).addReferenceIndex(refIndex++);
            }
            lastVert = nextVert;
        }
    }


    // Generate threading trees
    public void generateJunctionTrees() {
        Utils.validate(alreadyBuilt, "Assembly graph has not been constructed, please call BuildGraphIfNecessary() before trying to thread reads to the graph");
        // Adding handle vertex to support symbolic end alleles
        addVertex(SYMBOLIC_END_VETEX);
        SYMBOLIC_END_EDGE = addEdge(getReferenceSinkVertex(), SYMBOLIC_END_VETEX);

        readThreadingJunctionTrees = new HashMap<>();
        pending.values().stream().flatMap(Collection::stream).forEach(this::threadSequenceForJuncitonTree);
    }

    /**
     * Traverse all of the junction trees in the graph and remove branches supported by < minimumEdgeWeight edges recursively.
     * This will also handle pruning of trees that are uninformative (i.e. empty roots).
     *
     * @param minimumEdgeWeight minimum edge weight below which branches are removed
     */
    public void pruneJunctionTrees(final int minimumEdgeWeight) {
        readThreadingJunctionTrees.forEach((key, value) -> value.getRootNode().pruneNode(minimumEdgeWeight));
        readThreadingJunctionTrees = Maps.filterValues( readThreadingJunctionTrees, ThreadingTree::isEmptyTree);
    }

    /**
     * Takes a contiguous sequence of kmers and threads it through the graph, generating and subsequently adding to
     * juncition trees at each relevant node.
     *
     * @param seqForKmers
     */
    private void threadSequenceForJuncitonTree(final SequenceForKmers seqForKmers) {
        // Maybe handle this differently, the reference junction tree should be held seperatedly from everything else.
        if (seqForKmers.isRef) {
            return;
        }

        // List we will use to keep track of sequences
        JunctionTreeThreadingHelper nodeHelper = new JunctionTreeThreadingHelper();

        // Find the first kmer in the read that exists on the graph
        final int startPos = findStartForJunctionThreading(seqForKmers);
        if ( startPos == -1 ) {
            return;
        }

        final MultiDeBruijnVertex startingVertex = kmerToVertexMap.get(new Kmer(seqForKmers.sequence, startPos, kmerSize));

        // loop over all of the bases in sequence, extending the graph by one base at each point, as appropriate
        MultiDeBruijnVertex lastVertex = startingVertex;
        boolean hasToRediscoverKmer = false;
        for ( int i = startPos + 1; i <= seqForKmers.stop - kmerSize; i++ ) {
            MultiDeBruijnVertex vertex;
            if (!hasToRediscoverKmer) {
                vertex = extendJunctionThreadingByOne(lastVertex, seqForKmers.sequence, i, nodeHelper, true);
            } else {
                Kmer kmer = new Kmer(seqForKmers.sequence, i, kmerSize);
                vertex = kmerToVertexMap.get(kmer);
            }

            // If we missed the vertex, attempt to recover the path from the graph if there is no ambiguity
            if (vertex == null && !hasToRediscoverKmer) {
                final Set<MultiSampleEdge> outgoingEdges = outgoingEdgesOf(lastVertex);
                MultiDeBruijnVertex tentativeVertex = null;
                if (outgoingEdges.size()==1) {
                    tentativeVertex = getEdgeTarget(outgoingEdges.stream().findFirst().get());
                }
                // Loop over the tentative path until the end of the read, we have walked over kmersize bases, or there is another missing path in the graph respectively
                for (int j = i+1; j <= seqForKmers.stop - kmerSize &&
                                j <= i + kmerSize &&
                                tentativeVertex != null; j++) {
                    tentativeVertex = extendJunctionThreadingByOne(tentativeVertex, seqForKmers.sequence, j, null, false);
                }
                // If tentativeVertex is not null then we will be able to successfully traverse the read through the graph
                if (tentativeVertex != null) {
                    vertex = getEdgeTarget(outgoingEdges.stream().findFirst().get());
                }
            }

            // If for whatever reason vertex = null, then we have fallen off the corrected graph so we don't update anything
            if (vertex != null) {
                lastVertex = vertex;
                hasToRediscoverKmer = false;
            } else {
                nodeHelper.clear();
                hasToRediscoverKmer = true;
            }
        }

        // As a final step, if the last vetex happens to be the ref-stop vertex then we want to append a symbolic node to the junciton trees
        if (lastVertex == getReferenceSinkVertex()) {
            nodeHelper.addEdgeToJunctionTreeNodes(SYMBOLIC_END_EDGE);
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

    public Optional<ThreadingTree> getJunctionTreeForNode(MultiDeBruijnVertex vertex) {
        return Optional.ofNullable(readThreadingJunctionTrees.get(vertex));
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
        private final ThreadingNode rootNode;
        private final MultiDeBruijnVertex treeVertex; // this reference exists purely for toString() reasons

        private ThreadingTree(final MultiDeBruijnVertex treeVertex) {
            this.rootNode = new ThreadingNode(null);
            this.treeVertex = treeVertex;
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

        // Returns the junction choices as a list of suffixes to follow (this is used for writing human readable tests)
        @VisibleForTesting
        List<String> getPathsPresentAsBaseChoiceStrings() {
            return rootNode.getEdgesThroughNode().stream()
                    .map(path -> path.stream()
                                    .map(edge -> getEdgeTarget(edge).getSuffix())
                                    .map(b -> Character.toString((char)b.byteValue()))
                                    .collect(Collectors.joining()))
                    .collect(Collectors.toList());

        }

        public String toString() {
            return "ThreadingTree_at_"+treeVertex.getSequenceString()+"_with_evidence_"+rootNode.count;
        }

        // getter for the root node, TODO to possibly be replaced when the graph gets hidden from prying eyes.
        public ThreadingNode getRootNode() {
            return rootNode;
        }
    }

    // Linked node object for storing tree topography
    public class ThreadingNode {
        private Map<MultiSampleEdge, ThreadingNode> childrenNodes;
        private final MultiSampleEdge prevEdge; // This may be null if this node corresponds to the root of the graph
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

        // Recursively prunes nodes based on the provided threshold, removing branches without sufficient support.
        private void pruneNode(final int threshold) {
            childrenNodes = Maps.filterValues( childrenNodes, node -> node.getEvidenceCount() >= threshold);
            childrenNodes.forEach((edge, node) -> node.pruneNode(threshold));
        }

        // Returns a unique name based on the memory id that conforms to the restrictions placed on .dot file nodes
        public String getDotName() {
            return "TreadingNode_" + Integer.toHexString(hashCode());
        }

        // Getter for external tools to access a node
        public Map<MultiSampleEdge, ThreadingNode> getChildrenNodes() {
            return Collections.unmodifiableMap(childrenNodes);
        }

        // Returns true if there are no paths emanating from this node
        public boolean hasNoEvidence() {
            return childrenNodes.isEmpty();
        }

        // Return the count of total evidence supporting this node in the tree
        public int getEvidenceCount() {
            return count;
        }

        // Checks if this node is the symbolic end by determining if its previous edge ends on the symbolic edge
        public boolean isSymbolicEnd() {
            return prevEdge == SYMBOLIC_END_EDGE;
        }
    }

    /**
     * A convenient helper class designed to pull much of the direct interaction with junction trees in the threading code
     * into a single place.
     */
    private class JunctionTreeThreadingHelper {
        final List<ThreadingNode> trackedNodes = new ArrayList<>(3);

        // Helper method used to determine whether a vertex meets the criteria to hold a junction tree
        // The current criteria, if any outgoing edge from a particular vertex leads to a vertex with inDegree > 1, then it warrants a tree. Or if it is the reference start vertex.
        // NOTE: this check is necessary to handle the edge cases that may arise when a vertex has multiple exit paths but happens to lead to a vetex that needs a junction tree
        private boolean vertexWarrantsJunctionTree(final MultiDeBruijnVertex vertex) {
            return outgoingEdgesOf(vertex).stream().anyMatch(edge -> inDegreeOf(getEdgeTarget(edge)) > 1);
        }

        // Helper method that adds a single edge to all of the nodes in nodesToExtend.
        private void addEdgeToJunctionTreeNodes(final MultiSampleEdge outgoingEdge) {
            List<ThreadingNode> newNodes = trackedNodes.stream().map(n -> n.addEdge(outgoingEdge)).collect(Collectors.toList());
            trackedNodes.clear();
            trackedNodes.addAll(newNodes);
        }

        // removes all of the nodes currently on the tree for safety if there are gaps in a reads traversal
        private void clear() {
            trackedNodes.clear();
        }

        // Adds a junction tree prevVertex if the vertex warrants a junciton tree
        private void addTreeIfNecessary(MultiDeBruijnVertex prevVertex) {
            if (vertexWarrantsJunctionTree(prevVertex)) {
                trackedNodes.add(readThreadingJunctionTrees.computeIfAbsent(prevVertex, k -> new ThreadingTree(prevVertex)).getAndIncrementRootNode());
            }
        }
    }
}
