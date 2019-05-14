package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.primitives.Bytes;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.Kmer;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.Utils;
import org.jgrapht.EdgeFactory;
import org.spark_project.jetty.util.ArrayQueue;

import java.io.File;
import java.util.*;

/**
 * A graph that contains base sequence at each node
 */
public final class SeqGraph extends BaseGraph<SeqVertex, BaseEdge> {

    Map<Kmer, Pair<SeqVertex, Integer>> kmerToSeqMap;

    private final Logger logger = LogManager.getLogger(SeqGraph.class);

    private static final long serialVersionUID = 1l;

    @Override
    public SeqGraph clone() {
        return (SeqGraph) super.clone();
    }

    /**
     * Edge factory that creates non-reference multiplicity 1 edges
     */
    private static class MyEdgeFactory implements EdgeFactory<SeqVertex, BaseEdge> {
        @Override
        public BaseEdge createEdge(final SeqVertex sourceVertex, final SeqVertex targetVertex) {
            return new BaseEdge(false, 1);
        }
    }

    private static final boolean PRINT_SIMPLIFY_GRAPHS = false;

    /**
     * How many cycles of the graph simplifications algorithms will we run before
     * thinking something has gone wrong and throw an exception?
     */
    private static final int MAX_REASONABLE_SIMPLIFICATION_CYCLES = 100;

    /**
     * Construct an empty SeqGraph where we'll add nodes based on a kmer size of kmer
     *
     * The kmer size is purely information.  It is useful when converting a Debruijn graph -> SeqGraph
     * for us to track the kmer used to make the transformation.
     *
     * @param kmer kmer
     */
    public SeqGraph(final int kmer) {
        super(kmer, new MyEdgeFactory());
    }

    /**
     * Simplify this graph, merging vertices together and restructuring the graph in an
     * effort to minimize the number of overall vertices in the graph without changing
     * in any way the sequences implied by a complex enumeration of all paths through the graph.
     */
    public void simplifyGraph() {
        simplifyGraph(Integer.MAX_VALUE);
    }

    @VisibleForTesting
    void simplifyGraph(final int maxCycles) {
        // start off with one round of zipping of chains for performance reasons
        zipLinearChains();

        SeqGraph prevGraph = null;
        for( int i = 0; i < maxCycles; i++ ) {
            if ( i > MAX_REASONABLE_SIMPLIFICATION_CYCLES ) {
                logger.warn("Infinite loop detected in simpliciation routines.  Writing current graph to debugMeMark.dot");
                printGraph(new File("debugMeMark.dot"), 0);
                throw new IllegalStateException("Infinite loop detected in simplification routines for kmer graph " + getKmerSize());
            }

            final boolean didSomeWork = simplifyGraphOnce(i);
            if ( ! didSomeWork )
                // no simplification algorithm could run, so stop
            {
                break;
            }

            // we get five cycles before we start looking for changes in the graph
            // by cloning ourselves and then checking for any changes
            if ( i > 5 ) {
                // the previous graph and this graph have the same structure, so the simplification
                // algorithms are looping endless between states.  Just break and consider ourselves done
                if ( prevGraph != null && graphEquals(prevGraph, this) ) {
                    break;
                }

                prevGraph = clone();
            }
        }
    }

    /**
     * Run one full cycle of the graph simplification algorithms
     * @return true if any algorithms said they did some simplification
     */
    private boolean simplifyGraphOnce(final int iteration) {
        //logger.info("simplifyGraph iteration " + i);
        // iterate until we haven't don't anything useful
        printGraphSimplification(new File("simplifyGraph." + iteration + ".1.dot"));
        boolean didSomeWork = false;
        didSomeWork |= new MergeDiamonds(this).transformUntilComplete();
        didSomeWork |= new MergeTails(this).transformUntilComplete();
        printGraphSimplification(new File("simplifyGraph." + iteration + ".2.diamonds_and_tails.dot"));

        didSomeWork |= new SplitCommonSuffices(this).transformUntilComplete();
        printGraphSimplification(new File("simplifyGraph." + iteration + ".3.split_suffix.dot"));
        didSomeWork |= new MergeCommonSuffices(this).transformUntilComplete();
        printGraphSimplification(new File("simplifyGraph." + iteration + ".4.merge_suffix.dot"));

        didSomeWork |= zipLinearChains();
        return didSomeWork;
    }

    /**
     * Print simplication step of this graph, if PRINT_SIMPLIFY_GRAPHS is enabled
     * @param file the destination for the graph DOT file
     */
    private void printGraphSimplification(final File file) {
        if ( PRINT_SIMPLIFY_GRAPHS ) {
            subsetToNeighbors(getReferenceSourceVertex(), 5).printGraph(file, 0);
        }
    }

    /**
     * Zip up all of the simple linear chains present in this graph.
     *
     * Merges together all pairs of vertices in the graph v1 -> v2 into a single vertex v' containing v1 + v2 sequence
     *
     * Only works on vertices where v1's only outgoing edge is to v2 and v2's only incoming edge is from v1.
     *
     * If such a pair of vertices is found, they are merged and the graph is update.  Otherwise nothing is changed.
     *
     * @return true if any such pair of vertices could be found, false otherwise
     */
    public boolean zipLinearChains() {
        // create the list of start sites [doesn't modify graph yet]
        final Collection<SeqVertex> zipStarts = new LinkedList<>();
        for ( final SeqVertex source : vertexSet() ) {
            if ( isLinearChainStart(source) ) {
                zipStarts.add(source);
            }
        }

        if ( zipStarts.isEmpty() ) // nothing to do, as nothing could start a chain
        {
            return false;
        }

        // At this point, zipStarts contains all of the vertices in this graph that might start some linear
        // chain of vertices.  We walk through each start, building up the linear chain of vertices and then
        // zipping them up with mergeLinearChain, if possible
        boolean mergedOne = false;
        for ( final SeqVertex zipStart : zipStarts ) {
            final LinkedList<SeqVertex> linearChain = traceLinearChain(zipStart);

            // merge the linearized chain, recording if we actually did some useful work
            mergedOne |= mergeLinearChain(linearChain);
        }

        return mergedOne;
    }

    /**
     * Is source vertex potentially a start of a linear chain of vertices?
     *
     * We are a start of a zip chain if our out degree is 1 and either the
     * the vertex has no incoming connections or 2 or more (we must start a chain) or
     * we have exactly one incoming vertex and that one has out-degree > 1 (i.e., source's incoming
     * vertex couldn't be a start itself
     *
     * @param source a non-null vertex
     * @return true if source might start a linear chain
     */
    private boolean isLinearChainStart(final SeqVertex source) {
        return outDegreeOf(source) == 1
                && ( inDegreeOf(source) != 1
                     || outDegreeOf(incomingVerticesOf(source).iterator().next()) > 1 );
    }

    /**
     * Get all of the vertices in a linear chain of vertices starting at zipStart
     *
     * Build a list of vertices (in order) starting from zipStart such that each sequential pair of vertices
     * in the chain A and B can be zipped together.
     *
     * @param zipStart a vertex that starts a linear chain
     * @return a list of vertices that comprise a linear chain starting with zipStart.  The resulting
     *         list will always contain at least zipStart as the first element.
     */
    private LinkedList<SeqVertex> traceLinearChain(final SeqVertex zipStart) {
        final LinkedList<SeqVertex> linearChain = new LinkedList<>();
        linearChain.add(zipStart);

        boolean lastIsRef = isReferenceNode(zipStart); // remember because this calculation is expensive
        SeqVertex last = zipStart;
        while (true) {
            if ( outDegreeOf(last) != 1 )
                // cannot extend a chain from last if last has multiple outgoing branches
            {
                break;
            }

            // there can only be one (outgoing edge of last) by contract
            final SeqVertex target = getEdgeTarget(outgoingEdgeOf(last));

            if ( inDegreeOf(target) != 1 || last.equals(target) )
                // cannot zip up a target that has multiple incoming nodes or that's a cycle to the last node
            {
                break;
            }

            final boolean targetIsRef = isReferenceNode(target);
            if ( lastIsRef != targetIsRef ) // both our isRef states must be equal
            {
                break;
            }

            linearChain.add(target); // extend our chain by one

            // update our last state to be the current state, and continue
            last = target;
            lastIsRef = targetIsRef;
        }

        return linearChain;
    }

    /**
     * Merge a linear chain of vertices into a single combined vertex, and update this graph to such that
     * the incoming edges into the first element of the linearChain and the outgoing edges from linearChain.getLast()
     * all point to this new combined vertex.
     *
     * @param linearChain a non-empty chain of vertices that can be zipped up into a single vertex
     * @return true if we actually merged at least two vertices together
     */
    private boolean mergeLinearChain(final LinkedList<SeqVertex> linearChain) {
        Utils.validateArg(!linearChain.isEmpty(), () -> "BUG: cannot have linear chain with 0 elements but got " + linearChain);

        final SeqVertex first = linearChain.getFirst();
        final SeqVertex last = linearChain.getLast();

        if ( first == last ) {
            return false; // only one element in the chain, cannot be extended
        }

        // create the combined vertex, and add it to the graph
        // TODO -- performance problem -- can be optimized if we want

        final SeqVertex addedVertex = mergeLinearChainVertices(linearChain);
        addVertex(addedVertex);

        // update the incoming and outgoing edges to point to the new vertex
        for( final BaseEdge edge : outgoingEdgesOf(last) ) { addEdge(addedVertex, getEdgeTarget(edge), edge.copy()); }
        for( final BaseEdge edge : incomingEdgesOf(first) )  { addEdge(getEdgeSource(edge), addedVertex, edge.copy()); }

        removeAllVertices(linearChain);
        return true;
    }

    private static SeqVertex mergeLinearChainVertices(final Iterable<SeqVertex> vertices) {
        final List<byte[]> seqs = new LinkedList<>();
        for ( final SeqVertex v : vertices ) {
            seqs.add(v.getSequence());
        }
        final byte[] seqsCat = Bytes.concat(seqs.toArray(new byte[][]{}));
        return new SeqVertex( seqsCat );
    }


    /**
     * Method that that takes a sequence graph and fills in kmerToSeqMap with a map of Kmer to first base in the graph.
     *
     * This is accomplished by traversing the graph depth first backwards from each outDegree 0 node, splitting its running
     * kmer each time it reaches an edge
     */
    public void fillKmerMap() {
        // TODO get sinks will need to be altered to handle pathological loop structures
        kmerToSeqMap = new HashMap<>();

        for (SeqVertex sinkVertex : getSinks()) {
            KmerCrawler rootCrawler = new KmerCrawler(kmerSize, sinkVertex);
            Queue<KmerCrawler> runningCrawlers = new ArrayQueue<>();
            runningCrawlers.add(rootCrawler);

            while (!runningCrawlers.isEmpty()) {
                KmerCrawler topCrawler = runningCrawlers.poll();
                topCrawler.traverse(runningCrawlers, kmerToSeqMap, this);
            }
        }
    }


    /**
     * Helper object used to build kmer maps
     */
    static class KmerCrawler {
        SeqVertex startVertex;
        int currentSeqPos;

        byte[] rollingArray;
        int arrayPos;

        private KmerCrawler(final int kmerSize, final SeqVertex startVertex) {
            this.startVertex = startVertex;
            this.currentSeqPos = startVertex.sequence.length;
            this.rollingArray = new byte[kmerSize];
            this.arrayPos = kmerSize;
        }

        // Internal constructor intended for copying a kmerCrawler
        private KmerCrawler(final byte[] original, final int arrayPos, final SeqVertex startVertex) {
            this.startVertex = startVertex;
            this.currentSeqPos = startVertex.sequence.length;
            this.rollingArray = Arrays.copyOf(original, original.length);
            this.arrayPos = arrayPos;
        }

        private KmerCrawler splitCrawler(final SeqVertex startVertex) {
            return new KmerCrawler(this.rollingArray, arrayPos, startVertex);
        }

        /**
         * Walk all the way to an a start vertex, adding to the mapToOutput map the appropriate kmer at each stage.
         *
         * When an inDegree > 1 edge is encountered, we spawn off a copy of this KmerCralwer object for each edge beyond
         * the first pinting at the last base of each of the root vertexes. We then add each of those crawlers to the provided
         * waitingCrawlers queue.
         *
         * Upon a complete call of traverse we expect kmers to be generated for the entirety of one path forward from the
         * starting vertex already associated with this object. Kmers will be constructed from any suffix bases already
         * present in this crawler (which may have been created when the traversal split) and the bases seen from traversal.
         * If there are not enough accumulated bases to to construct a complete kmer then nothing is done.
         *
         * Each kmer has its associated starting position (SeqVertex and StartingBase index) stored in the provided
         * mapToOutput object.
         *
         * NOTE: this method makes the assumption that the sequence graph will have been constructed from unique kmers,
         *       thus making the assumption that there is only one graph position for each kmer valid.
         *
         * @param waitingCrawlers    Queue of spawned off crawlers that have not completed been traversed yet.
         * @param mapToOutput        Map of kmer, to start position to have resutls added into.
         */
        private void traverse(Queue<KmerCrawler> waitingCrawlers, Map<Kmer, Pair<SeqVertex, Integer>> mapToOutput, SeqGraph graph) {
            while (traverseOneEdge(waitingCrawlers, mapToOutput, graph)) {

            }
        }

        private boolean traverseOneEdge(Queue<KmerCrawler> waitingCrawlers, Map<Kmer, Pair<SeqVertex, Integer>> mapToOutput, SeqGraph graph) {
            // Check if we are at the end of the seq vertex
            if (--currentSeqPos < 0) {
                Set<SeqVertex> incoming = graph.incomingVerticesOf(startVertex);
                if (incoming.size() == 0) {
                    // Have hit the end of a source vertex and no longer have any inbound paths to traverse
                    return false;
                }
                final boolean[] seenFirst = {false};
                //TODO methodify this?
                incoming.forEach(vertex -> {
                    if (seenFirst[0]) {
                        waitingCrawlers.add(splitCrawler(vertex));
                    } else {
                        this.startVertex = vertex;
                        this.currentSeqPos = vertex.sequence.length;
                    }
                    seenFirst[0] = true;
                });

            // If we don't have to make a split, simply update the next step
            } else {
                // Update the position of the array (noting the rolling feature)
                rollingArray[(--arrayPos) % rollingArray.length] = startVertex.sequence[currentSeqPos];
                Kmer currentKmer = getKmerForRollingArray(rollingArray, arrayPos);

                // If we have traversed far enough to produce a full kmer
                if (currentKmer != null) {
                    // If the map already has this kmer, then we terminate because we have folded back in on an extant branch
                    if (mapToOutput.containsKey(currentKmer)) {
                        return false;
                    }

                    mapToOutput.put(currentKmer, new ImmutablePair<>(startVertex, currentSeqPos));
                }

            }
            return true;
        }


        /**
         * Takes the provided kmer array and the current (possibly negative) position and generates a Kmer object from it
         *
         * @param kmerArray  Array containing bases to be kmerized
         * @param currentPos Offset corresponding to the modulo offset for the start of the kmer in the array
         * @return A kmer corresponding to the bases wrapping starting from mod currentPos (or null if currentPos is > 0)
         */
        @VisibleForTesting
        public static Kmer getKmerForRollingArray(byte[] kmerArray, int currentPos) {
            // have not finished accumulating the first complete kmer into the array
            if (currentPos > 0) {
                return null;
            }

            // Copy kmerArray bases into a new array to use to build the output
            int offset = Math.floorMod(currentPos, kmerArray.length);
            byte[] correctedKmer = new byte[kmerArray.length];
            for(int i = 0; i < kmerArray.length; i++) {
                correctedKmer[i] = kmerArray[Math.floorMod(offset + i, kmerArray.length)];
            }
            return new Kmer(correctedKmer);
        }
    }

}
