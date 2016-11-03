package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.Kmer;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.BaseGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.KmerSearchableGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqGraph;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.SWPairwiseAlignment;
import org.jgrapht.EdgeFactory;

import java.io.File;
import java.util.*;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

/**
 * Note: not final but only intendent to be subclassed for testing.
 */
public class ReadThreadingGraph extends BaseGraph<MultiDeBruijnVertex, MultiSampleEdge> implements KmerSearchableGraph<MultiDeBruijnVertex,MultiSampleEdge> {

    private static final Logger logger = LogManager.getLogger(ReadThreadingGraph.class);

    private static final String ANONYMOUS_SAMPLE = "XXX_UNNAMED_XXX";
    private static final boolean WRITE_GRAPH = false;
    private static final boolean DEBUG_NON_UNIQUE_CALC = false;

    private static final int MAX_CIGAR_COMPLEXITY = 3;
    private static final long serialVersionUID = 1l;
    private int maxMismatchesInDanglingHead = -1;

    private boolean alreadyBuilt;

    private boolean startThreadingOnlyAtExistingVertex = false;

    /** for debugging info printing */
    private static int counter = 0;

    /**
     * Sequences added for read threading before we've actually built the graph
     */
    private final Map<String, List<SequenceForKmers>> pending = new LinkedHashMap<>();

    /**
     * A set of non-unique kmers that cannot be used as merge points in the graph
     */
    private Set<Kmer> nonUniqueKmers;

    /**
     * A map from kmers -> their corresponding vertex in the graph
     */
    private final Map<Kmer, MultiDeBruijnVertex> uniqueKmers = new LinkedHashMap<>();

    private final boolean debugGraphTransformations;
    private final byte minBaseQualityToUseInAssembly;

    private static final boolean INCREASE_COUNTS_BACKWARDS = true;
    private boolean increaseCountsThroughBranches = false; // this may increase the branches without bounds

    // --------------------------------------------------------------------------------
    // state variables, initialized in resetToInitialState()
    // --------------------------------------------------------------------------------
    private Kmer refSource;

    /**
     * Constructs an empty read-threading-grpah provided the kmerSize.
     * @param kmerSize 1 or greater.
     *
     * @throws IllegalArgumentException if (@code kmerSize) < 1.
     */
    public ReadThreadingGraph(final int kmerSize) {
        this(kmerSize, false, (byte)6, 1);
    }

    /**
     * Constructs a read-threading-graph for a string representation.
     *
     * <p>
     *     Note: only used for testing.
     * </p>
     * @param s the string representation of the graph {@code null}.
     */
    @VisibleForTesting
    protected ReadThreadingGraph(final int kmerSizeFromString, final EdgeFactory<MultiDeBruijnVertex, MultiSampleEdge> edgeFactory) {
        super(kmerSizeFromString, new MyEdgeFactory(1));
        debugGraphTransformations = false;
        minBaseQualityToUseInAssembly = 0;
    }

    @VisibleForTesting
    protected void setAlreadyBuilt() {
        alreadyBuilt = true;
    }

    @VisibleForTesting
    void setMaxMismatchesInDanglingHead(final int maxMismatchesInDanglingHead) {
        this.maxMismatchesInDanglingHead = maxMismatchesInDanglingHead;
    }

    /**
     * Return the collection of outgoing vertices that expand this vertex with a particular base.
     *
     * @param v original vertex.
     * @param b expanding base.
     * @return never null, but perhaps an empty set. You cannot assume that you can modify the result.
     */
    @VisibleForTesting
    Set<MultiDeBruijnVertex> getNextVertices(final MultiDeBruijnVertex v, final byte b) {
        Utils.nonNull(v, "the input vertex cannot be null");
        Utils.validateArg(vertexSet().contains(v), "the vertex must be present in the graph");
        final List<MultiDeBruijnVertex> result = new LinkedList<>();
        for (final MultiDeBruijnVertex w : outgoingVerticesOf(v)) {
            if (w.getSuffix() == b) {
                result.add(w);
            }
        }
        switch (result.size()) {
            case 0: return Collections.emptySet();
            case 1: return Collections.singleton(result.get(0));
            default:
                    return new HashSet<>(result);
        }
    }

    /**
     * Create a new ReadThreadingAssembler using kmerSize for matching
     * @param kmerSize must be >= 1
     */
    ReadThreadingGraph(final int kmerSize, final boolean debugGraphTransformations, final byte minBaseQualityToUseInAssembly, final int numPruningSamples) {
        super(kmerSize, new MyEdgeFactory(numPruningSamples));

        Utils.validateArg( kmerSize > 0, () -> "bad minkKmerSize " + kmerSize);

        this.debugGraphTransformations = debugGraphTransformations;
        this.minBaseQualityToUseInAssembly = minBaseQualityToUseInAssembly;

        resetToInitialState();
    }

    /**
     * Reset this assembler to its initial state, so we can create another assembly with a different set of reads
     */
    private void resetToInitialState() {
        pending.clear();
        nonUniqueKmers = null;
        uniqueKmers.clear();
        refSource = null;
        alreadyBuilt = false;
    }

    /**
     * Add the all bases in sequence to the graph
     * @param sequence a non-null sequence
     * @param isRef is this the reference sequence?
     */
    final void addSequence(final byte[] sequence, final boolean isRef) {
        addSequence("anonymous", sequence, isRef);
    }

    /**
     * Add all bases in sequence to this graph
     *
     * @see #addSequence(String, String, byte[], int, int, int, boolean) for full information
     */
    public final void addSequence(final String seqName, final byte[] sequence, final boolean isRef) {
        addSequence(seqName, sequence, 1, isRef);
    }

    /**
     * Add all bases in sequence to this graph
     *
     * @see #addSequence(String, String, byte[], int, int, int, boolean) for full information
     */
    public final void addSequence(final String seqName, final byte[] sequence, final int count, final boolean isRef) {
        addSequence(seqName, ANONYMOUS_SAMPLE, sequence, 0, sequence.length, count, isRef);
    }

    /**
     * Add bases in sequence to this graph
     *
     * @param seqName a useful seqName for this read, for debugging purposes
     * @param sequence non-null sequence of bases
     * @param start the first base offset in sequence that we should use for constructing the graph using this sequence, inclusive
     * @param stop the last base offset in sequence that we should use for constructing the graph using this sequence, exclusive
     * @param count the representative count of this sequence (to use as the weight)
     * @param isRef is this the reference sequence.
     */
    private void addSequence(final String seqName, final String sampleName, final byte[] sequence, final int start, final int stop, final int count, final boolean isRef) {
        // note that argument testing is taken care of in SequenceForKmers
        if ( alreadyBuilt ) {
            throw new IllegalStateException("Graph already built");
        }

        // get the list of sequences for this sample
        List<SequenceForKmers> sampleSequences = pending.get(sampleName);
        if ( sampleSequences == null ) { // need to create
            sampleSequences = new LinkedList<>();
            pending.put(sampleName, sampleSequences);
        }

        // add the new sequence to the list of sequences for sample
        sampleSequences.add(new SequenceForKmers(seqName, sequence, start, stop, count, isRef));
    }

    /**
     * Thread sequence seqForKmers through the current graph, updating the graph as appropriate
     * @param seqForKmers a non-null sequence
     */
    private void threadSequence(final SequenceForKmers seqForKmers) {
        final int uniqueStartPos = findStart(seqForKmers);
        if ( uniqueStartPos == -1 ) {
            return;
        }

        final MultiDeBruijnVertex startingVertex = getOrCreateKmerVertex(seqForKmers.sequence, uniqueStartPos);

        // increase the counts of all edges incoming into the starting vertex supported by going back in sequence
        if (INCREASE_COUNTS_BACKWARDS) {
            increaseCountsInMatchedKmers(seqForKmers, startingVertex, startingVertex.getSequence(), kmerSize - 2);
        }

        if ( debugGraphTransformations ) {
            startingVertex.addRead(seqForKmers.name);
        }

        // keep track of information about the reference source
        if ( seqForKmers.isRef ) {
            if ( refSource != null ) {
                throw new IllegalStateException("Found two refSources! prev: " + refSource + ", new: " + startingVertex);
            }
            refSource = new Kmer(seqForKmers.sequence, seqForKmers.start, kmerSize);
        }

        // loop over all of the bases in sequence, extending the graph by one base at each point, as appropriate
        MultiDeBruijnVertex vertex = startingVertex;
        for ( int i = uniqueStartPos + 1; i <= seqForKmers.stop - kmerSize; i++ ) {
            vertex = extendChainByOne(vertex, seqForKmers.sequence, i, seqForKmers.count, seqForKmers.isRef);
            if ( debugGraphTransformations ) {
                vertex.addRead(seqForKmers.name);
            }
        }
    }

    /**
     * Find vertex and its position in seqForKmers where we should start assembling seqForKmers
     *
     * @param seqForKmers the sequence we want to thread into the graph
     * @return the position of the starting vertex in seqForKmer, or -1 if it cannot find one
     */
    private int findStart(final SequenceForKmers seqForKmers) {
        if ( seqForKmers.isRef ) {
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
     * Checks whether a kmer can be the threading start based on the current threading start location policy.
     *
     * @see #setThreadingStartOnlyAtExistingVertex(boolean)
     * @see #getThreadingStartOnlyAtExistingVertex()
     *
     * @param kmer the query kmer.
     * @return {@code true} if we can start thread the sequence at this kmer, {@code false} otherwise.
     */
    private boolean isThreadingStart(final Kmer kmer) {
        Utils.nonNull(kmer);
        return startThreadingOnlyAtExistingVertex ? uniqueKmers.containsKey(kmer) : !nonUniqueKmers.contains(kmer);
    }

    /**
     * Changes the threading start location policy.
     *
     * @param value  {@code true} if threading will start only at existing vertices in the graph, {@code false} if
     *  it can start at any unique kmer.
     */
    public final void setThreadingStartOnlyAtExistingVertex(final boolean value) {
        startThreadingOnlyAtExistingVertex = value;
    }

    /**
     * Build the read threaded assembly graph if it hasn't already been constructed from the sequences that have
     * been added to the graph.
     */
    public final void buildGraphIfNecessary() {
        if ( alreadyBuilt ) {
            return;
        }

        // determine the kmer size we'll use, and capture the set of nonUniques for that kmer size
        final NonUniqueResult result = determineKmerSizeAndNonUniques(kmerSize, kmerSize);
        nonUniqueKmers = result.nonUniques;

        if ( DEBUG_NON_UNIQUE_CALC ) {
            logger.info("using " + kmerSize + " kmer size for this assembly with the following non-uniques");
        }

        // go through the pending sequences, and add them to the graph
        for ( final List<SequenceForKmers> sequencesForSample : pending.values() ) {
            for ( final SequenceForKmers sequenceForKmers : sequencesForSample ) {
                threadSequence(sequenceForKmers);
                if ( WRITE_GRAPH ) {
                    printGraph(new File("threading." + counter++ + '.' + sequenceForKmers.name.replace(" ", "_") + ".dot"), 0);
                }
            }

            // flush the single sample edge values from the graph
            for ( final MultiSampleEdge e : edgeSet() ) {
                e.flushSingleSampleMultiplicity();
            }
        }

        // clear
        pending.clear();
        alreadyBuilt = true;
        for (final MultiDeBruijnVertex v : uniqueKmers.values()) {
            v.setAdditionalInfo(v.getAdditionalInfo() + '+');
        }
    }


    @Override
    public boolean removeVertex(final MultiDeBruijnVertex V) {
        final boolean result = super.removeVertex(V);
        if (result) {
            final byte[] sequence = V.getSequence();
            final Kmer kmer = new Kmer(sequence);
            uniqueKmers.remove(kmer);
        }
        return result;
    }

    @Override
    public void removeSingletonOrphanVertices() {
        // Run through the graph and clean up singular orphaned nodes
        final Collection<MultiDeBruijnVertex> verticesToRemove = new LinkedList<>();
        for( final MultiDeBruijnVertex v : vertexSet() ) {
            if( inDegreeOf(v) == 0 && outDegreeOf(v) == 0 ) {
                verticesToRemove.add(v);
            }
        }
        this.removeVertex(null);
        removeAllVertices(verticesToRemove);
    }

    /**
     * Does the graph not have enough complexity?  We define low complexity as a situation where the number
     * of non-unique kmers is more than 20% of the total number of kmers.
     *
     * @return true if the graph has low complexity, false otherwise
     */
    public boolean isLowComplexity() {
        return nonUniqueKmers.size() * 4 > uniqueKmers.size();
    }

    @Override
    public ReadThreadingGraph clone() {
        return (ReadThreadingGraph) super.clone();
    }

    @VisibleForTesting
    void setIncreaseCountsThroughBranches(final boolean increaseCountsThroughBranches) {
        this.increaseCountsThroughBranches = increaseCountsThroughBranches;
    }

    /**
     * Edge factory that encapsulates the numPruningSamples assembly parameter
     */
    private static final class MyEdgeFactory implements EdgeFactory<MultiDeBruijnVertex, MultiSampleEdge> {
        final int numPruningSamples;

        private MyEdgeFactory(final int numPruningSamples) {
            this.numPruningSamples = numPruningSamples;
        }

        @Override
        public MultiSampleEdge createEdge(final MultiDeBruijnVertex sourceVertex, final MultiDeBruijnVertex targetVertex) {
            return new MultiSampleEdge(false, 1, numPruningSamples);
        }

        public MultiSampleEdge createEdge(final boolean isRef, final int multiplicity) {
            return new MultiSampleEdge(isRef, multiplicity, numPruningSamples);
        }
    }

    /**
     * Class to keep track of the important dangling chain merging data
     */
    static final class DanglingChainMergeHelper {
        final List<MultiDeBruijnVertex> danglingPath;
        final List<MultiDeBruijnVertex> referencePath;
        final byte[] danglingPathString;
        final byte[] referencePathString;
        final Cigar cigar;

        DanglingChainMergeHelper(final List<MultiDeBruijnVertex> danglingPath,
                                        final List<MultiDeBruijnVertex> referencePath,
                                        final byte[] danglingPathString,
                                        final byte[] referencePathString,
                                        final Cigar cigar) {
            this.danglingPath = danglingPath;
            this.referencePath = referencePath;
            this.danglingPathString = danglingPathString;
            this.referencePathString = referencePathString;
            this.cigar = cigar;
        }
    }

    /** structure that keeps track of the non-unique kmers for a given kmer size */
    private static final class NonUniqueResult {
        final Set<Kmer> nonUniques;

        private NonUniqueResult(final Set<Kmer> nonUniques) {
            this.nonUniques = nonUniques;
        }
    }

    /**
     * Keeps track of the information needed to add a sequence to the read threading assembly graph
     */
    static final class SequenceForKmers {
        final String name;
        final byte[] sequence;
        final int start;
        final int stop;
        final int count;
        final boolean isRef;

        /**
         * Create a new sequence for creating kmers
         */
        SequenceForKmers(final String name, final byte[] sequence, final int start, final int stop, final int count, final boolean ref) {
            Utils.nonNull(sequence, "Sequence is null ");
            Utils.validateArg( start >= 0, () -> "Invalid start " + start);
            Utils.validateArg( stop >= start, () -> "Invalid stop " + stop);
            Utils.validateArg( count > 0, "Invalid count " + count);

            this.name = name;
            this.sequence = sequence;
            this.start = start;
            this.stop = stop;
            this.count = count;
            this.isRef = ref;
        }
    }

    /**
     * Try to recover dangling tails
     *
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
     */
    public void recoverDanglingTails(final int pruneFactor, final int minDanglingBranchLength) {
        Utils.validateArg(pruneFactor >= 0, () -> "pruneFactor must be non-negative but was " + pruneFactor);
        Utils.validateArg(minDanglingBranchLength >= 0, () -> "minDanglingBranchLength must be non-negative but was " + minDanglingBranchLength);

        if ( ! alreadyBuilt ) {
            throw new IllegalStateException("recoverDanglingTails requires the graph be already built");
        }

        int attempted = 0;
        int nRecovered = 0;
        for ( final MultiDeBruijnVertex v : vertexSet() ) {
            if ( outDegreeOf(v) == 0 && ! isRefSink(v) ) {
                attempted++;
                nRecovered += recoverDanglingTail(v, pruneFactor, minDanglingBranchLength);
            }
        }

        logger.debug(String.format("Recovered %d of %d dangling tails", nRecovered, attempted));
    }


    /**
     * Try to recover dangling heads
     *
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
     */
    public void recoverDanglingHeads(final int pruneFactor, final int minDanglingBranchLength) {
        Utils.validateArg(pruneFactor >= 0, () -> "pruneFactor must be non-negative but was " + pruneFactor);
        Utils.validateArg(minDanglingBranchLength >= 0, () -> "minDanglingBranchLength must be non-negative but was " + minDanglingBranchLength);
        if ( ! alreadyBuilt ) {
            throw new IllegalStateException("recoverDanglingHeads requires the graph be already built");
        }

        // we need to build a list of dangling heads because that process can modify the graph (and otherwise generate
        // a ConcurrentModificationException if we do it while iterating over the vertexes)
        final Collection<MultiDeBruijnVertex> danglingHeads = vertexSet().stream()
                .filter(v -> inDegreeOf(v) == 0 && ! isRefSource(v))
                .collect(Collectors.toList());

        // now we can try to recover the dangling heads
        int attempted = 0;
        int nRecovered = 0;
        for ( final MultiDeBruijnVertex v : danglingHeads ) {
            attempted++;
            nRecovered += recoverDanglingHead(v, pruneFactor, minDanglingBranchLength);
        }

        logger.debug(String.format("Recovered %d of %d dangling heads", nRecovered, attempted));
    }

    /**
     * Attempt to attach vertex with out-degree == 0 to the graph
     *
     * @param vertex the vertex to recover
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
     * @return 1 if we successfully recovered the vertex and 0 otherwise
     */
    private int recoverDanglingTail(final MultiDeBruijnVertex vertex, final int pruneFactor, final int minDanglingBranchLength) {
        if ( outDegreeOf(vertex) != 0 ) {
            throw new IllegalStateException("Attempting to recover a dangling tail for " + vertex + " but it has out-degree > 0");
        }

        // generate the CIGAR string from Smith-Waterman between the dangling tail and reference paths
        final DanglingChainMergeHelper danglingTailMergeResult = generateCigarAgainstDownwardsReferencePath(vertex, pruneFactor, minDanglingBranchLength);

        // if the CIGAR is too complex (or couldn't be computed) then we do not allow the merge into the reference path
        if ( danglingTailMergeResult == null || ! cigarIsOkayToMerge(danglingTailMergeResult.cigar, false, true) ) {
            return 0;
        }

        // merge
        return mergeDanglingTail(danglingTailMergeResult);
    }

    /**
     * Attempt to attach vertex with in-degree == 0, or a vertex on its path, to the graph
     *
     * @param vertex the vertex to recover
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
     * @return 1 if we successfully recovered a vertex and 0 otherwise
     */
    private int recoverDanglingHead(final MultiDeBruijnVertex vertex, final int pruneFactor, final int minDanglingBranchLength) {
        if ( inDegreeOf(vertex) != 0 ) {
            throw new IllegalStateException("Attempting to recover a dangling head for " + vertex + " but it has in-degree > 0");
        }

        // generate the CIGAR string from Smith-Waterman between the dangling tail and reference paths
        final DanglingChainMergeHelper danglingHeadMergeResult = generateCigarAgainstUpwardsReferencePath(vertex, pruneFactor, minDanglingBranchLength);

        // if the CIGAR is too complex (or couldn't be computed) then we do not allow the merge into the reference path
        if ( danglingHeadMergeResult == null || ! cigarIsOkayToMerge(danglingHeadMergeResult.cigar, true, false) ) {
            return 0;
        }

        // merge
        return mergeDanglingHead(danglingHeadMergeResult);
    }

    /**
     * Determine whether the provided cigar is okay to merge into the reference path
     *
     * @param cigar    the cigar to analyze
     * @param requireFirstElementM if true, require that the first cigar element be an M operator in order for it to be okay
     * @param requireLastElementM  if true, require that the last cigar element be an M operator in order for it to be okay
     * @return true if it's okay to merge, false otherwise
     */
    @VisibleForTesting
    static boolean cigarIsOkayToMerge(final Cigar cigar, final boolean requireFirstElementM, final boolean requireLastElementM) {
        final List<CigarElement> elements = cigar.getCigarElements();
        final int numElements = elements.size();

        // don't allow more than a couple of different ops
        if ( numElements == 0 || numElements > MAX_CIGAR_COMPLEXITY ) {
            return false;
        }

        // the first element must be an M
        if ( requireFirstElementM && elements.get(0).getOperator() != CigarOperator.M ) {
            return false;
        }

        // the last element must be an M
        if ( requireLastElementM && elements.get(numElements - 1).getOperator() != CigarOperator.M ) {
            return false;
        }

        // note that there are checks for too many mismatches in the dangling branch later in the process

        return true;
    }

    /**
     * Actually merge the dangling tail if possible
     *
     * @param danglingTailMergeResult   the result from generating a Cigar for the dangling tail against the reference
     * @return 1 if merge was successful, 0 otherwise
     */
    @VisibleForTesting
    final int mergeDanglingTail(final DanglingChainMergeHelper danglingTailMergeResult) {

        final List<CigarElement> elements = danglingTailMergeResult.cigar.getCigarElements();
        final CigarElement lastElement = elements.get(elements.size() - 1);
        Utils.validateArg( lastElement.getOperator() == CigarOperator.M, "The last Cigar element must be an M");

        final int lastRefIndex = danglingTailMergeResult.cigar.getReferenceLength() - 1;
        final int matchingSuffix = Math.min(longestSuffixMatch(danglingTailMergeResult.referencePathString, danglingTailMergeResult.danglingPathString, lastRefIndex), lastElement.getLength());
        if ( matchingSuffix == 0 ) {
            return 0;
        }

        final int altIndexToMerge = Math.max(danglingTailMergeResult.cigar.getReadLength() - matchingSuffix - 1, 0);

        // there is an important edge condition that we need to handle here: Smith-Waterman correctly calculates that there is a
        // deletion, that deletion is left-aligned such that the LCA node is part of that deletion, and the rest of the dangling
        // tail is a perfect match to the suffix of the reference path.  In this case we need to push the reference index to merge
        // down one position so that we don't incorrectly cut a base off of the deletion.
        final boolean firstElementIsDeletion = elements.get(0).getOperator() == CigarOperator.D;
        final boolean mustHandleLeadingDeletionCase =  firstElementIsDeletion && (elements.get(0).getLength() + matchingSuffix == lastRefIndex + 1);
        final int refIndexToMerge = lastRefIndex - matchingSuffix + 1 + (mustHandleLeadingDeletionCase ? 1 : 0);

        // another edge condition occurs here: if Smith-Waterman places the whole tail into an insertion then it will try to
        // merge back to the LCA, which results in a cycle in the graph.  So we do not want to merge in such a case.
        if ( refIndexToMerge == 0 ) {
            return 0;
        }

        // it's safe to merge now
        addEdge(danglingTailMergeResult.danglingPath.get(altIndexToMerge), danglingTailMergeResult.referencePath.get(refIndexToMerge), ((MyEdgeFactory)getEdgeFactory()).createEdge(false, 1));

        return 1;
    }


    /**
     * calculates the longest suffix match between a sequence and a smaller kmer
     *
     * @param seq         the (reference) sequence
     * @param kmer        the smaller kmer sequence
     * @param seqStart    the index (inclusive) on seq to start looking backwards from
     * @return the longest matching suffix
     */
    @VisibleForTesting
    static int longestSuffixMatch(final byte[] seq, final byte[] kmer, final int seqStart) {
        for ( int len = 1; len <= kmer.length; len++ ) {
            final int seqI = seqStart - len + 1;
            final int kmerI = kmer.length - len;
            if ( seqI < 0 || seq[seqI] != kmer[kmerI] ) {
                return len - 1;
            }
        }
        return kmer.length;
    }

    /**
     * Actually merge the dangling head if possible
     *
     * @param danglingHeadMergeResult   the result from generating a Cigar for the dangling head against the reference
     * @return 1 if merge was successful, 0 otherwise
     */
    @VisibleForTesting
    int mergeDanglingHead(final DanglingChainMergeHelper danglingHeadMergeResult) {

        final List<CigarElement> elements = danglingHeadMergeResult.cigar.getCigarElements();
        final CigarElement firstElement = elements.get(0);
        Utils.validateArg( firstElement.getOperator() == CigarOperator.M, "The first Cigar element must be an M");

        final int indexesToMerge = bestPrefixMatch(danglingHeadMergeResult.referencePathString, danglingHeadMergeResult.danglingPathString, firstElement.getLength());
        if ( indexesToMerge <= 0 ) {
            return 0;
        }

        // we can't push back the reference path
        if ( indexesToMerge >= danglingHeadMergeResult.referencePath.size() - 1 ) {
            return 0;
        }

        // but we can manipulate the dangling path if we need to
        if ( indexesToMerge >= danglingHeadMergeResult.danglingPath.size() &&
                ! extendDanglingPathAgainstReference(danglingHeadMergeResult, indexesToMerge - danglingHeadMergeResult.danglingPath.size() + 2) ) {
            return 0;
        }

        addEdge(danglingHeadMergeResult.referencePath.get(indexesToMerge+1), danglingHeadMergeResult.danglingPath.get(indexesToMerge), ((MyEdgeFactory)getEdgeFactory()).createEdge(false, 1));

        return 1;
    }

    /**
     * Generates the CIGAR string from the Smith-Waterman alignment of the dangling path (where the
     * provided vertex is the sink) and the reference path.
     *
     * @param vertex   the sink of the dangling chain
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @return a SmithWaterman object which can be null if no proper alignment could be generated
     */
    @VisibleForTesting
    final DanglingChainMergeHelper generateCigarAgainstDownwardsReferencePath(final MultiDeBruijnVertex vertex, final int pruneFactor, final int minDanglingBranchLength) {
        final int minTailPathLength = Math.max(1, minDanglingBranchLength); // while heads can be 0, tails absolutely cannot

        // find the lowest common ancestor path between this vertex and the diverging master path if available
        final List<MultiDeBruijnVertex> altPath = findPathUpwardsToLowestCommonAncestor(vertex, pruneFactor);
        if ( altPath == null || isRefSource(altPath.get(0)) || altPath.size() < minTailPathLength + 1 ) // add 1 to include the LCA
        {
            return null;
        }

        // now get the reference path from the LCA
        final List<MultiDeBruijnVertex> refPath = getReferencePath(altPath.get(0), TraversalDirection.downwards, Optional.ofNullable(incomingEdgeOf(altPath.get(1))));

        // create the Smith-Waterman strings to use
        final byte[] refBases = getBasesForPath(refPath, false);
        final byte[] altBases = getBasesForPath(altPath, false);

        // run Smith-Waterman to determine the best alignment (and remove trailing deletions since they aren't interesting)
        final SWPairwiseAlignment alignment = new SWPairwiseAlignment(refBases, altBases, SWPairwiseAlignment.STANDARD_NGS, SWPairwiseAlignment.OverhangStrategy.LEADING_INDEL);
        return new DanglingChainMergeHelper(altPath, refPath, altBases, refBases, AlignmentUtils.removeTrailingDeletions(alignment.getCigar()));
    }

    /**
     * Generates the CIGAR string from the Smith-Waterman alignment of the dangling path (where the
     * provided vertex is the source) and the reference path.
     *
     * @param vertex   the source of the dangling head
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @return a SmithWaterman object which can be null if no proper alignment could be generated
     */
    @VisibleForTesting
    final DanglingChainMergeHelper generateCigarAgainstUpwardsReferencePath(final MultiDeBruijnVertex vertex, final int pruneFactor, final int minDanglingBranchLength) {

        // find the highest common descendant path between vertex and the reference source if available
        final List<MultiDeBruijnVertex> altPath = findPathDownwardsToHighestCommonDescendantOfReference(vertex, pruneFactor);
        if ( altPath == null || isRefSink(altPath.get(0)) || altPath.size() < minDanglingBranchLength + 1 ) // add 1 to include the LCA
        {
            return null;
        }

        // now get the reference path from the LCA
        final List<MultiDeBruijnVertex> refPath = getReferencePath(altPath.get(0), TraversalDirection.upwards, Optional.empty());

        // create the Smith-Waterman strings to use
        final byte[] refBases = getBasesForPath(refPath, true);
        final byte[] altBases = getBasesForPath(altPath, true);

        // run Smith-Waterman to determine the best alignment (and remove trailing deletions since they aren't interesting)
        final SWPairwiseAlignment alignment = new SWPairwiseAlignment(refBases, altBases, SWPairwiseAlignment.STANDARD_NGS, SWPairwiseAlignment.OverhangStrategy.LEADING_INDEL);
        return new DanglingChainMergeHelper(altPath, refPath, altBases, refBases, AlignmentUtils.removeTrailingDeletions(alignment.getCigar()));
    }

    /**
     * Finds the path upwards in the graph from this vertex to the first diverging node, including that (lowest common ancestor) vertex.
     * Note that nodes are excluded if their pruning weight is less than the pruning factor.
     *
     * @param vertex   the original vertex
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @return the path if it can be determined or null if this vertex either doesn't merge onto another path or
     *  has an ancestor with multiple incoming edges before hitting the reference path
     */
    private List<MultiDeBruijnVertex> findPathUpwardsToLowestCommonAncestor(final MultiDeBruijnVertex vertex, final int pruneFactor) {
        return findPath(vertex, pruneFactor, v -> inDegreeOf(v) != 1 || outDegreeOf(v) >= 2, v -> outDegreeOf(v) > 1, v -> incomingEdgeOf(v), e -> getEdgeSource(e));
    }

    /**
     * Finds the path downwards in the graph from this vertex to the reference sequence, including the highest common descendant vertex.
     * However note that the path is reversed so that this vertex ends up at the end of the path.
     * Also note that nodes are excluded if their pruning weight is less than the pruning factor.
     *
     * @param vertex   the original vertex
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @return the path if it can be determined or null if this vertex either doesn't merge onto the reference path or
     *  has a descendant with multiple outgoing edges before hitting the reference path
     */
    private List<MultiDeBruijnVertex> findPathDownwardsToHighestCommonDescendantOfReference(final MultiDeBruijnVertex vertex, final int pruneFactor) {
        return findPath(vertex, pruneFactor, v -> isReferenceNode(v) || outDegreeOf(v) != 1, v -> isReferenceNode(v), v -> outgoingEdgeOf(v), e -> getEdgeTarget(e));
    }

    private  List<MultiDeBruijnVertex> findPath(final MultiDeBruijnVertex vertex, final int pruneFactor,
                                            final Predicate<MultiDeBruijnVertex> done,
                                            final Predicate<MultiDeBruijnVertex> returnPath,
                                            final Function<MultiDeBruijnVertex, MultiSampleEdge> nextEdge,
                                            final Function<MultiSampleEdge, MultiDeBruijnVertex> nextNode){
        final LinkedList<MultiDeBruijnVertex> path = new LinkedList<>();

        MultiDeBruijnVertex v = vertex;
        while ( ! done.test(v) ) {
            final MultiSampleEdge edge = nextEdge.apply(v);
            // if it has too low a weight, don't use it (or previous vertexes) for the path
            if ( edge.getPruningMultiplicity() < pruneFactor ) {
                path.clear();
            }// otherwise it is safe to use
            else {
                path.addFirst(v);
            }
            v = nextNode.apply(edge);
        }
        path.addFirst(v);

        return returnPath.test(v) ? path : null;
    }

    private enum TraversalDirection {
        downwards,
        upwards
    }

    /**
     * Finds the path in the graph from this vertex to the reference sink, including this vertex
     *
     * @param start   the reference vertex to start from
     * @param direction describes which direction to move in the graph (i.e. down to the reference sink or up to the source)
     * @param blacklistedEdges edges to ignore in the traversal down; useful to exclude the non-reference dangling paths
     * @return the path (non-null, non-empty)
     */
    private List<MultiDeBruijnVertex> getReferencePath(final MultiDeBruijnVertex start,
                                                       final TraversalDirection direction,
                                                       final Optional<MultiSampleEdge> blacklistedEdge) {

        final List<MultiDeBruijnVertex> path = new ArrayList<>();

        MultiDeBruijnVertex v = start;
        while ( v != null ) {
            path.add(v);
            v = (direction == TraversalDirection.downwards ? getNextReferenceVertex(v, true, blacklistedEdge) : getPrevReferenceVertex(v));
        }

        return path;
    }

    /**
     * The base sequence for the given path.
     *
     * @param path the list of vertexes that make up the path
     * @param expandSource if true and if we encounter a source node, then expand (and reverse) the character sequence for that node
     * @return  non-null sequence of bases corresponding to the given path
     */
    @VisibleForTesting
    byte[] getBasesForPath(final List<MultiDeBruijnVertex> path, final boolean expandSource) {
        Utils.nonNull(path, "Path cannot be null");

        final StringBuilder sb = new StringBuilder();
        for ( final MultiDeBruijnVertex v : path ) {
            if ( expandSource && isSource(v) ) {
                final String seq = v.getSequenceString();
                sb.append(new StringBuilder(seq).reverse().toString());
            } else {
                sb.append((char)v.getSuffix());
            }
        }

        return sb.toString().getBytes();
    }

    /**
     * Finds the index of the best extent of the prefix match between the provided paths, for dangling head merging.
     * Assumes that path1.length >= maxIndex and path2.length >= maxIndex.
     *
     * @param path1  the first path
     * @param path2  the second path
     * @param maxIndex the maximum index to traverse (not inclusive)
     * @return the index of the ideal prefix match or -1 if it cannot find one, must be less than maxIndex
     */
    private int bestPrefixMatch(final byte[] path1, final byte[] path2, final int maxIndex) {
        final int maxMismatches = getMaxMismatches(maxIndex);
        int mismatches = 0;
        int index = 0;
        int lastGoodIndex = -1;
        while ( index < maxIndex ) {
            if ( path1[index] != path2[index] ) {
                if ( ++mismatches > maxMismatches ) {
                    return -1;
                }
                lastGoodIndex = index;
            }
            index++;
        }
        // if we got here then we hit the max index
        return lastGoodIndex;
    }

    /**
     * Determine the maximum number of mismatches permitted on the branch.
     * Unless it's preset (e.g. by unit tests) it should be the length of the branch divided by the kmer size.
     *
     * @param lengthOfDanglingBranch  the length of the branch itself
     * @return positive integer
     */
    private int getMaxMismatches(final int lengthOfDanglingBranch) {
        return maxMismatchesInDanglingHead > 0 ? maxMismatchesInDanglingHead : Math.max(1, (lengthOfDanglingBranch / kmerSize));
    }

    private boolean extendDanglingPathAgainstReference(final DanglingChainMergeHelper danglingHeadMergeResult, final int numNodesToExtend) {

        final int indexOfLastDanglingNode = danglingHeadMergeResult.danglingPath.size() - 1;
        final int indexOfRefNodeToUse = indexOfLastDanglingNode + numNodesToExtend;
        if ( indexOfRefNodeToUse >= danglingHeadMergeResult.referencePath.size() ) {
            return false;
        }

        final MultiDeBruijnVertex danglingSource = danglingHeadMergeResult.danglingPath.remove(indexOfLastDanglingNode);
        final StringBuilder sb = new StringBuilder();
        final byte[] refSourceSequence = danglingHeadMergeResult.referencePath.get(indexOfRefNodeToUse).getSequence();
        for ( int i = 0; i < numNodesToExtend; i++ ) {
            sb.append((char) refSourceSequence[i]);
        }
        sb.append(danglingSource.getSequenceString());
        final byte[] sequenceToExtend = sb.toString().getBytes();

        // clean up the source and edge
        final MultiSampleEdge sourceEdge = outgoingEdgeOf(danglingSource);
        MultiDeBruijnVertex prevV = getEdgeTarget(sourceEdge);
        removeEdge(danglingSource, prevV);

        // extend the path
        for ( int i = numNodesToExtend; i > 0; i-- ) {
            final MultiDeBruijnVertex newV = new MultiDeBruijnVertex(Arrays.copyOfRange(sequenceToExtend, i, i + kmerSize));
            addVertex(newV);
            final MultiSampleEdge newE = addEdge(newV, prevV);
            newE.setMultiplicity(sourceEdge.getMultiplicity());
            danglingHeadMergeResult.danglingPath.add(newV);
            prevV = newV;
        }

        return true;
    }

    /**
     * Compute the smallest kmer size >= minKmerSize and <= maxKmerSize that has no non-unique kmers
     * among all sequences added to the current graph.  Will always return a result for maxKmerSize if
     * all smaller kmers had non-unique kmers.
     *
     * @param minKmerSize the minimum kmer size to consider when constructing the graph
     * @param maxKmerSize the maximum kmer size to consider
     * @return a non-null NonUniqueResult
     */
    private NonUniqueResult determineKmerSizeAndNonUniques(final int minKmerSize, final int maxKmerSize) {
        final Collection<SequenceForKmers> withNonUniques = getAllPendingSequences();
        final Set<Kmer> nonUniqueKmers = new HashSet<>();

        // go through the sequences and determine which kmers aren't unique within each read
        for (int kmerSize = minKmerSize ; kmerSize <= maxKmerSize; kmerSize++) {
            // clear out set of non-unique kmers
            nonUniqueKmers.clear();

            // loop over all sequences that have non-unique kmers in them from the previous iterator
            final Iterator<SequenceForKmers> it = withNonUniques.iterator();
            while ( it.hasNext() ) {
                final SequenceForKmers sequenceForKmers = it.next();

                // determine the non-unique kmers for this sequence
                final Collection<Kmer> nonUniquesFromSeq = determineNonUniqueKmers(sequenceForKmers, kmerSize);
                if ( nonUniquesFromSeq.isEmpty() ) {
                    // remove this sequence from future consideration
                    it.remove();
                } else {
                    // keep track of the non-uniques for this kmerSize, and keep it in the list of sequences that have non-uniques
                    nonUniqueKmers.addAll(nonUniquesFromSeq);
                }
            }

            if ( nonUniqueKmers.isEmpty() )
                // this kmerSize produces no non-unique sequences, so go ahead and use it for our assembly
            {
                break;
            }
        }

        // necessary because the loop breaks with kmerSize = max + 1
        return new NonUniqueResult(nonUniqueKmers);
    }

    /**
     * Get the collection of all sequences for kmers across all samples in no particular order
     * @return non-null Collection
     */
    private Collection<SequenceForKmers> getAllPendingSequences() {
        return pending.values().stream().flatMap(oneSampleWorth -> oneSampleWorth.stream()).collect(Collectors.toList());
    }

    /**
     * Get the collection of non-unique kmers from sequence for kmer size kmerSize
     * @param seqForKmers a sequence to get kmers from
     * @param kmerSize the size of the kmers
     * @return a non-null collection of non-unique kmers in sequence
     */
    static Collection<Kmer> determineNonUniqueKmers(final SequenceForKmers seqForKmers, final int kmerSize) {
        // count up occurrences of kmers within each read
        final Set<Kmer> allKmers = new LinkedHashSet<>();
        final List<Kmer> nonUniqueKmers = new ArrayList<>();
        final int stopPosition = seqForKmers.stop - kmerSize;
        for (int i = 0; i <= stopPosition; i++) {
            final Kmer kmer = new Kmer(seqForKmers.sequence, i, kmerSize);
            if (!allKmers.add(kmer)) {
                nonUniqueKmers.add(kmer);
            }
        }
        return nonUniqueKmers;
    }

    @Override
    public SeqGraph toSequenceGraph() {
        buildGraphIfNecessary();
        return super.toSequenceGraph();
    }

    private void increaseCountsInMatchedKmers(final SequenceForKmers seqForKmers,
                                              final MultiDeBruijnVertex vertex,
                                              final byte[] originalKmer,
                                              final int offset) {
        if ( offset == -1 ) {
            return;
        }

        for ( final MultiSampleEdge edge : incomingEdgesOf(vertex) ) {
            final MultiDeBruijnVertex prev = getEdgeSource(edge);
            final byte suffix = prev.getSuffix();
            final byte seqBase = originalKmer[offset];
            if ( suffix == seqBase && (increaseCountsThroughBranches || inDegreeOf(vertex) == 1) ) {
                edge.incMultiplicity(seqForKmers.count);
                increaseCountsInMatchedKmers(seqForKmers, prev, originalKmer, offset-1);
            }
        }
    }

    /**
     * Get the vertex for the kmer in sequence starting at start
     * @param sequence the sequence
     * @param start the position of the kmer start
     * @return a non-null vertex
     */
    private MultiDeBruijnVertex getOrCreateKmerVertex(final byte[] sequence, final int start) {
        final Kmer kmer = new Kmer(sequence, start, kmerSize);
        final MultiDeBruijnVertex vertex = getUniqueKmerVertex(kmer, true);
        return ( vertex != null ) ? vertex : createVertex(kmer);
    }

    /**
     * Get the unique vertex for kmer, or null if not possible.
     *
     * @param allowRefSource if true, we will allow kmer to match the reference source vertex
     * @return a vertex for kmer, or null if it's not unique
     */
    private MultiDeBruijnVertex getUniqueKmerVertex(final Kmer kmer, final boolean allowRefSource) {
        if ( ! allowRefSource && kmer.equals(refSource) ) {
            return null;
        }

        return uniqueKmers.get(kmer);
    }


    /**
     * Create a new vertex for kmer.  Add it to the uniqueKmers map if appropriate.
     *
     * kmer must not have a entry in unique kmers, or an error will be thrown
     *
     * @param kmer the kmer we want to create a vertex for
     * @return the non-null created vertex
     */
    private MultiDeBruijnVertex createVertex(final Kmer kmer) {
        final MultiDeBruijnVertex newVertex = new MultiDeBruijnVertex(kmer.bases());
        final int prevSize = vertexSet().size();
        addVertex(newVertex);

        // make sure we aren't adding duplicates (would be a bug)
        if ( vertexSet().size() != prevSize + 1) {
            throw new IllegalStateException("Adding vertex " + newVertex + " to graph didn't increase the graph size");
        }

        // add the vertex to the unique kmer map, if it is in fact unique
        if ( ! nonUniqueKmers.contains(kmer) && ! uniqueKmers.containsKey(kmer) ) // TODO -- not sure this last test is necessary
        {
            uniqueKmers.put(kmer, newVertex);
        }

        return newVertex;
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
    private MultiDeBruijnVertex extendChainByOne(final MultiDeBruijnVertex prevVertex, final byte[] sequence, final int kmerStart, final int count, final boolean isRef) {
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
        final MultiDeBruijnVertex uniqueMergeVertex = getUniqueKmerVertex(kmer, false);

        if ( isRef && uniqueMergeVertex != null ) {
            throw new IllegalStateException("Found a unique vertex to merge into the reference graph " + prevVertex + " -> " + uniqueMergeVertex);
        }

        // either use our unique merge vertex, or create a new one in the chain
        final MultiDeBruijnVertex nextVertex = uniqueMergeVertex == null ? createVertex(kmer) : uniqueMergeVertex;
        addEdge(prevVertex, nextVertex, ((MyEdgeFactory)getEdgeFactory()).createEdge(isRef, count));
        return nextVertex;
    }

    /**
     * Add the given read to the sequence graph.  Ultimately the read will get sent through addSequence(), but first
     * this method ensures we only use high quality bases and accounts for reduced reads, etc.
     *
     * @param read a non-null read
     */
    @VisibleForTesting
    void addRead(final GATKRead read, final SAMFileHeader header) {
        final byte[] sequence = read.getBases();
        final byte[] qualities = read.getBaseQualities();

        int lastGood = -1;
        for( int end = 0; end <= sequence.length; end++ ) {
            if ( end == sequence.length || ! baseIsUsableForAssembly(sequence[end], qualities[end]) ) {
                // the first good base is at lastGood, can be -1 if last base was bad
                final int start = lastGood;
                // the stop base is end - 1 (if we're not at the end of the sequence)
                final int len = end - start;

                if ( start != -1 && len >= kmerSize ) {
                    // if the sequence is long enough to get some value out of, add it to the graph
                    final String name = read.getName() + '_' + start + '_' + end;
                    addSequence(name, ReadUtils.getSampleName(read, header), read.getBases(), start, end, 1, false);
                }

                lastGood = -1; // reset the last good base
            } else if ( lastGood == -1 ) {
                lastGood = end; // we're at a good base, the last good one is us
            }
        }
    }

    /**
     * Determines whether a base can safely be used for assembly.
     * Currently disallows Ns and/or those with low quality
     *
     * @param base  the base under consideration
     * @param qual  the quality of that base
     * @return true if the base can be used for assembly, false otherwise
     */
    private boolean baseIsUsableForAssembly(final byte base, final byte qual) {
        return base != BaseUtils.Base.N.base && qual >= minBaseQualityToUseInAssembly;
    }

    /**
     * Get the set of non-unique kmers in this graph.  For debugging purposes
     * @return a non-null set of kmers
     */
    @VisibleForTesting
    Set<Kmer> getNonUniqueKmers() {
        return nonUniqueKmers;
    }

    @Override
    public String toString() {
        return "ReadThreadingAssembler{kmerSize=" + kmerSize + '}';
    }


    @Override
    public MultiDeBruijnVertex findKmer(final Kmer k) {
        return uniqueKmers.get(k);
    }


}