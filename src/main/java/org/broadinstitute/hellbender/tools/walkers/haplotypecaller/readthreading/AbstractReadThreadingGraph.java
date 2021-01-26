package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.Kmer;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.BaseGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.KmerSearchableGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignment;
import org.jgrapht.EdgeFactory;

import java.io.File;
import java.util.*;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

/**
 * Read threading graph class intended to contain duplicated code between {@link ReadThreadingGraph} and {@link JunctionTreeLinkedDeBruijnGraph}.
 */
public abstract class AbstractReadThreadingGraph extends BaseGraph<MultiDeBruijnVertex, MultiSampleEdge> implements KmerSearchableGraph<MultiDeBruijnVertex, MultiSampleEdge> {
    private static final long serialVersionUID = 1l;
    private static final String ANONYMOUS_SAMPLE = "XXX_UNNAMED_XXX";
    private static final boolean WRITE_GRAPH = false;
    private static final boolean DEBUG_NON_UNIQUE_CALC = false;
    private static final int MAX_CIGAR_COMPLEXITY = 3;
    private static final boolean INCREASE_COUNTS_BACKWARDS = true;
    private int minMatchingBasesToDangingEndRecovery = -1;

    /**
     * for debugging info printing
     */
    private static int counter = 0;
    /**
     * Sequences added for read threading before we've actually built the graph
     */
    protected final Map<String, List<SequenceForKmers>> pending = new LinkedHashMap<>();
    /**
     * A map from kmers -> their corresponding vertex in the graph
     */
    protected final Map<Kmer, MultiDeBruijnVertex> kmerToVertexMap = new LinkedHashMap<>();
    protected final boolean debugGraphTransformations;
    protected final byte minBaseQualityToUseInAssembly;
    protected List<MultiDeBruijnVertex> referencePath = null;
    protected boolean alreadyBuilt = false;
    // --------------------------------------------------------------------------------
    // state variables, initialized in setToInitialState()
    // --------------------------------------------------------------------------------
    private Kmer refSource = null;
    private boolean startThreadingOnlyAtExistingVertex = false;
    private int maxMismatchesInDanglingHead = -1; // this argument exists purely for testing purposes in constructing helpful tests and is currently not hooked up
    private boolean increaseCountsThroughBranches = false; // this may increase the branches without bounds

    protected enum TraversalDirection {
        downwards,
        upwards
    }

    AbstractReadThreadingGraph(int kmerSize, EdgeFactory<MultiDeBruijnVertex, MultiSampleEdge> edgeFactory) {
        super(kmerSize, edgeFactory);
        debugGraphTransformations = false;
        minBaseQualityToUseInAssembly = 0;
    }

    /**
     * Create a new ReadThreadingAssembler using kmerSize for matching
     *
     * @param kmerSize must be >= 1
     */
    public AbstractReadThreadingGraph(final int kmerSize, final boolean debugGraphTransformations, final byte minBaseQualityToUseInAssembly, final int numPruningSamples, final int numDanglingMatchingPrefixBases) {
        super(kmerSize, new MyEdgeFactory(numPruningSamples));

        Utils.validateArg(kmerSize > 0, () -> "bad minkKmerSize " + kmerSize);

        this.debugGraphTransformations = debugGraphTransformations;
        this.minBaseQualityToUseInAssembly = minBaseQualityToUseInAssembly;
        this.minMatchingBasesToDangingEndRecovery = numDanglingMatchingPrefixBases;
    }

    /**
     * Checks whether a kmer can be the threading start based on the current threading start location policy.
     *
     * @param kmer the query kmer.
     * @return {@code true} if we can start thread the sequence at this kmer, {@code false} otherwise.
     * @see #setThreadingStartOnlyAtExistingVertex(boolean)
     */
    protected abstract boolean isThreadingStart(final Kmer kmer, final boolean startThreadingOnlyAtExistingVertex);

    // get the next kmerVertex for ChainExtension and validate if necessary.
    protected abstract MultiDeBruijnVertex getNextKmerVertexForChainExtension(final Kmer kmer, final boolean isRef, final MultiDeBruijnVertex prevVertex);

    // perform any necessary preprocessing on the graph (such as non-unique kmer determination) before the graph is constructed
    protected abstract void preprocessReads();

    // heuristic to decide if the graph should be thown away based on low complexity/poor assembly
    public abstract boolean isLowQualityGraph();

    // whether reads are needed after graph construction
    protected abstract boolean shouldRemoveReadsAfterGraphConstruction();

    // Method that will be called immediately before haplotype finding in the event there are alteations that must be made to the graph based on implementation
    public abstract void postProcessForHaplotypeFinding(File debugGraphOutputPath, Locatable refHaplotype);

    /**
     * Define the behavior for how the graph should keep track of a potentially new kmer.
     *
     * @param kmer      (potentially) new kmer to track
     * @param newVertex corresponding vertex for that kmer
     */
    protected abstract void trackKmer(Kmer kmer, MultiDeBruijnVertex newVertex);

    /**
     * Determine whether the provided cigar is okay to merge into the reference path
     *
     * @param cigar                the cigar to analyze
     * @param requireFirstElementM if true, require that the first cigar element be an M operator in order for it to be okay
     * @param requireLastElementM  if true, require that the last cigar element be an M operator in order for it to be okay
     * @return true if it's okay to merge, false otherwise
     */
    @VisibleForTesting
    static boolean cigarIsOkayToMerge(final Cigar cigar, final boolean requireFirstElementM, final boolean requireLastElementM) {
        final List<CigarElement> elements = cigar.getCigarElements();
        final int numElements = elements.size();

        // don't allow more than a couple of different ops
        if (numElements == 0 || numElements > MAX_CIGAR_COMPLEXITY) {
            return false;
        }

        // the first element must be an M
        if (requireFirstElementM && elements.get(0).getOperator() != CigarOperator.M) {
            return false;
        }

        // the last element must be an M
        if (requireLastElementM && elements.get(numElements - 1).getOperator() != CigarOperator.M) {
            return false;
        }

        // note that there are checks for too many mismatches in the dangling branch later in the process

        return true;
    }

    /**
     * calculates the longest suffix match between a sequence and a smaller kmer
     *
     * @param seq      the (reference) sequence
     * @param kmer     the smaller kmer sequence
     * @param seqStart the index (inclusive) on seq to start looking backwards from
     * @return the longest matching suffix
     */
    @VisibleForTesting
    static int longestSuffixMatch(final byte[] seq, final byte[] kmer, final int seqStart) {
        for (int len = 1; len <= kmer.length; len++) {
            final int seqI = seqStart - len + 1;
            final int kmerI = kmer.length - len;
            if (seqI < 0 || seq[seqI] != kmer[kmerI]) {
                return len - 1;
            }
        }
        return kmer.length;
    }

    @VisibleForTesting
    protected void setAlreadyBuilt() {
        alreadyBuilt = true;
    }

    @VisibleForTesting
    void setMinMatchingBasesToDangingEndRecovery(final int minMatchingBasesToDangingEndRecovery) {
        this.minMatchingBasesToDangingEndRecovery = minMatchingBasesToDangingEndRecovery;
    }

    @VisibleForTesting
    void setMaxMismatchesInDanglingHead(final int maxMismatchesInDanglingHead) {
        this.maxMismatchesInDanglingHead = maxMismatchesInDanglingHead;
    }

    /**
     * Find vertex and its position in seqForKmers where we should start assembling seqForKmers
     *
     * @param seqForKmers the sequence we want to thread into the graph
     * @return the position of the starting vertex in seqForKmer, or -1 if it cannot find one
     */
    protected int findStart(final SequenceForKmers seqForKmers) {
        if (seqForKmers.isRef) {
            return 0;
        }

        for (int i = seqForKmers.start; i < seqForKmers.stop - kmerSize; i++) {
            final Kmer kmer1 = new Kmer(seqForKmers.sequence, i, kmerSize);
            if (isThreadingStart(kmer1, startThreadingOnlyAtExistingVertex)) {
                return i;
            }
        }

        return -1;
    }

    /**
     * Add the all bases in sequence to the graph
     *
     * @param sequence a non-null sequence
     * @param isRef    is this the reference sequence?
     */
    @VisibleForTesting
    public final void addSequence(final byte[] sequence, final boolean isRef) {
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
     * @param seqName  a useful seqName for this read, for debugging purposes
     * @param sequence non-null sequence of bases
     * @param start    the first base offset in sequence that we should use for constructing the graph using this sequence, inclusive
     * @param stop     the last base offset in sequence that we should use for constructing the graph using this sequence, exclusive
     * @param count    the representative count of this sequence (to use as the weight)
     * @param isRef    is this the reference sequence.
     */
    protected void addSequence(final String seqName, final String sampleName, final byte[] sequence, final int start, final int stop, final int count, final boolean isRef) {
        // note that argument testing is taken care of in SequenceForKmers
        Utils.validate(!alreadyBuilt, "Attempting to add sequence to a graph that has already been built");

        // get the list of sequences for this sample
        List<SequenceForKmers> sampleSequences = pending.computeIfAbsent(sampleName, s -> new LinkedList<>());

        // add the new sequence to the list of sequences for sample
        sampleSequences.add(new SequenceForKmers(seqName, sequence, start, stop, count, isRef));
    }

    /**
     * Thread sequence seqForKmers through the current graph, updating the graph as appropriate
     *
     * @param seqForKmers a non-null sequence
     */
    private void threadSequence(final SequenceForKmers seqForKmers) {
        final int startPos = findStart(seqForKmers);
        if (startPos == -1) {
            return;
        }

        final MultiDeBruijnVertex startingVertex = getOrCreateKmerVertex(seqForKmers.sequence, startPos);

        // increase the counts of all edges incoming into the starting vertex supported by going back in sequence
        if (INCREASE_COUNTS_BACKWARDS) {
            increaseCountsInMatchedKmers(seqForKmers, startingVertex, startingVertex.getSequence(), kmerSize - 2);
        }

        if (debugGraphTransformations) {
            startingVertex.addRead(seqForKmers.name);
        }

        // keep track of information about the reference source
        if (seqForKmers.isRef) {
            if (refSource != null) {
                throw new IllegalStateException("Found two refSources! prev: " + refSource + ", new: " + startingVertex);
            }
            referencePath = new ArrayList<>(seqForKmers.sequence.length - kmerSize);
            referencePath.add(startingVertex);
            refSource = new Kmer(seqForKmers.sequence, seqForKmers.start, kmerSize);
        }

        // loop over all of the bases in sequence, extending the graph by one base at each point, as appropriate
        MultiDeBruijnVertex vertex = startingVertex;
        for (int i = startPos + 1; i <= seqForKmers.stop - kmerSize; i++) {
            vertex = extendChainByOne(vertex, seqForKmers.sequence, i, seqForKmers.count, seqForKmers.isRef);
            if (seqForKmers.isRef) {
                referencePath.add(vertex);
            }
            if (debugGraphTransformations) {
                vertex.addRead(seqForKmers.name);
            }
        }
        if (seqForKmers.isRef) {
            referencePath = Collections.unmodifiableList(referencePath);
        }
    }

    /**
     * Changes the threading start location policy.
     *
     * @param value {@code true} if threading will start only at existing vertices in the graph, {@code false} if
     *              it can start at any unique kmer.
     */
    public final void setThreadingStartOnlyAtExistingVertex(final boolean value) {
        startThreadingOnlyAtExistingVertex = value;
    }

    /**
     * Build the read threaded assembly graph if it hasn't already been constructed from the sequences that have
     * been added to the graph.
     */
    public final void buildGraphIfNecessary() {
        if (alreadyBuilt) {
            return;
        }

        // Capture the set of non-unique kmers for the given kmer size (if applicable)
        preprocessReads();

        if (DEBUG_NON_UNIQUE_CALC) {
            ReadThreadingGraph.logger.info("using " + kmerSize + " kmer size for this assembly with the following non-uniques");
        }

        // go through the pending sequences, and add them to the graph
        for (final List<SequenceForKmers> sequencesForSample : pending.values()) {
            for (final SequenceForKmers sequenceForKmers : sequencesForSample) {
                threadSequence(sequenceForKmers);
                if (WRITE_GRAPH) {
                    printGraph(new File("threading." + counter++ + '.' + sequenceForKmers.name.replace(" ", "_") + ".dot"), 0);
                }
            }

            // flush the single sample edge values from the graph
            for (final MultiSampleEdge e : edgeSet()) {
                e.flushSingleSampleMultiplicity();
            }
        }

        // clear the pending reads pile to conserve memory
        if (shouldRemoveReadsAfterGraphConstruction()) {
            pending.clear();
        }
        alreadyBuilt = true;
        for (final MultiDeBruijnVertex v : kmerToVertexMap.values()) {
            v.setAdditionalInfo(v.getAdditionalInfo() + '+');
        }
    }

    @Override
    public boolean removeVertex(final MultiDeBruijnVertex V) {
        final boolean result = super.removeVertex(V);
        if (result) {
            final byte[] sequence = V.getSequence();
            final Kmer kmer = new Kmer(sequence);
            kmerToVertexMap.remove(kmer);
        }
        return result;
    }

    @Override
    public void removeSingletonOrphanVertices() {
        // Run through the graph and clean up singular orphaned nodes
        final Collection<MultiDeBruijnVertex> verticesToRemove = new LinkedList<>();
        for (final MultiDeBruijnVertex v : vertexSet()) {
            if (inDegreeOf(v) == 0 && outDegreeOf(v) == 0) {
                verticesToRemove.add(v);
            }
        }
        removeVertex(null);
        removeAllVertices(verticesToRemove);
    }

    @VisibleForTesting
    void setIncreaseCountsThroughBranches(final boolean increaseCountsThroughBranches) {
        this.increaseCountsThroughBranches = increaseCountsThroughBranches;
    }

    /**
     * Try to recover dangling tails
     *
     * @param pruneFactor             the prune factor to use in ignoring chain pieces if edge multiplicity is < pruneFactor
     * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
     * @param recoverAll              recover even branches with forks
     * @param aligner
     */
    public void recoverDanglingTails(final int pruneFactor, final int minDanglingBranchLength, final boolean recoverAll, final SmithWatermanAligner aligner) {
        Utils.validateArg(pruneFactor >= 0, () -> "pruneFactor must be non-negative but was " + pruneFactor);
        Utils.validateArg(minDanglingBranchLength >= 0, () -> "minDanglingBranchLength must be non-negative but was " + minDanglingBranchLength);

        if (!alreadyBuilt) {
            throw new IllegalStateException("recoverDanglingTails requires the graph be already built");
        }

        int attempted = 0;
        int nRecovered = 0;
        for (final MultiDeBruijnVertex v : vertexSet()) {
            if (outDegreeOf(v) == 0 && !isRefSink(v)) {
                attempted++;
                nRecovered += recoverDanglingTail(v, pruneFactor, minDanglingBranchLength, recoverAll, aligner);
            }
        }

        ReadThreadingGraph.logger.debug(String.format("Recovered %d of %d dangling tails", nRecovered, attempted));
    }

    /**
     * Try to recover dangling heads
     *
     * @param pruneFactor             the prune factor to use in ignoring chain pieces if edge multiplicity is < pruneFactor
     * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
     * @param recoverAll              recover even branches with forks
     * @param aligner
     */
    public void recoverDanglingHeads(final int pruneFactor, final int minDanglingBranchLength, final boolean recoverAll, final SmithWatermanAligner aligner) {
        Utils.validateArg(pruneFactor >= 0, () -> "pruneFactor must be non-negative but was " + pruneFactor);
        Utils.validateArg(minDanglingBranchLength >= 0, () -> "minDanglingBranchLength must be non-negative but was " + minDanglingBranchLength);
        if (!alreadyBuilt) {
            throw new IllegalStateException("recoverDanglingHeads requires the graph be already built");
        }

        // we need to build a list of dangling heads because that process can modify the graph (and otherwise generate
        // a ConcurrentModificationException if we do it while iterating over the vertexes)
        final Collection<MultiDeBruijnVertex> danglingHeads = vertexSet().stream()
                .filter(v -> inDegreeOf(v) == 0 && !isRefSource(v))
                .collect(Collectors.toList());

        // now we can try to recover the dangling heads
        int attempted = 0;
        int nRecovered = 0;
        for (final MultiDeBruijnVertex v : danglingHeads) {
            attempted++;
            nRecovered += recoverDanglingHead(v, pruneFactor, minDanglingBranchLength, recoverAll, aligner);
        }

        ReadThreadingGraph.logger.debug(String.format("Recovered %d of %d dangling heads", nRecovered, attempted));
    }

    /**
     * Attempt to attach vertex with out-degree == 0 to the graph
     *
     * @param vertex                  the vertex to recover
     * @param pruneFactor             the prune factor to use in ignoring chain pieces if edge multiplicity is < pruneFactor
     * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
     * @param aligner
     * @return 1 if we successfully recovered the vertex and 0 otherwise
     */
    private int recoverDanglingTail(final MultiDeBruijnVertex vertex, final int pruneFactor, final int minDanglingBranchLength, final boolean recoverAll, final SmithWatermanAligner aligner) {
        if (outDegreeOf(vertex) != 0) {
            throw new IllegalStateException("Attempting to recover a dangling tail for " + vertex + " but it has out-degree > 0");
        }

        // generate the CIGAR string from Smith-Waterman between the dangling tail and reference paths
        final DanglingChainMergeHelper danglingTailMergeResult = generateCigarAgainstDownwardsReferencePath(vertex, pruneFactor, minDanglingBranchLength, recoverAll, aligner);

        // if the CIGAR is too complex (or couldn't be computed) then we do not allow the merge into the reference path
        if (danglingTailMergeResult == null || !cigarIsOkayToMerge(danglingTailMergeResult.cigar, false, true)) {
            return 0;
        }

        // merge
        return mergeDanglingTail(danglingTailMergeResult);
    }

    /**
     * Finds the index of the best extent of the prefix match between the provided paths, for dangling head merging.
     * Requires that at a minimum there are at least #getMinMatchingBases() matches between the reference and the read
     * at the end in order to emit an alignment offset.
     *
     * @param cigarElements cigar elements corresponding to the alignment between path1 and path2
     * @param path1  the first path
     * @param path2  the second path
     * @return an integer pair object where the key is the offset into path1 and the value is offset into path2 (both -1 if no path is found)
     */
    private Pair<Integer, Integer> bestPrefixMatch(final List<CigarElement> cigarElements, final byte[] path1, final byte[] path2) {
        final int minMismatchingBases = getMinMatchingBases();

        int refIdx = cigarElements.stream().mapToInt(ce -> ce.getOperator().consumesReferenceBases()? ce.getLength() : 0).sum() - 1;
        int readIdx = path2.length - 1;

        // NOTE: this only works when the last cigar element has a sufficient number of M bases, so no indels within min-mismatches of the edge.
        for (final CigarElement ce : Lists.reverse(cigarElements)) {
            if (!(ce.getOperator().consumesReadBases() && ce.getOperator().consumesReferenceBases())) {
                break;
            }
            for (int j = 0; j < ce.getLength(); j++, refIdx--, readIdx--) {
                if (path1[refIdx] != path2[readIdx]) {
                    break;
                }
            }
        }

        final int matches = path2.length - 1 - readIdx;
        return matches < minMismatchingBases ? Pair.of(-1,-1) : Pair.of(readIdx, refIdx);
    }

    /**
     * The minimum number of matches to be considered allowable for recovering dangling ends
     */
    private int getMinMatchingBases() {
        return minMatchingBasesToDangingEndRecovery;
    }

    /**
     * Attempt to attach vertex with in-degree == 0, or a vertex on its path, to the graph
     *
     * @param vertex                  the vertex to recover
     * @param pruneFactor             the prune factor to use in ignoring chain pieces if edge multiplicity is < pruneFactor
     * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
     * @param recoverAll              recover even branches with forks
     * @param aligner
     * @return 1 if we successfully recovered a vertex and 0 otherwise
     */
    private int recoverDanglingHead(final MultiDeBruijnVertex vertex, final int pruneFactor, final int minDanglingBranchLength, final boolean recoverAll, final SmithWatermanAligner aligner) {
        if (inDegreeOf(vertex) != 0) {
            throw new IllegalStateException("Attempting to recover a dangling head for " + vertex + " but it has in-degree > 0");
        }

        // generate the CIGAR string from Smith-Waterman between the dangling tail and reference paths
        final DanglingChainMergeHelper danglingHeadMergeResult = generateCigarAgainstUpwardsReferencePath(vertex, pruneFactor, minDanglingBranchLength, recoverAll, aligner);

        // if the CIGAR is too complex (or couldn't be computed) then we do not allow the merge into the reference path
        if (danglingHeadMergeResult == null || !cigarIsOkayToMerge(danglingHeadMergeResult.cigar, true, false)) {
            return 0;
        }

        // merge
        return minMatchingBasesToDangingEndRecovery >= 0 ? mergeDanglingHead(danglingHeadMergeResult) : mergeDanglingHeadLegacy(danglingHeadMergeResult);
    }

    /**
     * Actually merge the dangling tail if possible
     *
     * @param danglingTailMergeResult the result from generating a Cigar for the dangling tail against the reference
     * @return 1 if merge was successful, 0 otherwise
     */
    @VisibleForTesting
    final int mergeDanglingTail(final DanglingChainMergeHelper danglingTailMergeResult) {

        final List<CigarElement> elements = danglingTailMergeResult.cigar.getCigarElements();
        final CigarElement lastElement = elements.get(elements.size() - 1);
        Utils.validateArg(lastElement.getOperator() == CigarOperator.M, "The last Cigar element must be an M");

        final int lastRefIndex = danglingTailMergeResult.cigar.getReferenceLength() - 1;
        final int matchingSuffix = Math.min(longestSuffixMatch(danglingTailMergeResult.referencePathString, danglingTailMergeResult.danglingPathString, lastRefIndex), lastElement.getLength());
        if (minMatchingBasesToDangingEndRecovery >= 0 ? matchingSuffix < minMatchingBasesToDangingEndRecovery : matchingSuffix == 0 ) {
            return 0;
        }

        final int altIndexToMerge = Math.max(danglingTailMergeResult.cigar.getReadLength() - matchingSuffix - 1, 0);

        // there is an important edge condition that we need to handle here: Smith-Waterman correctly calculates that there is a
        // deletion, that deletion is left-aligned such that the LCA node is part of that deletion, and the rest of the dangling
        // tail is a perfect match to the suffix of the reference path.  In this case we need to push the reference index to merge
        // down one position so that we don't incorrectly cut a base off of the deletion.
        final boolean firstElementIsDeletion = elements.get(0).getOperator() == CigarOperator.D;
        final boolean mustHandleLeadingDeletionCase = firstElementIsDeletion && (elements.get(0).getLength() + matchingSuffix == lastRefIndex + 1);
        final int refIndexToMerge = lastRefIndex - matchingSuffix + 1 + (mustHandleLeadingDeletionCase ? 1 : 0);

        // another edge condition occurs here: if Smith-Waterman places the whole tail into an insertion then it will try to
        // merge back to the LCA, which results in a cycle in the graph.  So we do not want to merge in such a case.
        if (refIndexToMerge == 0) {
            return 0;
        }

        // it's safe to merge now
        addEdge(danglingTailMergeResult.danglingPath.get(altIndexToMerge), danglingTailMergeResult.referencePath.get(refIndexToMerge), ((MyEdgeFactory) getEdgeFactory()).createEdge(false, 1));

        return 1;
    }

    /**
     * Actually merge the dangling head if possible, this is the old codepath that does not handle indels
     *
     * @param danglingHeadMergeResult   the result from generating a Cigar for the dangling head against the reference
     * @return 1 if merge was successful, 0 otherwise
     */
    @VisibleForTesting
    int mergeDanglingHeadLegacy(final DanglingChainMergeHelper danglingHeadMergeResult) {

        final List<CigarElement> elements = danglingHeadMergeResult.cigar.getCigarElements();
        final CigarElement firstElement = elements.get(0);
        Utils.validateArg(firstElement.getOperator() == CigarOperator.M, "The first Cigar element must be an M");

        final int indexesToMerge = bestPrefixMatchLegacy(danglingHeadMergeResult.referencePathString, danglingHeadMergeResult.danglingPathString, firstElement.getLength());
        if (indexesToMerge <= 0) {
            return 0;
        }

        // we can't push back the reference path
        if (indexesToMerge >= danglingHeadMergeResult.referencePath.size() - 1) {
            return 0;
        }

        // but we can manipulate the dangling path if we need to
        if (indexesToMerge >= danglingHeadMergeResult.danglingPath.size() &&
                !extendDanglingPathAgainstReference(danglingHeadMergeResult, indexesToMerge - danglingHeadMergeResult.danglingPath.size() + 2)) {
            return 0;
        }

        addEdge(danglingHeadMergeResult.referencePath.get(indexesToMerge + 1), danglingHeadMergeResult.danglingPath.get(indexesToMerge), ((MyEdgeFactory) getEdgeFactory()).createEdge(false, 1));

        return 1;
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

        final Pair<Integer, Integer> indexesToMerge = bestPrefixMatch(elements, danglingHeadMergeResult.referencePathString, danglingHeadMergeResult.danglingPathString);
        if ( indexesToMerge.getKey() <= 0 ||  indexesToMerge.getValue() <= 0 ) {
            return 0;
        }

        // we can't push back the reference path
        if ( indexesToMerge.getKey() >= danglingHeadMergeResult.referencePath.size() - 1 ) {
            return 0;
        }

        // but we can manipulate the dangling path if we need to
        if ( indexesToMerge.getValue() >= danglingHeadMergeResult.danglingPath.size() &&
                ! extendDanglingPathAgainstReference(danglingHeadMergeResult, indexesToMerge.getValue() - danglingHeadMergeResult.danglingPath.size() + 2) ) {
            return 0;
        }

        addEdge(danglingHeadMergeResult.referencePath.get(indexesToMerge.getKey()), danglingHeadMergeResult.danglingPath.get(indexesToMerge.getValue()), ((MyEdgeFactory)getEdgeFactory()).createEdge(false, 1));

        return 1;
    }

    /**
     * Generates the CIGAR string from the Smith-Waterman alignment of the dangling path (where the
     * provided vertex is the sink) and the reference path.
     *
     * @param aligner
     * @param vertex      the sink of the dangling chain
     * @param pruneFactor the prune factor to use in ignoring chain pieces if edge multiplicity is < pruneFactor
     * @param recoverAll  recover even branches with forks
     * @return a SmithWaterman object which can be null if no proper alignment could be generated
     */
    @VisibleForTesting
    DanglingChainMergeHelper generateCigarAgainstDownwardsReferencePath(final MultiDeBruijnVertex vertex, final int pruneFactor, final int minDanglingBranchLength, final boolean recoverAll, SmithWatermanAligner aligner) {
        final int minTailPathLength = Math.max(1, minDanglingBranchLength); // while heads can be 0, tails absolutely cannot

        // find the lowest common ancestor path between this vertex and the diverging master path if available
        final List<MultiDeBruijnVertex> altPath = findPathUpwardsToLowestCommonAncestor(vertex, pruneFactor, !recoverAll);
        if (altPath == null || isRefSource(altPath.get(0)) || altPath.size() < minTailPathLength + 1) // add 1 to include the LCA
        {
            return null;
        }

        // now get the reference path from the LCA
        final List<MultiDeBruijnVertex> refPath = getReferencePath(altPath.get(0), TraversalDirection.downwards, Optional.ofNullable(getHeaviestIncomingEdge(altPath.get(1))));

        // create the Smith-Waterman strings to use
        final byte[] refBases = getBasesForPath(refPath, false);
        final byte[] altBases = getBasesForPath(altPath, false);

        // run Smith-Waterman to determine the best alignment (and remove trailing deletions since they aren't interesting)
        final SmithWatermanAlignment alignment = aligner.align(refBases, altBases, SmithWatermanAligner.STANDARD_NGS, SWOverhangStrategy.LEADING_INDEL);
        return new DanglingChainMergeHelper(altPath, refPath, altBases, refBases, AlignmentUtils.removeTrailingDeletions(alignment.getCigar()));
    }

    /**
     * Generates the CIGAR string from the Smith-Waterman alignment of the dangling path (where the
     * provided vertex is the source) and the reference path.
     *
     * @param aligner
     * @param vertex      the source of the dangling head
     * @param pruneFactor the prune factor to use in ignoring chain pieces if edge multiplicity is < pruneFactor
     * @param recoverAll  recover even branches with forks
     * @return a SmithWaterman object which can be null if no proper alignment could be generated
     */
    @VisibleForTesting
    DanglingChainMergeHelper generateCigarAgainstUpwardsReferencePath(final MultiDeBruijnVertex vertex, final int pruneFactor, final int minDanglingBranchLength, final boolean recoverAll, SmithWatermanAligner aligner) {

        // find the highest common descendant path between vertex and the reference source if available
        final List<MultiDeBruijnVertex> altPath = findPathDownwardsToHighestCommonDescendantOfReference(vertex, pruneFactor, !recoverAll);
        if (altPath == null || isRefSink(altPath.get(0)) || altPath.size() < minDanglingBranchLength + 1) // add 1 to include the LCA
        {
            return null;
        }

        // now get the reference path from the LCA
        final List<MultiDeBruijnVertex> refPath = getReferencePath(altPath.get(0), TraversalDirection.upwards, Optional.empty());

        // create the Smith-Waterman strings to use
        final byte[] refBases = getBasesForPath(refPath, true);
        final byte[] altBases = getBasesForPath(altPath, true);

        // run Smith-Waterman to determine the best alignment (and remove trailing deletions since they aren't interesting)
        final SmithWatermanAlignment alignment = aligner.align(refBases, altBases, SmithWatermanAligner.STANDARD_NGS, SWOverhangStrategy.LEADING_INDEL);
        return new DanglingChainMergeHelper(altPath, refPath, altBases, refBases, AlignmentUtils.removeTrailingDeletions(alignment.getCigar()));
    }

    /**
     * Finds the path upwards in the graph from this vertex to the first diverging node, including that (lowest common ancestor) vertex.
     * Note that nodes are excluded if their pruning weight is less than the pruning factor.
     *
     * @param vertex         the original vertex
     * @param pruneFactor    the prune factor to use in ignoring chain pieces if edge multiplicity is < pruneFactor
     * @param giveUpAtBranch stop trying to find a path if a vertex with multiple incoming or outgoing edge is found
     * @return the path if it can be determined or null if this vertex either doesn't merge onto another path or
     * has an ancestor with multiple incoming edges before hitting the reference path
     */
    private List<MultiDeBruijnVertex> findPathUpwardsToLowestCommonAncestor(final MultiDeBruijnVertex vertex, final int pruneFactor, final boolean giveUpAtBranch) {
        return giveUpAtBranch ? findPath(vertex, pruneFactor, v -> inDegreeOf(v) != 1 || outDegreeOf(v) >= 2, v -> outDegreeOf(v) > 1, v -> incomingEdgeOf(v), e -> getEdgeSource(e))
                : findPath(vertex, pruneFactor, v -> hasIncidentRefEdge(v) || inDegreeOf(v) == 0, v -> outDegreeOf(v) > 1 && hasIncidentRefEdge(v), this::getHeaviestIncomingEdge, e -> getEdgeSource(e));
    }

    private boolean hasIncidentRefEdge(final MultiDeBruijnVertex v) {
        for (final MultiSampleEdge edge : incomingEdgesOf(v)) {
            if (edge.isRef()) {
                return true;
            }
        }
        return false;
    }

    private MultiSampleEdge getHeaviestIncomingEdge(final MultiDeBruijnVertex v) {
        final Set<MultiSampleEdge> incomingEdges = incomingEdgesOf(v);
        return incomingEdges.size() == 1 ? incomingEdges.iterator().next() :
                incomingEdges.stream().max(Comparator.comparingInt(MultiSampleEdge::getMultiplicity)).get();
    }

    /**
     * Finds the path downwards in the graph from this vertex to the reference sequence, including the highest common descendant vertex.
     * However note that the path is reversed so that this vertex ends up at the end of the path.
     * Also note that nodes are excluded if their pruning weight is less than the pruning factor.
     *
     * @param vertex         the original vertex
     * @param pruneFactor    the prune factor to use in ignoring chain pieces
     * @param giveUpAtBranch stop trying to find a path if a vertex with multiple incoming or outgoing edge is found
     * @return the path if it can be determined or null if this vertex either doesn't merge onto the reference path or
     * has a descendant with multiple outgoing edges before hitting the reference path
     */
    private List<MultiDeBruijnVertex> findPathDownwardsToHighestCommonDescendantOfReference(final MultiDeBruijnVertex vertex, final int pruneFactor, final boolean giveUpAtBranch) {
        return giveUpAtBranch ? findPath(vertex, pruneFactor, v -> isReferenceNode(v) || outDegreeOf(v) != 1, v -> isReferenceNode(v), v -> outgoingEdgeOf(v), e -> getEdgeTarget(e))
                : findPath(vertex, pruneFactor, v -> isReferenceNode(v) || outDegreeOf(v) == 0, v -> isReferenceNode(v), this::getHeaviestOutgoingEdge, e -> getEdgeTarget(e));
    }

    private MultiSampleEdge getHeaviestOutgoingEdge(final MultiDeBruijnVertex v) {
        final Set<MultiSampleEdge> outgoing = outgoingEdgesOf(v);
        return outgoing.size() == 1 ? outgoing.iterator().next() :
                outgoing.stream().max(Comparator.comparingInt(MultiSampleEdge::getMultiplicity)).get();
    }

    /**
     * Finds a path starting from a given vertex and satisfying various predicates
     *
     * @param vertex      the original vertex
     * @param pruneFactor the prune factor to use in ignoring chain pieces if edge multiplicity is < pruneFactor
     * @param done        test for whether a vertex is at the end of the path
     * @param returnPath  test for whether to return a found path based on its terminal vertex
     * @param nextEdge    function on vertices returning the next edge in the path
     * @param nextNode    function of edges returning the next vertex in the path
     * @return a path, if one satisfying all predicates is found, {@code null} otherwise
     */
    private List<MultiDeBruijnVertex> findPath(final MultiDeBruijnVertex vertex, final int pruneFactor,
                                               final Predicate<MultiDeBruijnVertex> done,
                                               final Predicate<MultiDeBruijnVertex> returnPath,
                                               final Function<MultiDeBruijnVertex, MultiSampleEdge> nextEdge,
                                               final Function<MultiSampleEdge, MultiDeBruijnVertex> nextNode) {
        final LinkedList<MultiDeBruijnVertex> path = new LinkedList<>();
        final Set<MultiDeBruijnVertex> visitedNodes = new HashSet<>(); // This code is necessary to

        MultiDeBruijnVertex v = vertex;
        while (!done.test(v)) {
            final MultiSampleEdge edge = nextEdge.apply(v);
            // if it has too low a weight, don't use it (or previous vertexes) for the path
            if (edge.getPruningMultiplicity() < pruneFactor) {
                // save the previously visited notes to protect us from riding forever 'neath the streets of boston
                visitedNodes.addAll(path);
                path.clear();
            } else {
                path.addFirst(v);
            }
            v = nextNode.apply(edge);
            // Check that we aren't stuck in a loop
            if (path.contains(v) || visitedNodes.contains(v)) {
                System.err.println("Dangling End recovery killed because of a loop (findPath)");
                return null;
            }
        }
        path.addFirst(v);

        return returnPath.test(v) ? path : null;
    }

    /**
     * Finds the path in the graph from this vertex to the reference sink, including this vertex
     *
     * @param start           the reference vertex to start from
     * @param direction       describes which direction to move in the graph (i.e. down to the reference sink or up to the source)
     * @param blacklistedEdge edge to ignore in the traversal down; useful to exclude the non-reference dangling paths
     * @return the path (non-null, non-empty)
     */
    protected List<MultiDeBruijnVertex> getReferencePath(final MultiDeBruijnVertex start,
                                                         final TraversalDirection direction,
                                                         final Optional<MultiSampleEdge> blacklistedEdge) {

        final List<MultiDeBruijnVertex> path = new ArrayList<>();

        MultiDeBruijnVertex v = start;
        while (v != null) {
            path.add(v);
            v = (direction == TraversalDirection.downwards ? getNextReferenceVertex(v, true, blacklistedEdge) : getPrevReferenceVertex(v));
        }

        return path;
    }

    /**
     * The base sequence for the given path.
     *
     * @param path         the list of vertexes that make up the path
     * @param expandSource if true and if we encounter a source node, then expand (and reverse) the character sequence for that node
     * @return non-null sequence of bases corresponding to the given path
     */
    @VisibleForTesting
    byte[] getBasesForPath(final List<MultiDeBruijnVertex> path, final boolean expandSource) {
        Utils.nonNull(path, "Path cannot be null");

        final StringBuilder sb = new StringBuilder();
        for (final MultiDeBruijnVertex v : path) {
            if (expandSource && isSource(v)) {
                final String seq = v.getSequenceString();
                sb.append(new StringBuilder(seq).reverse().toString());
            } else {
                sb.append((char) v.getSuffix());
            }
        }

        return sb.toString().getBytes();
    }

    /**
     * Finds the index of the best extent of the prefix match between the provided paths, for dangling head merging.
     * Assumes that path1.length >= maxIndex and path2.length >= maxIndex.
     *
     * @param path1    the first path
     * @param path2    the second path
     * @param maxIndex the maximum index to traverse (not inclusive)
     * @return the index of the ideal prefix match or -1 if it cannot find one, must be less than maxIndex
     */
    private int bestPrefixMatchLegacy(final byte[] path1, final byte[] path2, final int maxIndex) {
        final int maxMismatches = getMaxMismatchesLegacy(maxIndex);
        int mismatches = 0;
        int index = 0;
        int lastGoodIndex = -1;
        while (index < maxIndex) {
            if (path1[index] != path2[index]) {
                if (++mismatches > maxMismatches) {
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
     * NOTE: this method is only used for dangling heads and not tails.
     * 
     * Determine the maximum number of mismatches permitted on the branch.
     * Unless it's preset (e.g. by unit tests) it should be the length of the branch divided by the kmer size.
     *
     * @param lengthOfDanglingBranch the length of the branch itself
     * @return positive integer
     */
    private int getMaxMismatchesLegacy(final int lengthOfDanglingBranch) {
        return maxMismatchesInDanglingHead > 0 ? maxMismatchesInDanglingHead : Math.max(1, (lengthOfDanglingBranch / kmerSize));
    }

    private boolean extendDanglingPathAgainstReference(final DanglingChainMergeHelper danglingHeadMergeResult, final int numNodesToExtend) {

        final int indexOfLastDanglingNode = danglingHeadMergeResult.danglingPath.size() - 1;
        final int indexOfRefNodeToUse = indexOfLastDanglingNode + numNodesToExtend;
        if (indexOfRefNodeToUse >= danglingHeadMergeResult.referencePath.size()) {
            return false;
        }

        final MultiDeBruijnVertex danglingSource = danglingHeadMergeResult.danglingPath.remove(indexOfLastDanglingNode);
        final StringBuilder sb = new StringBuilder();
        final byte[] refSourceSequence = danglingHeadMergeResult.referencePath.get(indexOfRefNodeToUse).getSequence();
        for (int i = 0; i < numNodesToExtend; i++) {
            sb.append((char) refSourceSequence[i]);
        }
        sb.append(danglingSource.getSequenceString());
        final byte[] sequenceToExtend = sb.toString().getBytes();

        // clean up the source and edge
        final MultiSampleEdge sourceEdge = getHeaviestOutgoingEdge(danglingSource);
        MultiDeBruijnVertex prevV = getEdgeTarget(sourceEdge);
        removeEdge(danglingSource, prevV);

        // extend the path
        for (int i = numNodesToExtend; i > 0; i--) {
            final MultiDeBruijnVertex newV = new MultiDeBruijnVertex(Arrays.copyOfRange(sequenceToExtend, i, i + kmerSize));
            addVertex(newV);
            final MultiSampleEdge newE = addEdge(newV, prevV);
            newE.setMultiplicity(sourceEdge.getMultiplicity());
            danglingHeadMergeResult.danglingPath.add(newV);
            prevV = newV;
        }

        return true;
    }

    private void increaseCountsInMatchedKmers(final SequenceForKmers seqForKmers,
                                              final MultiDeBruijnVertex vertex,
                                              final byte[] originalKmer,
                                              final int offset) {
        if (offset == -1) {
            return;
        }

        for (final MultiSampleEdge edge : incomingEdgesOf(vertex)) {
            final MultiDeBruijnVertex prev = getEdgeSource(edge);
            final byte suffix = prev.getSuffix();
            final byte seqBase = originalKmer[offset];
            if (suffix == seqBase && (increaseCountsThroughBranches || inDegreeOf(vertex) == 1)) {
                edge.incMultiplicity(seqForKmers.count);
                increaseCountsInMatchedKmers(seqForKmers, prev, originalKmer, offset - 1);
            }
        }
    }

    /**
     * Get the vertex for the kmer in sequence starting at start
     *
     * @param sequence the sequence
     * @param start    the position of the kmer start
     * @return a non-null vertex
     */
    private MultiDeBruijnVertex getOrCreateKmerVertex(final byte[] sequence, final int start) {
        final Kmer kmer = new Kmer(sequence, start, kmerSize);
        final MultiDeBruijnVertex vertex = getKmerVertex(kmer, true);
        return (vertex != null) ? vertex : createVertex(kmer);
    }

    /**
     * Get the unique vertex for kmer, or null if not possible.
     *
     * @param allowRefSource if true, we will allow kmer to match the reference source vertex
     * @return a vertex for kmer, or null (either because it doesn't exist or is non-unique for graphs that have such a distinction)
     */
    protected MultiDeBruijnVertex getKmerVertex(final Kmer kmer, final boolean allowRefSource) {
        if (!allowRefSource && kmer.equals(refSource)) {
            return null;
        }

        return kmerToVertexMap.get(kmer);
    }

    /**
     * Create a new vertex for kmer.  Add it to the kmerToVertexMap map if appropriate.
     *
     * @param kmer the kmer we want to create a vertex for
     * @return the non-null created vertex
     */
    private MultiDeBruijnVertex createVertex(final Kmer kmer) {
        final MultiDeBruijnVertex newVertex = new MultiDeBruijnVertex(kmer.bases());
        final int prevSize = vertexSet().size();
        addVertex(newVertex);

        // make sure we aren't adding duplicates (would be a bug)
        if (vertexSet().size() != prevSize + 1) {
            throw new IllegalStateException("Adding vertex " + newVertex + " to graph didn't increase the graph size");
        }
        trackKmer(kmer, newVertex);

        return newVertex;
    }

    /**
     * Workhorse routine of the assembler.  Given a sequence whose last vertex is anchored in the graph, extend
     * the graph one bp according to the bases in sequence.
     *
     * @param prevVertex a non-null vertex where sequence was last anchored in the graph
     * @param sequence   the sequence we're threading through the graph
     * @param kmerStart  the start of the current kmer in graph we'd like to add
     * @param count      the number of observations of this kmer in graph (can be > 1 for GGA)
     * @param isRef      is this the reference sequence?
     * @return a non-null vertex connecting prevVertex to in the graph based on sequence
     */
    protected MultiDeBruijnVertex extendChainByOne(final MultiDeBruijnVertex prevVertex, final byte[] sequence, final int kmerStart, final int count, final boolean isRef) {
        final Set<MultiSampleEdge> outgoingEdges = outgoingEdgesOf(prevVertex);

        final int nextPos = kmerStart + kmerSize - 1;
        for (final MultiSampleEdge outgoingEdge : outgoingEdges) {
            final MultiDeBruijnVertex target = getEdgeTarget(outgoingEdge);
            if (target.getSuffix() == sequence[nextPos]) {
                // we've got a match in the chain, so simply increase the count of the edge by 1 and continue
                outgoingEdge.incMultiplicity(count);
                return target;
            }
        }

        // none of our outgoing edges had our unique suffix base, so we check for an opportunity to merge back in
        final Kmer kmer = new Kmer(sequence, kmerStart, kmerSize);
        final MultiDeBruijnVertex mergeVertex = getNextKmerVertexForChainExtension(kmer, isRef, prevVertex);

        // either use our merge vertex, or create a new one in the chain
        final MultiDeBruijnVertex nextVertex = mergeVertex == null ? createVertex(kmer) : mergeVertex;
        addEdge(prevVertex, nextVertex, ((MyEdgeFactory) getEdgeFactory()).createEdge(isRef, count));
        return nextVertex;
    }

    /**
     * Add a read to the sequence graph.  Finds maximal consecutive runs of bases with sufficient quality
     * and applies {@see addSequence} to these subreads if they are longer than the kmer size.
     *
     * @param read a non-null read
     */
    @VisibleForTesting
    void addRead(final GATKRead read, final SAMFileHeader header) {
        final byte[] sequence = read.getBases();
        final byte[] qualities = read.getBaseQualities();

        int lastGood = -1;
        for (int end = 0; end <= sequence.length; end++) {
            if (end == sequence.length || !baseIsUsableForAssembly(sequence[end], qualities[end])) {
                // the first good base is at lastGood, can be -1 if last base was bad
                final int start = lastGood;
                // the stop base is end - 1 (if we're not at the end of the sequence)
                final int len = end - start;

                if (start != -1 && len >= kmerSize) {
                    // if the sequence is long enough to get some value out of, add it to the graph
                    final String name = read.getName() + '_' + start + '_' + end;
                    addSequence(name, ReadUtils.getSampleName(read, header), sequence, start, end, 1, false);
                }

                lastGood = -1; // reset the last good base
            } else if (lastGood == -1) {
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
    protected boolean baseIsUsableForAssembly(final byte base, final byte qual) {
        return base != BaseUtils.Base.N.base && qual >= minBaseQualityToUseInAssembly;
    }

    @Override
    public MultiDeBruijnVertex findKmer(final Kmer k) {
        return kmerToVertexMap.get(k);
    }

    /**
     * Edge factory that encapsulates the numPruningSamples assembly parameter
     */
    protected static final class MyEdgeFactory implements EdgeFactory<MultiDeBruijnVertex, MultiSampleEdge> {
        final int numPruningSamples;

        MyEdgeFactory(final int numPruningSamples) {
            this.numPruningSamples = numPruningSamples;
        }

        @Override
        public MultiSampleEdge createEdge(final MultiDeBruijnVertex sourceVertex, final MultiDeBruijnVertex targetVertex) {
            return new MultiSampleEdge(false, 1, numPruningSamples);
        }

        MultiSampleEdge createEdge(final boolean isRef, final int multiplicity) {
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
            Utils.validateArg(start >= 0, () -> "Invalid start " + start);
            Utils.validateArg(stop >= start, () -> "Invalid stop " + stop);
            Utils.validateArg(count > 0, "Invalid count " + count);

            this.name = name;
            this.sequence = sequence;
            this.start = start;
            this.stop = stop;
            this.count = count;
            isRef = ref;
        }
    }
}
