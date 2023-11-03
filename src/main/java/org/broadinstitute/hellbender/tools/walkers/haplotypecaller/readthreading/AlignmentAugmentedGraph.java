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
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
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

public class AlignmentAugmentedGraph extends BaseGraph<AlignmentAugmentedKmerVertex, MultiSampleEdge> {
    private static final long serialVersionUID = 1l;
    private static final String ANONYMOUS_SAMPLE = "XXX_UNNAMED_XXX";


    /**
     * Sequences added for read threading before we've actually built the graph
     */
    protected final Map<String, List<SequenceForKmers>> pending = new LinkedHashMap<>();

    // TODO: this should be from Kmers to a *list* of compatible vertices
    /**
     * A map from kmers -> their corresponding vertex in the graph
     */
    protected final Map<Kmer, MultiDeBruijnVertex> kmerToVertexMap = new LinkedHashMap<>();

    // TODO: need method for merging dangling ends



    protected final boolean debugGraphTransformations;
    protected final byte minBaseQualityToUseInAssembly;
    protected List<AlignmentAugmentedKmerVertex> referencePath = null;
    protected boolean alreadyBuilt = false;

    // --------------------------------------------------------------------------------
    // state variables, initialized in setToInitialState()
    // --------------------------------------------------------------------------------
    private Kmer refSource = null;
    private boolean startThreadingOnlyAtExistingVertex = false;
    private boolean increaseCountsThroughBranches = false; // this may increase the branches without bounds

    protected enum TraversalDirection {
        downwards,
        upwards
    }

    AlignmentAugmentedGraph(int kmerSize, EdgeFactory<MultiDeBruijnVertex, MultiSampleEdge> edgeFactory) {
        super(kmerSize, edgeFactory);
        debugGraphTransformations = false;
        minBaseQualityToUseInAssembly = 0;
    }

    /**
     * Create a new ReadThreadingAssembler using kmerSize for matching
     *
     * @param kmerSize must be >= 1
     */
    public AlignmentAugmentedGraph(final int kmerSize, final boolean debugGraphTransformations, final byte minBaseQualityToUseInAssembly, final int numPruningSamples, final int numDanglingMatchingPrefixBases) {
        super(kmerSize, new MyEdgeFactory(numPruningSamples));

        Utils.validateArg(kmerSize > 0, () -> "bad minkKmerSize " + kmerSize);

        this.debugGraphTransformations = debugGraphTransformations;
        this.minBaseQualityToUseInAssembly = minBaseQualityToUseInAssembly;
        this.minMatchingBasesToDanglingEndRecovery = numDanglingMatchingPrefixBases;
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
        increaseCountsInMatchedKmers(seqForKmers, startingVertex, startingVertex.getSequence(), kmerSize - 2);

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
     * Build the read threaded assembly graph if it hasn't already been constructed from the sequences that have
     * been added to the graph.
     */
    public final void buildGraphIfNecessary() {
        if (alreadyBuilt) {
            return;
        }

        // Capture the set of non-unique kmers for the given kmer size (if applicable)
        preprocessReads();

        // go through the pending sequences, and add them to the graph
        for (final List<SequenceForKmers> sequencesForSample : pending.values()) {
            for (final SequenceForKmers sequenceForKmers : sequencesForSample) {
                threadSequence(sequenceForKmers);
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

    /**
     * Try to recover dangling tails
     *
     * @param pruneFactor             the prune factor to use in ignoring chain pieces if edge multiplicity is < pruneFactor
     * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
     * @param recoverAll              recover even branches with forks
     */
    public void recoverDanglingTails(final int pruneFactor, final int minDanglingBranchLength, final boolean recoverAll, final SmithWatermanAligner aligner, final SWParameters danglingTailSWParameters) {
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
                nRecovered += recoverDanglingTail(v, pruneFactor, minDanglingBranchLength, recoverAll, aligner, danglingTailSWParameters);
            }
        }

        ReadThreadingGraph.logger.debug(String.format("Recovered %d of %d dangling tails", nRecovered, attempted));
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
            v = (direction == generateCigarAgainstDownwardsReferencePath()TraversalDirection.downwards ? getNextReferenceVertex(v, true, blacklistedEdge) : getPrevReferenceVertex(v));
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
