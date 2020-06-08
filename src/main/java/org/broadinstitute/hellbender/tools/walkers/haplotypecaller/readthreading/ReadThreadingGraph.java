package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.Kmer;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqGraph;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.jgrapht.EdgeFactory;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Note: not final but only intended to be subclassed for testing.
 */
public class ReadThreadingGraph extends AbstractReadThreadingGraph {

    protected static final Logger logger = LogManager.getLogger(ReadThreadingGraph.class);

    private static final long serialVersionUID = 1l;

    /**
     * A set of non-unique kmers that cannot be used as merge points in the graph
     */
    protected Set<Kmer> nonUniqueKmers;

    /**
     * Constructs an empty read-threading-grpah provided the kmerSize.
     * @param kmerSize 1 or greater.
     *
     * @throws IllegalArgumentException if (@code kmerSize) < 1.
     */
    public ReadThreadingGraph(final int kmerSize) {
        this(kmerSize, false, (byte)6, 1, -1);
    }

    /**
     * Constructs a read-threading-graph for a string representation.
     *
     * <p>
     *     Note: only used for testing.
     * </p>
     */
    @VisibleForTesting
    protected ReadThreadingGraph(final int kmerSizeFromString, final EdgeFactory<MultiDeBruijnVertex, MultiSampleEdge> edgeFactory) {
        super(kmerSizeFromString, new MyEdgeFactory(1));
        nonUniqueKmers = null;
    }

    /**
     * Create a new ReadThreadingAssembler using kmerSize for matching
     * @param kmerSize must be >= 1
     */
    ReadThreadingGraph(final int kmerSize, final boolean debugGraphTransformations, final byte minBaseQualityToUseInAssembly, final int numPruningSamples, final int numDanglingMatchingPrefixBases) {
        super(kmerSize, debugGraphTransformations, minBaseQualityToUseInAssembly, numPruningSamples, numDanglingMatchingPrefixBases);
        nonUniqueKmers = null;
    }

    /**
     * Since we want to duplicate non-unique kmers in the graph code we must determine what those kmers are
     */
    @Override
    protected void preprocessReads() {
        nonUniqueKmers = determineNonUniques(kmerSize, getAllPendingSequences());
    }

    @Override
    protected boolean shouldRemoveReadsAfterGraphConstruction() {
        return true;
    }

    /**
     * Does the graph not have enough complexity?  We define low complexity as a situation where the number
     * of non-unique kmers is more than 20% of the total number of kmers.
     *
     * @return true if the graph has low complexity, false otherwise
     */
    @Override
    public boolean isLowQualityGraph() {
        return nonUniqueKmers.size() * 4 > kmerToVertexMap.size();
    }

    // only add the new kmer to the map if it exists and isn't in our non-unique kmer list
    @Override
    protected void trackKmer(final Kmer kmer, final MultiDeBruijnVertex newVertex) {
        if ( ! nonUniqueKmers.contains(kmer) && ! kmerToVertexMap.containsKey(kmer) ) {
            kmerToVertexMap.put(kmer, newVertex);
        }
    }

    @Override
    public ReadThreadingGraph clone() {
        return (ReadThreadingGraph) super.clone();
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
        return startThreadingOnlyAtExistingVertex ? kmerToVertexMap.containsKey(kmer) : !nonUniqueKmers.contains(kmer);
    }

    /**
     * Compute the smallest kmer size >= minKmerSize and <= maxKmerSize that has no non-unique kmers
     * among all sequences added to the current graph.  Will always return a result for maxKmerSize if
     * all smaller kmers had non-unique kmers.
     *
     * @param kmerSize the kmer size to check for non-unique kmers of
     * @return a non-null NonUniqueResult
     */
    private static Set<Kmer> determineNonUniques(final int kmerSize, Collection<SequenceForKmers> withNonUniques) {
        final Set<Kmer> nonUniqueKmers = new HashSet<>();

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

        return nonUniqueKmers;
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

    /**
     * Get the set of non-unique kmers in this graph.  For debugging purposes
     * @return a non-null set of kmers
     */
    @VisibleForTesting
    Set<Kmer> getNonUniqueKmers() {
        return nonUniqueKmers;
    }

    @Override
    protected MultiDeBruijnVertex getNextKmerVertexForChainExtension(final Kmer kmer, final boolean isRef, final MultiDeBruijnVertex prevVertex) {
        final MultiDeBruijnVertex uniqueMergeVertex = getKmerVertex(kmer, false);

        Utils.validate(!(isRef && uniqueMergeVertex != null), "Found a unique vertex to merge into the reference graph " + prevVertex + " -> " + uniqueMergeVertex);

        return uniqueMergeVertex;
    }

    @Override
    public void postProcessForHaplotypeFinding(final File debugGraphOutputPath, final Locatable refHaplotype) {
        return; // There is no processing required for this graph so simply return
    }

    @Override
    public String toString() {
        return "ReadThreadingAssembler{kmerSize=" + kmerSize + '}';
    }

}