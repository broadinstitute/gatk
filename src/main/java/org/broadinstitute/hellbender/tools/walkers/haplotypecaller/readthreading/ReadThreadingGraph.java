package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
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
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignment;
import org.jgrapht.EdgeFactory;

import java.io.File;
import java.util.*;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

/**
 * Note: not final but only intended to be subclassed for testing.
 */
public class ReadThreadingGraph extends ReadThreadingGraphInterface {

    protected static final Logger logger = LogManager.getLogger(ReadThreadingGraph.class);

    private static final long serialVersionUID = 1l;

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
    }

    /**
     * Create a new ReadThreadingAssembler using kmerSize for matching
     * @param kmerSize must be >= 1
     */
    ReadThreadingGraph(final int kmerSize, final boolean debugGraphTransformations, final byte minBaseQualityToUseInAssembly, final int numPruningSamples) {
        super(kmerSize, debugGraphTransformations, minBaseQualityToUseInAssembly, numPruningSamples);
    }

    /**
     * Reset this assembler to its initial state, so we can create another assembly with a different set of reads
     */
    @Override
    protected void resetToInitialState() {
        pending.clear();
        nonUniqueKmers = null;
        uniqueKmers.clear();
        refSource = null;
        alreadyBuilt = false;
    }

    @Override
    protected void removePendingSequencesIfNecessary() {
        pending.clear();
    }

    /**
     * Does the graph not have enough complexity?  We define low complexity as a situation where the number
     * of non-unique kmers is more than 20% of the total number of kmers.
     *
     * @return true if the graph has low complexity, false otherwise
     */
    @Override
    public boolean isLowComplexity() {
        return nonUniqueKmers.size() * 4 > uniqueKmers.size();
    }

    @Override
    public ReadThreadingGraph clone() {
        return (ReadThreadingGraph) super.clone();
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
        final Collection<SequenceForKmers> withNonUniques = getAllPendingSequences();
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
     * Determines whether a base can safely be used for assembly.
     * Currently disallows Ns and/or those with low quality
     *
     * @param base  the base under consideration
     * @param qual  the quality of that base
     * @return true if the base can be used for assembly, false otherwise
     */
    @Override
    protected boolean baseIsUsableForAssembly(final byte base, final byte qual) {
        return base != BaseUtils.Base.N.base && qual >= minBaseQualityToUseInAssembly;
    }

    @Override
    public String toString() {
        return "ReadThreadingAssembler{kmerSize=" + kmerSize + '}';
    }


}