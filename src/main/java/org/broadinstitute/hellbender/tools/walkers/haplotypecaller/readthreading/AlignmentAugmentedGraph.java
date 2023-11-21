package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import com.google.api.client.util.Lists;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.Kmer;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.AdaptiveChainPruner;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.ChainPruner;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;

import java.io.File;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public class AlignmentAugmentedGraph extends AbstractReadThreadingGraph {

    // TODO: magic constant, make adjustable
    private static final int MAX_BRANCH_VERTICES = 100;

    public AlignmentAugmentedGraph(final int kmerSize, final boolean debugGraphTransformations, final byte minBaseQualityToUseInAssembly) {
        super(kmerSize, debugGraphTransformations, minBaseQualityToUseInAssembly, 1, -1);
    }

    public void doAssembly(final Haplotype refHaplotype, final Iterable<GATKRead> reads, final ChainPruner pruner) {
        final Set<Kmer> decisionKmers = findDecisionKmers(refHaplotype, reads, pruner);

    }

    private Set<Kmer> findDecisionKmers(Haplotype refHaplotype, Iterable<GATKRead> reads, ChainPruner pruner) {
        // make a regular read threading graph in order to identify decision vertices
        final PlainDeBruijnGraph initialGraph = new PlainDeBruijnGraph(kmerSize, minBaseQualityToUseInAssembly);
        initialGraph.addSequence("ref", refHaplotype.getBases(), 1, true);
        reads.forEach(read -> initialGraph.addRead(read, null));
        initialGraph.buildGraphIfNecessary();
        pruner.pruneLowWeightChains(initialGraph);

        final List<MultiDeBruijnVertex> branchVertices = initialGraph.vertexSet().stream()
                .filter(v -> initialGraph.outDegreeOf(v) > 1).toList();

        final List<Pair<MultiDeBruijnVertex, Integer>> decisionVertices = new ArrayList<>();
        for (final MultiDeBruijnVertex vertex : branchVertices) {
            final List<MultiSampleEdge> edges = Lists.newArrayList(initialGraph.outgoingEdgesOf(vertex));
            final int[] multiplicities = edges.stream().mapToInt(MultiSampleEdge::getMultiplicity).toArray();
            final int branchiness = (int) MathUtils.sum(multiplicities) - MathUtils.arrayMax(multiplicities);
            // TODO: if there are more than 2 outgoing edges this gives minor edges more branchiness than they deserve
            edges.forEach(edge -> decisionVertices.add(Pair.of(initialGraph.getEdgeTarget(edge), branchiness)));
        }

        return decisionVertices.stream()
                .sorted(Comparator.comparingInt(Pair<MultiDeBruijnVertex, Integer>::getRight).reversed())
                .limit(MAX_BRANCH_VERTICES)
                .map(Pair::getLeft)
                .map(vertex -> new Kmer(vertex.getSequence()))
                .collect(Collectors.toSet());
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
    protected void preprocessReads() {
        throw new UnsupportedOperationException("HAven't implemented yet.");
    }

    // heuristic to decide if the graph should be thrown away based on low complexity/poor assembly
    @Override
    public boolean isLowQualityGraph() {
        return false;
    }

    // whether reads are needed after graph construction
    @Override
    protected boolean shouldRemoveReadsAfterGraphConstruction() {
        return false;
    }

    // Method that will be called immediately before haplotype finding in the event there are alteations that must be made to the graph based on implementation
    @Override
    public void postProcessForHaplotypeFinding(File debugGraphOutputPath, Locatable refHaplotype) {
        throw new UnsupportedOperationException("HAven't implemented yet.");
    }

    /**
     * Define the behavior for how the graph should keep track of a potentially new kmer.
     *
     * @param kmer      (potentially) new kmer to track
     * @param newVertex corresponding vertex for that kmer
     */
    protected abstract void trackKmer(Kmer kmer, MultiDeBruijnVertex newVertex);

    @Override
    public void recoverDanglingTails(final int pruneFactor, final int minDanglingBranchLength, final boolean recoverAll,
                                     final SmithWatermanAligner aligner, final SWParameters danglingTailSWParameters) {
        throw new UnsupportedOperationException("Not implemented yet but it's going to be different.");
    }

    @Override
    public void recoverDanglingHeads(final int pruneFactor, final int minDanglingBranchLength, final boolean recoverAll,
                                     final SmithWatermanAligner aligner, final SWParameters danglingTailSWParameters) {
        throw new UnsupportedOperationException("Not implemented yet but it's going to be different.");
    }



}
