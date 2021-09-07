package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadThreadingAssemblerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.MultiDeBruijnVertex;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingGraph;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignmentConstants;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanJavaAligner;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public final class ChainPrunerUnitTest extends GATKBaseTest {
    private static final SWParameters DANGLING_END_SW_PARAMETERS = SmithWatermanAlignmentConstants.STANDARD_NGS;

    @DataProvider(name = "pruneLowWeightChainsData")
    public Object[][] makePruneChainsData() {
        List<Object[]> tests = new ArrayList<>();

        final SeqVertex v1 = new SeqVertex("A");
        final SeqVertex v2 = new SeqVertex("C");
        final SeqVertex v3 = new SeqVertex("G");
        final SeqVertex v4 = new SeqVertex("T");
        final SeqVertex v5 = new SeqVertex("AA");
        final SeqVertex v6 = new SeqVertex("CC");

        for ( final int edgeWeight : Arrays.asList(1, 2, 3) ) {
            for ( final int pruneFactor : Arrays.asList(1, 2, 3, 4) ) {
                for ( final boolean isRef : Arrays.asList(true, false)) {
                    // just an isolated chain
                    final int nExpected = edgeWeight < pruneFactor && ! isRef ? 3 : 0;
                    SeqGraph graph = new SeqGraph(11);
                    graph.addVertices(v1, v2, v3);
                    graph.addEdges(() -> new BaseEdge(isRef, edgeWeight), v1, v2, v3);
                    tests.add(new Object[]{"combinatorial", graph, pruneFactor, nExpected > 0 ? Collections.emptySet() : graph.vertexSet(), OptionalInt.of(1)});
                }
            }
        }

        { // connects to ref chain
            SeqGraph graph = new SeqGraph(11);
            graph.addVertices(v1, v2, v3);
            graph.addVertices(v4, v5);
            graph.addEdges(() -> new BaseEdge(true, 1), v4, v5);
            graph.addEdges(() -> new BaseEdge(false, 1), v4, v1, v2, v3, v5);
            tests.add(new Object[]{"bad internal branch", graph, 2, new HashSet<>(Arrays.asList(v4, v5)), OptionalInt.of(2)});
        }

        { // has bad cycle
            SeqGraph graph = new SeqGraph(11);
            graph.addVertices(v1, v2, v3, v4);
            graph.addEdges(() -> new BaseEdge(false, 1), v4, v1, v2, v3, v1);
            // note that we'll remove v4 because it's low weight
            tests.add(new Object[]{"has bad cycle", graph, 2, Collections.emptySet(), OptionalInt.of(2)});
        }

        { // has good cycle
            SeqGraph graph = new SeqGraph(111);
            graph.addVertices(v1, v2, v3, v4);
            graph.addEdges(() -> new BaseEdge(false, 3), v4, v1, v2, v3, v1);
            // note that we'll remove v4 because it's low weight
            tests.add(new Object[]{"has good cycle", graph, 2, graph.vertexSet(), OptionalInt.of(2)});
        }

        { // has branch
            SeqGraph graph = new SeqGraph(11);
            graph.addVertices(v1, v2, v3, v4, v5, v6);
            graph.addEdges(() -> new BaseEdge(false, 1), v1, v2, v3, v4, v6);
            graph.addEdges(() -> new BaseEdge(false, 1), v1, v2, v3, v5, v6);
            tests.add(new Object[]{"has two bad branches", graph, 2, Collections.emptySet(), OptionalInt.of(3)});
        }

        { // middle vertex above threshold => no one can be removed
            SeqGraph graph = new SeqGraph(11);
            graph.addVertices(v1, v2, v3, v4, v5);
            graph.addEdges(() -> new BaseEdge(false, 1), v1, v2);
            graph.addEdges(() -> new BaseEdge(false, 3), v2, v3);
            graph.addEdges(() -> new BaseEdge(false, 1), v3, v4, v5);
            tests.add(new Object[]{"middle vertex above factor", graph, 2, graph.vertexSet(), OptionalInt.of(1)});
        }

        { // the branching node has value > pruneFactor
            SeqGraph graph = new SeqGraph(11);
            graph.addVertices(v1, v2, v3, v4, v5, v6);
            graph.addEdges(() -> new BaseEdge(false, 3), v1, v2);
            graph.addEdges(() -> new BaseEdge(false, 3), v2, v3);
            graph.addEdges(() -> new BaseEdge(false, 1), v3, v4, v6);
            graph.addEdges(() -> new BaseEdge(false, 3), v2, v5, v6);
            tests.add(new Object[]{"branch node greater than pruneFactor", graph, 2, graph.vertexSet(), OptionalInt.of(3)});
        }

        { // A single isolated chain with weights all below pruning should be pruned
            SeqGraph graph = new SeqGraph(11);
            graph.addVertices(v1, v2, v3, v4, v5);
            graph.addEdges(() -> new BaseEdge(false, 1), v1, v2, v3);
            graph.addEdges(() -> new BaseEdge(false, 5), v4, v5);
            tests.add(new Object[]{"isolated chain", graph, 2, new LinkedHashSet<>(Arrays.asList(v4, v5)), OptionalInt.of(2)});
        }

        { // A chain with weights all below pruning should be pruned, even if it connects to another good chain
            SeqGraph graph = new SeqGraph(11);
            graph.addVertices(v1, v2, v3, v4, v5, v6);
            graph.addEdges(() -> new BaseEdge(false, 1), v1, v2, v3, v5);
            graph.addEdges(() -> new BaseEdge(false, 5), v4, v5, v6);
            tests.add(new Object[]{"bad chain branching into good one", graph, 2, new HashSet<>(Arrays.asList(v4, v5, v6)), OptionalInt.of(3)});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "pruneLowWeightChainsData")
    public void testPruneLowWeightChains(final String name, final SeqGraph graph, final int pruneFactor, final Set<SeqVertex> remainingVertices, final OptionalInt chainCountBeforePruning) {
        final Set<SeqVertex> copy = new HashSet<>(remainingVertices);
        final ChainPruner<SeqVertex, BaseEdge> pruner = new LowWeightChainPruner<>(pruneFactor);
        chainCountBeforePruning.ifPresent(n -> Assert.assertEquals(pruner.findAllChains(graph).size(), n));
        pruner.pruneLowWeightChains(graph);
        Assert.assertEquals(graph.vertexSet(), copy);
    }

    /**
     * Comprehensive test of adaptive pruning -- given an alt haplotype and a ref haplotype
     * @param kmerSize
     * @param ref           reference haplotype
     * @param alt           alt haplotype
     * @param altFraction   alt allele fraction to simulate somatic, mitochondrial etc variants
     * @param errorRate     substitution error rate of simulated reads
     * @param depthPerAlignmentStart    number of reads starting at each position.  Note that holding this constant yields
     *                                  low coverage at the beginning of the graph and high in the middle and end, simulating
     *                                  the leading edge of an exome target, for example
     * @param logOddsThreshold  log-10 odds threshold for pruning chains
     */
    @Test(dataProvider = "chainPrunerData")
    public void testAdaptivePruning(final int kmerSize, final byte[] ref, final byte[] alt, final double altFraction, final double errorRate, final int depthPerAlignmentStart, final double logOddsThreshold) {
        final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(kmerSize + FastMath.round(10000*(errorRate + altFraction))));
        final ReadThreadingGraph graph = new ReadThreadingGraph(kmerSize);
        graph.addSequence(ref, true);
        final List<byte[]> reads = IntStream.range(0, ref.length)
                .mapToObj(start -> IntStream.range(0, depthPerAlignmentStart).mapToObj(n -> generateReadWithErrors(rng.nextDouble() < altFraction ? alt : ref, start, errorRate, rng)))
                .flatMap(s -> s).collect(Collectors.toList());

        reads.forEach(read -> graph.addSequence(read, false));


        // note: these are the steps in ReadThreadingAssembler::createGraph
        graph.buildGraphIfNecessary();

        final ChainPruner<MultiDeBruijnVertex, MultiSampleEdge> pruner = new AdaptiveChainPruner<>(0.001,
                logOddsThreshold, ReadThreadingAssemblerArgumentCollection.DEFAULT_PRUNING_SEEDING_LOG_ODDS_THRESHOLD, 50);
        pruner.pruneLowWeightChains(graph);

        final SmithWatermanAligner aligner = SmithWatermanJavaAligner.getInstance();
        graph.recoverDanglingTails(1, 3, false, aligner, DANGLING_END_SW_PARAMETERS);
        graph.recoverDanglingHeads(1, 3, false, aligner, DANGLING_END_SW_PARAMETERS);
        graph.removePathsNotConnectedToRef();

        final SeqGraph seqGraph = graph.toSequenceGraph();
        seqGraph.zipLinearChains();
        seqGraph.removeSingletonOrphanVertices();
        seqGraph.removeVerticesNotConnectedToRefRegardlessOfEdgeDirection();
        seqGraph.simplifyGraph();
        seqGraph.removePathsNotConnectedToRef();
        seqGraph.simplifyGraph();

        final List<KBestHaplotype<SeqVertex, BaseEdge>> bestPaths = new GraphBasedKBestHaplotypeFinder<>(seqGraph).findBestHaplotypes(10);

        final OptionalInt altIndex = IntStream.range(0, bestPaths.size()).filter(n -> bestPaths.get(n).haplotype().basesMatch(alt)).findFirst();
        Assert.assertTrue(altIndex.isPresent());

        // ref path should not be pruned even if all reads are alt
        final OptionalInt refIndex = IntStream.range(0, bestPaths.size()).filter(n -> bestPaths.get(n).haplotype().basesMatch(ref)).findFirst();
        Assert.assertTrue(refIndex.isPresent());

        // the haplotype score is the sum of the log-10 of all branching fractions, so the alt haplotype score should come out to
        // around the log-10 of the allele fraction up to some fudge factor, assuming we didn't do any dumb pruning
        Assert.assertEquals(bestPaths.get(altIndex.getAsInt()).score(), Math.log10(altFraction), 0.5);
        Assert.assertTrue(bestPaths.size() < 15);
    }

    // test that in graph with good path A -> B -> C and bad edges A -> D -> C and D -> B that the adjacency of bad edges --
    // such that when bad edges meet the multiplicities do not indicate an error - does not harm pruning.
    // we test with and without a true variant path A -> E -> C
    @Test
    public void testAdaptivePruningWithAdjacentBadEdges() {
        final int goodMultiplicity = 1000;
        final int variantMultiplicity = 50;
        final int badMultiplicity = 5;

        final SeqVertex source = new SeqVertex("source");
        final SeqVertex sink = new SeqVertex("sink");
        final SeqVertex A = new SeqVertex("A");
        final SeqVertex B = new SeqVertex("B");
        final SeqVertex C = new SeqVertex("C");
        final SeqVertex D = new SeqVertex("D");
        final SeqVertex E = new SeqVertex("E");


        for (boolean variantPresent : new boolean[] {false, true}) {
            final SeqGraph graph = new SeqGraph(20);

            graph.addVertices(source, A, B, C, D, sink);
            graph.addEdges(() -> new BaseEdge(true, goodMultiplicity), source, A, B, C, sink);
            graph.addEdges(() -> new BaseEdge(false, badMultiplicity), A, D, C);
            graph.addEdges(() -> new BaseEdge(false, badMultiplicity), D, B);

            if (variantPresent) {
                graph.addVertices(E);
                graph.addEdges(() -> new BaseEdge(false, variantMultiplicity), A, E, C);
            }

            final ChainPruner<SeqVertex, BaseEdge> pruner = new AdaptiveChainPruner<>(0.01, 2.0,
                    ReadThreadingAssemblerArgumentCollection.DEFAULT_PRUNING_SEEDING_LOG_ODDS_THRESHOLD, 50);
            pruner.pruneLowWeightChains(graph);

            Assert.assertFalse(graph.containsVertex(D));
            if (variantPresent) {
                Assert.assertTrue(graph.containsVertex(E));
            }
        }
    }

    // test that in graph with good path A -> B -> C and bad edges A -> D and E -> C with a bubble with edges F, G between D and E
    // that the bad bubble does not harm pruning.
    // we test with and without a true variant path A -> H -> C
    @Test
    public void testAdaptivePruningWithBadBubble() {
        final int goodMultiplicity = 1000;
        final int variantMultiplicity = 50;
        final int badMultiplicity = 5;

        final SeqVertex source = new SeqVertex("source");
        final SeqVertex sink = new SeqVertex("sink");
        final SeqVertex A = new SeqVertex("A");
        final SeqVertex B = new SeqVertex("B");
        final SeqVertex C = new SeqVertex("C");
        final SeqVertex D = new SeqVertex("D");
        final SeqVertex E = new SeqVertex("E");
        final SeqVertex F = new SeqVertex("F");
        final SeqVertex G = new SeqVertex("G");
        final SeqVertex H = new SeqVertex("H");


        for (boolean variantPresent : new boolean[] {false, true}) {
            final SeqGraph graph = new SeqGraph(20);

            graph.addVertices(source, A, B, C, D, E, F, G, sink);
            graph.addEdges(() -> new BaseEdge(true, goodMultiplicity), source, A, B, C, sink);
            graph.addEdges(() -> new BaseEdge(false, badMultiplicity), A, D);
            graph.addEdges(() -> new BaseEdge(false, badMultiplicity), D, F, E);
            graph.addEdges(() -> new BaseEdge(false, badMultiplicity), D, G, E);
            graph.addEdges(() -> new BaseEdge(false, badMultiplicity), E, C);

            if (variantPresent) {
                graph.addVertices(H);
                graph.addEdges(() -> new BaseEdge(false, variantMultiplicity), A, H, C);
            }

            final ChainPruner<SeqVertex, BaseEdge> pruner = new AdaptiveChainPruner<>(0.01, ReadThreadingAssemblerArgumentCollection.DEFAULT_PRUNING_LOG_ODDS_THRESHOLD,
                    ReadThreadingAssemblerArgumentCollection.DEFAULT_PRUNING_SEEDING_LOG_ODDS_THRESHOLD, 50);
            pruner.pruneLowWeightChains(graph);

            Assert.assertFalse(graph.containsVertex(D));
            if (variantPresent) {
                Assert.assertTrue(graph.containsVertex(H));
            }
        }
    }

    @DataProvider(name = "chainPrunerData")
    public Object[][] getChainPrunerData() {
        final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(9));
        final int refLength = 100;
        final int leftSNVPosition = 15;
        final int middleSNVPosition = refLength / 2;
        final int rightSNVPosition = refLength - leftSNVPosition;

        final byte[] ref = new byte[refLength];
        IntStream.range(0, refLength).forEach(n -> ref[n] = BaseUtils.baseIndexToSimpleBase(rng.nextInt(4)));
        ref[leftSNVPosition] = 'A';
        ref[middleSNVPosition] = 'G';
        ref[rightSNVPosition] = 'T';

        final byte[] leftSNV = Arrays.copyOf(ref, refLength);
        leftSNV[leftSNVPosition] = 'G';

        final byte[] middleSNV = Arrays.copyOf(ref, refLength);
        middleSNV[middleSNVPosition] = 'T';

        final byte[] rightSNV = Arrays.copyOf(ref, refLength);
        rightSNV[rightSNVPosition] = 'A';

        // kmer size, ref bases, alt bases, alt fraction, base error rate, depth per start, log odds threshold, max unpruned variants
        return new Object[][] {
                { 10, ref, leftSNV, 0.5, 0.001, 20, ReadThreadingAssemblerArgumentCollection.DEFAULT_PRUNING_LOG_ODDS_THRESHOLD},
                { 10, ref, leftSNV, 1.0, 0.001, 20, ReadThreadingAssemblerArgumentCollection.DEFAULT_PRUNING_LOG_ODDS_THRESHOLD},
                { 10, ref, middleSNV, 0.1, 0.001, 5, ReadThreadingAssemblerArgumentCollection.DEFAULT_PRUNING_LOG_ODDS_THRESHOLD},
                { 25, ref, middleSNV, 0.1, 0.001, 5, ReadThreadingAssemblerArgumentCollection.DEFAULT_PRUNING_LOG_ODDS_THRESHOLD},
                { 25, ref, middleSNV, 0.01, 0.001, 1000, ReadThreadingAssemblerArgumentCollection.DEFAULT_PRUNING_LOG_ODDS_THRESHOLD},   // note the extreme depth -- this would confuse non-adaptive pruning
                { 10, ref, rightSNV, 0.1, 0.001, 2, ReadThreadingAssemblerArgumentCollection.DEFAULT_PRUNING_LOG_ODDS_THRESHOLD}
        };
    }

    private static byte[] generateReadWithErrors(final byte[] sequence, final int start, final int end, final double errorRate, final RandomGenerator rng) {
        final double adjustedErrorRate = errorRate * 4 / 3; // one time in four a random base won't be an error
        final byte[] result = new byte[end - start];
        IntStream.range(start, end).forEach(n -> {
            result[n - start] = rng.nextDouble() > adjustedErrorRate ? sequence[n] : BaseUtils.baseIndexToSimpleBase(rng.nextInt(4));
        });
        return result;
    }

    private static byte[] generateReadWithErrors(final byte[] sequence, final int start, final double errorRate, final RandomGenerator rng) {
        return generateReadWithErrors(sequence, start, sequence.length, errorRate, rng);
    }
}
