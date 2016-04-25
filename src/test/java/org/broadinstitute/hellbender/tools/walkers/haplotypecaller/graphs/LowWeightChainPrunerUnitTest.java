package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public final class LowWeightChainPrunerUnitTest extends BaseTest {
    @DataProvider(name = "pruneChainsData")
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
                    { // just an isolated chain
                        final int nExpected = edgeWeight < pruneFactor && ! isRef ? 3 : 0;
                        SeqGraph graph = new SeqGraph(11);
                        graph.addVertices(v1, v2, v3);
                        graph.addEdges(() -> new BaseEdge(isRef, edgeWeight), v1, v2, v3);
                        tests.add(new Object[]{"combinatorial", graph, pruneFactor, nExpected > 0 ? Collections.emptySet() : graph.vertexSet()});
                    }
                }
            }
        }

        { // connects to ref chain
            SeqGraph graph = new SeqGraph(11);
            graph.addVertices(v1, v2, v3);
            graph.addVertices(v4, v5);
            graph.addEdges(() -> new BaseEdge(true, 1), v4, v5);
            graph.addEdges(() -> new BaseEdge(false, 1), v4, v1, v2, v3, v5);
            tests.add(new Object[]{"bad internal branch", graph, 2, new HashSet<>(Arrays.asList(v4, v5))});
        }

        { // has bad cycle
            SeqGraph graph = new SeqGraph(11);
            graph.addVertices(v1, v2, v3, v4);
            graph.addEdges(() -> new BaseEdge(false, 1), v4, v1, v2, v3, v1);
            // note that we'll remove v4 because it's low weight
            tests.add(new Object[]{"has bad cycle", graph, 2, Collections.emptySet()});
        }

        { // has good cycle
            SeqGraph graph = new SeqGraph(111);
            graph.addVertices(v1, v2, v3, v4);
            graph.addEdges(() -> new BaseEdge(false, 3), v4, v1, v2, v3, v1);
            // note that we'll remove v4 because it's low weight
            tests.add(new Object[]{"has good cycle", graph, 2, graph.vertexSet()});
        }

        { // has branch
            SeqGraph graph = new SeqGraph(11);
            graph.addVertices(v1, v2, v3, v4, v5, v6);
            graph.addEdges(() -> new BaseEdge(false, 1), v1, v2, v3, v4, v6);
            graph.addEdges(() -> new BaseEdge(false, 1), v1, v2, v3, v5, v6);
            tests.add(new Object[]{"has two bad branches", graph, 2, Collections.emptySet()});
        }

        { // middle vertex above threshold => no one can be removed
            SeqGraph graph = new SeqGraph(11);
            graph.addVertices(v1, v2, v3, v4, v5);
            graph.addEdges(() -> new BaseEdge(false, 1), v1, v2);
            graph.addEdges(() -> new BaseEdge(false, 3), v2, v3);
            graph.addEdges(() -> new BaseEdge(false, 1), v3, v4, v5);
            tests.add(new Object[]{"middle vertex above factor", graph, 2, graph.vertexSet()});
        }

        { // the branching node has value > pruneFactor
            SeqGraph graph = new SeqGraph(11);
            graph.addVertices(v1, v2, v3, v4, v5, v6);
            graph.addEdges(() -> new BaseEdge(false, 3), v1, v2);
            graph.addEdges(() -> new BaseEdge(false, 3), v2, v3);
            graph.addEdges(() -> new BaseEdge(false, 1), v3, v4, v6);
            graph.addEdges(() -> new BaseEdge(false, 3), v2, v5, v6);
            tests.add(new Object[]{"branch node greater than pruneFactor", graph, 2, graph.vertexSet()});
        }

        { // A single isolated chain with weights all below pruning should be pruned
            SeqGraph graph = new SeqGraph(11);
            graph.addVertices(v1, v2, v3, v4, v5);
            graph.addEdges(() -> new BaseEdge(false, 1), v1, v2, v3);
            graph.addEdges(() -> new BaseEdge(false, 5), v4, v5);
            tests.add(new Object[]{"isolated chain", graph, 2, new LinkedHashSet<>(Arrays.asList(v4, v5))});
        }

        { // A chain with weights all below pruning should be pruned, even if it connects to another good chain
            SeqGraph graph = new SeqGraph(11);
            graph.addVertices(v1, v2, v3, v4, v5, v6);
            graph.addEdges(() -> new BaseEdge(false, 1), v1, v2, v3, v5);
            graph.addEdges(() -> new BaseEdge(false, 5), v4, v5, v6);
            tests.add(new Object[]{"bad chain branching into good one", graph, 2, new HashSet<>(Arrays.asList(v4, v5, v6))});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "pruneChainsData", enabled = true)
    public void testPruneChains(final String name, final SeqGraph graph, final int pruneFactor, final Set<SeqVertex> remainingVertices) {
        final Set<SeqVertex> copy = new HashSet<>(remainingVertices);
//        graph.printGraph(new File("in.dot"), 0);
        final LowWeightChainPruner<SeqVertex, BaseEdge> pruner = new LowWeightChainPruner<>(pruneFactor);
        pruner.pruneLowWeightChains(graph);
//        graph.printGraph(new File("out.dot"), 0);
        Assert.assertEquals(graph.vertexSet(), copy);
    }
}
