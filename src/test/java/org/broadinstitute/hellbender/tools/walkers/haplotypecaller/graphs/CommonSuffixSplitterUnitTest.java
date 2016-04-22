package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public final class CommonSuffixSplitterUnitTest extends BaseTest {
    private final static boolean DEBUG = false;

    @DataProvider(name = "SplitData")
    public Object[][] makeSplitData() {
        return CommonSuffixMergerUnitTest.makeSplitMergeData(-1);
    }

    @Test(dataProvider = "SplitData", enabled = !DEBUG)
    public void testSplit(final CommonSuffixMergerUnitTest.SplitMergeData data) {
        final boolean expectedMerge = ! data.commonSuffix.isEmpty() && data.graph.inDegreeOf(data.v) > 1;

        final SeqGraph original = data.graph.clone();
//        original.printGraph(new File("original.dot"), 0);
        final boolean succeed = CommonSuffixSplitter.split(data.graph, data.v);
//        data.graph.printGraph(new File("actual.dot"), 0);
        Assert.assertEquals(succeed, expectedMerge, "Not excepted merge success/fail result");
        if ( succeed ) {
            Assert.assertEquals(data.graph.incomingVerticesOf(data.v).iterator().next().getSequenceString(), data.commonSuffix, "Common suffix not computed correctly");
        }

        CommonSuffixMergerUnitTest.assertSameHaplotypes(String.format("suffixSplit.%s.%d", data.commonSuffix, data.graph.vertexSet().size()), data.graph, original);
    }

    @Test(enabled = !DEBUG)
    public void testSplitPrevHaveMultipleEdges() {
        final SeqGraph original = new SeqGraph(11);
        final SeqVertex v1 = new SeqVertex("A");
        final SeqVertex v2 = new SeqVertex("A");
        final SeqVertex v3 = new SeqVertex("A");
        final SeqVertex v4 = new SeqVertex("A");

        original.addVertices(v1, v2, v3, v4);
        original.addEdges(v1, v3);

        Assert.assertFalse(CommonSuffixSplitter.split(original, v3), "Cannot split graph with only one vertex");

        original.addEdges(v2, v3);
        original.addEdges(v2, v4);

        Assert.assertFalse(CommonSuffixSplitter.split(original, v3), "Cannot split graph with multiple outgoing edges from middle nodes");
    }

    @Test(enabled = !DEBUG)
    public void testSplitNoCycles() {
        final SeqGraph original = new SeqGraph(11);
        final SeqVertex v1 = new SeqVertex("A");
        final SeqVertex v2 = new SeqVertex("AC");
        final SeqVertex v3 = new SeqVertex("TC");
        final SeqVertex v4 = new SeqVertex("G");

        original.addVertices(v1, v2, v3, v4);
        original.addEdges(v1, v3, v4);
        original.addEdges(v1, v2, v4);

        Assert.assertTrue(CommonSuffixSplitter.split(original.clone(), v4), "Should be able to split pre-cycle graph");

        original.addEdges(v4, v4);
        Assert.assertFalse(CommonSuffixSplitter.split(original, v4), "Cannot split graph with a cycle of the bottom list");
    }

    @Test(timeOut = 10000, enabled = !DEBUG)
    public void testSplitComplexCycle() {
        final SeqGraph original = new SeqGraph(11);
        final SeqVertex r1 = new SeqVertex("ACTG");
        final SeqVertex r2 = new SeqVertex("ATGC");
        final SeqVertex cat1 = new SeqVertex("CAT");
        final SeqVertex cat2 = new SeqVertex("CAT");
        final SeqVertex c1 = new SeqVertex("C");
        final SeqVertex c2 = new SeqVertex("C");

        original.addVertices(r1, r2, cat1, cat2, c1, c2);
        original.addEdges(r1, cat1, c1, cat2, c1);
        original.addEdges(r2, c2, cat2);

        //original.printGraph(new File("testSplitComplexCycle.dot"), 0);

        for ( final SeqVertex v : Arrays.asList(cat2) ) { // original.vertexSet() ) {
            final SeqGraph graph = original.clone();
            final boolean success = CommonSuffixSplitter.split(graph, v);
            if ( success ) graph.printGraph(new File("testSplitComplexCycle.fail.dot"), 0);
            Assert.assertFalse(success, "Shouldn't be able to split any vertices but CommonSuffixSplitter says it could for " + v);
        }
    }

    @Test(timeOut = 10000)
    public void testSplitInfiniteCycleFailure() {
        final SeqGraph original = new SeqGraph(11);
        final SeqVertex v1 = new SeqVertex("GC");
        final SeqVertex v2 = new SeqVertex("X");
        final SeqVertex v3 = new SeqVertex("N");
        final SeqVertex v4 = new SeqVertex("C");

        original.addVertices(v1, v2, v3, v4);
        original.addEdge(v1, v2, new BaseEdge(false, 12));
        original.addEdge(v2, v3, new BaseEdge(false, 23));
        original.addEdge(v3, v4, new BaseEdge(false, 34));
        original.addEdge(v4, v2, new BaseEdge(false, 42));

//        original.printGraph(new File("testSplitInfiniteCycleFailure.dot"), 0);

        final SeqGraph graph = original.clone();
        final boolean success = CommonSuffixSplitter.split(graph, v2);
        Assert.assertTrue(success);

        for ( final SeqVertex v : graph.vertexSet() ) {
//            graph.printGraph(new File("testSplitInfiniteCycleFailure.first_split.dot"), 0);
            final boolean success2 = CommonSuffixSplitter.split(graph.clone(), v);
//            if ( success2 ) graph.printGraph(new File("testSplitInfiniteCycleFailure.fail.dot"), 0);
            Assert.assertFalse(success2, "Shouldn't be able to split any vertices but CommonSuffixSplitter says it could for " + v);
        }
    }
}
