package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public final class BaseGraphUnitTest extends BaseTest {
    SeqGraph graph;
    SeqVertex v1, v2, v3, v4, v5;

    @BeforeMethod
    public void setUp() {
        graph = new SeqGraph(11);

        v1 = new SeqVertex("A");
        v2 = new SeqVertex("C");
        v3 = new SeqVertex("C");
        v4 = new SeqVertex("C");
        v5 = new SeqVertex("C");

        graph.addVertices(v1, v2, v3, v4, v5);
        graph.addEdge(v1, v2);
        graph.addEdge(v2, v4);
        graph.addEdge(v3, v2);
        graph.addEdge(v2, v3);
        graph.addEdge(v4, v5);
    }

    @Test
    public void sequenceGraph(){
        final SeqGraph seqGraph = graph.toSequenceGraph();
        Assert.assertEquals(seqGraph.vertexSet().size(), 5);
        Assert.assertEquals(seqGraph.edgeSet().size(), 5);
    }
    @Test
    public void sourcesAndSinksAndNeighbors(){
        Assert.assertEquals(graph.getSources(), Collections.singleton(v1));
        Assert.assertEquals(graph.getSinks(), Collections.singleton(v5));
        Assert.assertEquals(graph.neighboringVerticesOf(v2), new HashSet<>(Arrays.asList(v1, v3, v4)));
        Assert.assertTrue(graph.hasCycles());
        Assert.assertTrue(graph.containsAllVertices(Arrays.asList(v1, v3, v4)));
        Assert.assertEquals(graph.subsetToNeighbors(v1, 0).vertexSet(), new HashSet<>(Arrays.asList(v1)));
        Assert.assertEquals(graph.subsetToNeighbors(v1, 1).vertexSet(), new HashSet<>(Arrays.asList(v1, v2)));
        Assert.assertEquals(graph.subsetToNeighbors(v1, 2).vertexSet(), new HashSet<>(Arrays.asList(v1, v2, v3, v4)));
        Assert.assertEquals(graph.subsetToNeighbors(v1, 3).vertexSet(), new HashSet<>(Arrays.asList(v1, v2, v3, v4, v5)));
        Assert.assertEquals(graph.subsetToNeighbors(v1, 3), graph);
        Assert.assertNull(graph.getReferenceSourceVertex());
    }

    @Test
    public void testIncomingAndOutgoingVertices() throws Exception {
        assertVertexSetEquals(graph.outgoingVerticesOf(v1), v2);
        assertVertexSetEquals(graph.incomingVerticesOf(v1));

        assertVertexSetEquals(graph.outgoingVerticesOf(v2), v3, v4);
        assertVertexSetEquals(graph.incomingVerticesOf(v2), v1, v3);

        assertVertexSetEquals(graph.outgoingVerticesOf(v3), v2);
        assertVertexSetEquals(graph.incomingVerticesOf(v3), v2);

        assertVertexSetEquals(graph.outgoingVerticesOf(v4), v5);
        assertVertexSetEquals(graph.incomingVerticesOf(v4), v2);

        assertVertexSetEquals(graph.outgoingVerticesOf(v5));
        assertVertexSetEquals(graph.incomingVerticesOf(v5), v4);
    }

    @Test
    public void testRemoveSingletonOrphanVertices() throws Exception {
        // all vertices in graph are connected
        final List<SeqVertex> kept = new LinkedList<>(graph.vertexSet());
        final SeqVertex rm1 = new SeqVertex("CAGT");
        final SeqVertex rm2 = new SeqVertex("AGTC");
        graph.addVertices(rm1, rm2);
        Assert.assertEquals(graph.vertexSet().size(), kept.size() + 2);
        final BaseEdge rm12e = new BaseEdge(false, 1);
        graph.addEdge(rm1, rm2, rm12e);

        final SeqGraph original = graph.clone();
        graph.removeSingletonOrphanVertices();
        Assert.assertTrue(BaseGraph.graphEquals(original, graph), "Graph with disconnected component but edges between components shouldn't be modified");

        graph.removeEdge(rm12e); // now we should be able to remove rm1 and rm2
        graph.removeSingletonOrphanVertices();
        Assert.assertTrue(graph.vertexSet().containsAll(kept));
        Assert.assertFalse(graph.containsVertex(rm1));
        Assert.assertFalse(graph.containsVertex(rm2));
    }

    @Test
    public void testRemoveSingletonOrphanVerticesOnSingleRefNode() throws Exception {
        final SeqGraph original = new SeqGraph(11);
        original.addVertex(v1);
        original.removeSingletonOrphanVertices();
        Assert.assertTrue(original.containsVertex(v1));
        Assert.assertEquals(original.vertexSet().size(), 1);
    }

    @Test
    public void testIsRefSourceAndSink() throws Exception {

        final SeqGraph g = new SeqGraph(11);
        Assert.assertEquals(g.getKmerSize(), 11);
        g.addVertex(v1);
        Assert.assertTrue(g.isRefSource(v1));
        Assert.assertTrue(g.isRefSink(v1));
        Assert.assertTrue(g.isReferenceNode(v1));

        g.addVertices(v2, v3, v4, v5);
        g.addEdge(v1, v2);
        g.addEdge(v2, v3);
        final BaseEdge refEdge = new BaseEdge(true, 1);

        final int edgeCountWayBefore = g.edgeSet().size();
        g.addOrUpdateEdge(v3, v4, refEdge);
        final int edgeCountBefore = g.edgeSet().size();
        Assert.assertEquals(edgeCountWayBefore+1, edgeCountBefore);
        g.addOrUpdateEdge(v3, v4, refEdge);
        final int edgeCountAfter = g.edgeSet().size();
        Assert.assertEquals(edgeCountBefore, edgeCountAfter);

        final BaseEdge v4v5 = g.addEdge(v4, v5);
        Assert.assertEquals(g.incomingEdgeOf(v5), v4v5);
        Assert.assertEquals(g.outgoingEdgeOf(v4), v4v5);

        Assert.assertFalse(g.isRefSource(v1));
        Assert.assertFalse(g.isRefSink(v1));
        Assert.assertFalse(g.isReferenceNode(v1));

        Assert.assertFalse(g.isRefSource(v2));
        Assert.assertFalse(g.isRefSink(v2));
        Assert.assertFalse(g.isReferenceNode(v2));

        Assert.assertTrue(g.isRefSource(v3));
        Assert.assertFalse(g.isRefSink(v3));
        Assert.assertTrue(g.isReferenceNode(v3));

        Assert.assertFalse(g.isRefSource(v4));
        Assert.assertTrue(g.isRefSink(v4));
        Assert.assertTrue(g.isReferenceNode(v4));

        Assert.assertFalse(g.isRefSource(v5));
        Assert.assertFalse(g.isRefSink(v5));
        Assert.assertFalse(g.isReferenceNode(v5));

        Assert.assertNull(g.getPrevReferenceVertex(null));
        Assert.assertNull(g.getPrevReferenceVertex(v1));
        Assert.assertEquals(g.getPrevReferenceVertex(v4), v3);

        Assert.assertNull(g.getNextReferenceVertex(null));
        Assert.assertNull(g.getNextReferenceVertex(v1));
        Assert.assertEquals(g.getNextReferenceVertex(v3), v4);

        Assert.assertEquals(g.getNextReferenceVertex(v4, false, Optional.of(v4v5)), null);
        Assert.assertEquals(g.getNextReferenceVertex(v4, true, Optional.of(v4v5)), null);

        Assert.assertEquals(g.subsetToRefSource(10), g);
        g.toString(); //just checking non-blowup

    }

    @Test
    public void testRemovePathsNotConnectedToRef() throws Exception {
        final SeqGraph graph = new SeqGraph(11);

        SeqVertex src = new SeqVertex("A");
        SeqVertex end = new SeqVertex("A");
        SeqVertex g1 = new SeqVertex("C");
        SeqVertex g2 = new SeqVertex("G");
        SeqVertex g3 = new SeqVertex("T");
        SeqVertex g4 = new SeqVertex("AA");
        SeqVertex g5 = new SeqVertex("AA");
        SeqVertex g6 = new SeqVertex("AA");
        SeqVertex g8 = new SeqVertex("AA");
        SeqVertex g7 = new SeqVertex("AA");
        SeqVertex b1 = new SeqVertex("CC");
        SeqVertex b2 = new SeqVertex("GG");
        SeqVertex b3 = new SeqVertex("TT");
        SeqVertex b4 = new SeqVertex("AAA");
        SeqVertex b5 = new SeqVertex("CCC");
        SeqVertex b6 = new SeqVertex("GGG");
        SeqVertex b7 = new SeqVertex("AAAA");
        SeqVertex b8 = new SeqVertex("GGGG");
        SeqVertex b9 = new SeqVertex("CCCC");

        graph.addVertices(src, end, g1, g2, g3, g4, g5, g6, g7, g8);
        graph.addEdges(() -> new BaseEdge(true, 1), src, g1, g2, g4, end);
        graph.addEdges(src, g1, g5, g6, g7, end);
        graph.addEdges(src, g1, g5, g8, g7, end);
        graph.addEdges(src, g1, g3, end);

        // the current state of the graph is the good one
        final SeqGraph good = graph.clone();

        // now add the bads to the graph
        graph.addVertices(b1, b2, b3, b4, b5, b6, b7, b8, b9);
        graph.addEdges(src, b1); // source -> b1 is dead
        graph.addEdges(b6, src); // x -> source is bad
        graph.addEdges(g4, b2); // off random vertex is bad
        graph.addEdges(g3, b3, b4); // two vertices that don't connect to end are bad
        graph.addEdges(end, b5); // vertex off end is bad
        graph.addEdges(g3, b7, b8, b7); // cycle is bad
        graph.addEdges(g3, b9, b9); // self-cycle is bad

        final boolean debug = false;
        if ( debug ) good.printGraph(new File("expected.dot"), 0);
        if ( debug ) graph.printGraph(new File("bad.dot"), 0);
        graph.removePathsNotConnectedToRef();
        if ( debug ) graph.printGraph(new File("actual.dot"), 0);

        Assert.assertTrue(BaseGraph.graphEquals(graph, good), "Failed to remove exactly the bad nodes");
    }

    @Test
    public void testRemoveVerticesNotConnectedToRefRegardlessOfEdgeDirection() throws Exception {
        final SeqGraph graph = new SeqGraph(11);

        SeqVertex src = new SeqVertex("A");
        SeqVertex end = new SeqVertex("A");
        SeqVertex g1 = new SeqVertex("C");
        SeqVertex g2 = new SeqVertex("G");
        SeqVertex g3 = new SeqVertex("T");
        SeqVertex g4 = new SeqVertex("AA");
        SeqVertex g5 = new SeqVertex("AA");
        SeqVertex g6 = new SeqVertex("AA");
        SeqVertex g8 = new SeqVertex("AA");
        SeqVertex g7 = new SeqVertex("AA");
        SeqVertex gPrev = new SeqVertex("AA");
        SeqVertex gPrev1 = new SeqVertex("AA");
        SeqVertex gPrev2 = new SeqVertex("AA");
        SeqVertex gAfter = new SeqVertex("AA");
        SeqVertex gAfter1 = new SeqVertex("AA");
        SeqVertex gAfter2 = new SeqVertex("AA");
        SeqVertex b1 = new SeqVertex("CC");
        SeqVertex b2 = new SeqVertex("GG");
        SeqVertex b3 = new SeqVertex("TT");
        SeqVertex b4 = new SeqVertex("AAA");
        SeqVertex b5 = new SeqVertex("CCC");
        SeqVertex b6 = new SeqVertex("GGG");

        graph.addVertices(src, end, g1, g2, g3, g4, g5, g6, g7, g8, gPrev, gPrev1, gPrev2, gAfter, gAfter1, gAfter2);
        graph.addEdges(() -> new BaseEdge(true, 1), src, g1, g2, g4, end);
        graph.addEdges(src, g1, g5, g6, g7, end);
        graph.addEdges(src, g1, g5, g8, g7, end);
        graph.addEdges(src, g1, g3, end);

        // these should be kept, but are in the wrong direction
        graph.addEdges(gPrev, src);
        graph.addEdges(gPrev1, gPrev2, src);
        graph.addEdges(end, gAfter);
        graph.addEdges(end, gAfter1, gAfter2);

        // the current state of the graph is the good one
        final SeqGraph good = graph.clone();

        // now add the bads to the graph
        graph.addVertices(b1, b2, b3, b4, b5, b6);
        graph.addEdges(b2, b3); // b2 -> b3
        graph.addEdges(b4, b5, b4); // cycle
        graph.addEdges(b6, b6); // isolated self cycle

        final boolean debug = false;
        if ( debug ) good.printGraph(new File("expected.dot"), 0);
        if ( debug ) graph.printGraph(new File("bad.dot"), 0);
        graph.removeVerticesNotConnectedToRefRegardlessOfEdgeDirection();
        if ( debug ) graph.printGraph(new File("actual.dot"), 0);

        Assert.assertTrue(BaseGraph.graphEquals(graph, good), "Failed to remove exactly the bad nodes");
    }

    @Test
    public void testPrintEmptyGraph() throws Exception {
        final File tmp = createTempFile("tmp", "dot");
        new SeqGraph(11).printGraph(tmp, 10);
        new TestGraph().printGraph(tmp, 10);
    }

    @Test
    public void testComplexGraph() throws Exception {
        final File tmp = createTempFile("tmp", "dot");
        graph.printGraph(tmp, 10);
    }

    private void assertVertexSetEquals(final Collection<SeqVertex> actual, final SeqVertex ... expected) {
        final Set<SeqVertex> actualSet = new HashSet<>(actual);
        Assert.assertEquals(actualSet.size(), actual.size(), "Duplicate elements found in vertex list");
        final Set<SeqVertex> expectedSet = expected == null ? Collections.emptySet() : new HashSet<>(Arrays.asList(expected));
        Assert.assertEquals(actualSet, expectedSet);
    }
}
