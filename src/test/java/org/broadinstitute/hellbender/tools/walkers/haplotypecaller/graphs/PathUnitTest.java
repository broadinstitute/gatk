package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import htsjdk.samtools.Cigar;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.BaseEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.Path;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqVertex;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public final class PathUnitTest extends BaseTest {
    @Test
    public void testAlignReallyLongDeletion() {
        final String ref = "ATGGTGGCTCATACCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGAACATCACCTGAGGCCAGGAGTTCAAAACCAGCCTGGCTAACATAGCAAAACCCCATCTCTAATGAAAATACAAAAATTAGCTGGGTGTGGTGGTGTCCGCCTGTAGTCCCAGCTACTCAGGAGACTAAGGCATGAGAATCACTTGAACCCAGGATGCAGAGGCTGTAGTGAGCCGAGATTGCACCACGGCTGCACTCCAGCCTGGGCAACAGAGCGAGACTCTGTCTCAAATAAAATAGCGTAACGTAACATAACATAACATAACATAACATAACATAACATAACATAACATAACATAACATAACACAACAACAAAATAAAATAACATAAATCATGTTGTTAGGAAAAAAATCAGTTATGCAGCTACATGCTATTTACAAGAGATATACCTTAAAATATAAGACACAGAGGCCGGGCGCGGTAGCTCATGCCTGTAATCCCAGCACTTTGGGAGGCTGAGGCAAGCGGATCATGAGGTCAGGAGATCGAGACCATCC";
        final String hap = "ATGGTGGCTCATACCTGTAATCCCAGCACTTTGGGAGGCTGAGGCAAGCGGATCATGAGGTCAGGAGATCGAGACCATCCT";

        final SeqGraph graph = new SeqGraph(11);
        final SeqVertex v = new SeqVertex(hap);
        graph.addVertex(v);
        final Path<SeqVertex,BaseEdge> path = new Path<>(v, graph);
        final Cigar cigar = path.calculateCigar(ref.getBytes());
        Assert.assertNull(cigar, "Should have failed gracefully");
    }

    @Test
    public void testMakePath() {
        final SeqGraph g = new SeqGraph(3);
        final SeqVertex v1 = new SeqVertex("a");
        final SeqVertex v2 = new SeqVertex("b");
        final SeqVertex v3 = new SeqVertex("c");
        final SeqVertex v4 = new SeqVertex("d");
        g.addVertex(v1);   //source
        g.addVertex(v2);
        g.addVertex(v3);
        g.addVertex(v4);  //sink
        final BaseEdge e1 = g.addEdge(v1, v2);
        final BaseEdge e2 = g.addEdge(v2, v3);
        g.addEdge(v3, v4);
        final Path<SeqVertex,BaseEdge> path = new Path<>(v2, g);
        final Path<SeqVertex,BaseEdge> path1 = new Path<>(path, e2);
        final Path<SeqVertex,BaseEdge> path2 = new Path<>(e1, path1);
        Assert.assertEquals(path.length(), 0);
        Assert.assertEquals(path1.length(), 1);
        Assert.assertEquals(path2.length(), 2);

        Assert.assertTrue(path2.containsVertex(v1));
        Assert.assertFalse(path2.pathsAreTheSame(path1));
        Assert.assertTrue(path1.pathsAreTheSame(path1));
        path2.toString();//just test not blowing up - we dont' make any claims about the toString code

        Assert.assertEquals(path.getVertices().size(), 1);
        Assert.assertEquals(path.getLastVertex(), v2);
        Assert.assertEquals(path.getFirstVertex(), v2);
        Assert.assertEquals(path1.getFirstVertex(), v2);
        Assert.assertEquals(path1.getLastVertex(), v3);
        Assert.assertEquals(path.getBases(), "b".getBytes());
        Assert.assertEquals(path2.getBases(), "abc".getBytes());

    }

}
