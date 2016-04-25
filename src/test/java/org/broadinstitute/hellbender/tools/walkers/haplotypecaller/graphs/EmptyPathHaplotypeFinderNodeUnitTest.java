package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Collections;

public final class EmptyPathHaplotypeFinderNodeUnitTest {

    @Test
    public void testCount(){
        final SeqGraph g = new SeqGraph(5);
        final SeqVertex v = new SeqVertex("acgt");
        g.addVertex(v);
        final EmptyPathHaplotypeFinderNode f = new EmptyPathHaplotypeFinderNode(g, v);
        Assert.assertEquals(f.getCount(), 1);
    }

    @Test
    public void testsubFinderLabels(){
        final SeqGraph g = new SeqGraph(5);
        final SeqVertex v = new SeqVertex("acgt");
        g.addVertex(v);
        final EmptyPathHaplotypeFinderNode f = new EmptyPathHaplotypeFinderNode(g, v);
        Assert.assertEquals(f.subFinderLabels(), Collections.emptySet());
    }

    @Test
    public void testID(){
        final SeqGraph g = new SeqGraph(5);
        final SeqVertex v1 = new SeqVertex("acgt");
        final SeqVertex v2 = new SeqVertex("fred");
        g.addVertices(v1, v2);
        final EmptyPathHaplotypeFinderNode f1 = new EmptyPathHaplotypeFinderNode(g, v1);
        final EmptyPathHaplotypeFinderNode f2 = new EmptyPathHaplotypeFinderNode(g, v1);
        final EmptyPathHaplotypeFinderNode f3 = new EmptyPathHaplotypeFinderNode(g, v2);
        Assert.assertEquals(f1.id(), f2.id());
        Assert.assertNotEquals(f1.id(), f3.id());
        Assert.assertNotEquals(f2.id(), f3.id());
    }

    @Test
    public void testLabel(){
        final SeqGraph g = new SeqGraph(5);
        final SeqVertex v = new SeqVertex("acgt");
        g.addVertex(v);
        final EmptyPathHaplotypeFinderNode f = new EmptyPathHaplotypeFinderNode(g, v);
        Assert.assertEquals(f.label(), "acgt");
    }

    @Test
    public void testKzero(){
        final SeqGraph g = new SeqGraph(5);
        final SeqVertex v = new SeqVertex("acgt");
        g.addVertex(v);
        final EmptyPathHaplotypeFinderNode f = new EmptyPathHaplotypeFinderNode(g, v);
        final KBestHaplotype kBest = f.getKBest(0);
        Assert.assertEquals(kBest.bases(), "acgt".getBytes());
        Assert.assertEquals(kBest.graph(), g);
        Assert.assertEquals(kBest.score(), 0.0);
        Assert.assertEquals(kBest.rank(), 0);
        Assert.assertEquals(kBest.head(), v);
        Assert.assertEquals(kBest.tail(), null);
        Assert.assertEquals(kBest.isReference(), true);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testKneg(){
        final SeqGraph g = new SeqGraph(5);
        final SeqVertex v = new SeqVertex("acgt");
        final EmptyPathHaplotypeFinderNode f = new EmptyPathHaplotypeFinderNode(g, v);
        f.getKBest(-1);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testKpos(){
        final SeqGraph g = new SeqGraph(5);
        final SeqVertex v = new SeqVertex("acgt");
        final EmptyPathHaplotypeFinderNode f = new EmptyPathHaplotypeFinderNode(g, v);
        f.getKBest(1);
    }

    @Test
    public void testReference(){
        final SeqGraph g = new SeqGraph(5);
        final SeqVertex v = new SeqVertex("acgt");
        g.addVertex(v);
        final EmptyPathHaplotypeFinderNode f = new EmptyPathHaplotypeFinderNode(g, v);
        Assert.assertEquals(f.isReference(), true);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testReferenceError(){
        final SeqGraph g = new SeqGraph(5);
        final SeqVertex v = new SeqVertex("acgt");
        final EmptyPathHaplotypeFinderNode f = new EmptyPathHaplotypeFinderNode(g, v);
        f.isReference(); //blows up because v is not in g        //TODO: is this the intended semantics?
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testScoreNullBases(){
        final SeqGraph g = new SeqGraph(5);
        final SeqVertex v = new SeqVertex("acgt");
        g.addVertex(v);
        final EmptyPathHaplotypeFinderNode finder = new EmptyPathHaplotypeFinderNode(g, v);
        finder.score(null, 1, 1);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testScoreNegOffset(){
        final SeqGraph g = new SeqGraph(5);
        final SeqVertex v = new SeqVertex("acgt");
        g.addVertex(v);
        final EmptyPathHaplotypeFinderNode finder = new EmptyPathHaplotypeFinderNode(g, v);
        finder.score(new byte[]{'a', 'c', 'g', 't'}, -1, 1);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testScoreNegLength(){
        final SeqGraph g = new SeqGraph(5);
        final SeqVertex v = new SeqVertex("acgt");
        g.addVertex(v);
        final EmptyPathHaplotypeFinderNode finder = new EmptyPathHaplotypeFinderNode(g, v);
        finder.score(new byte[]{'a', 'c', 'g', 't'}, 1, -1);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testScoreEnd(){
        final SeqGraph g = new SeqGraph(5);
        final SeqVertex v = new SeqVertex("acgt");
        g.addVertex(v);
        final EmptyPathHaplotypeFinderNode finder = new EmptyPathHaplotypeFinderNode(g, v);
        finder.score(new byte[]{'a', 'c', 'g', 't'}, 2, 3);
    }

    @Test
    public void testScore(){
        final SeqGraph g = new SeqGraph(5);
        final SeqVertex v = new SeqVertex("acgt");
        g.addVertex(v);
        final EmptyPathHaplotypeFinderNode finder = new EmptyPathHaplotypeFinderNode(g, v);
        Assert.assertTrue(Double.isNaN(finder.score("foo".getBytes(), 1, 2)));
    }

    @Test
    public void testScore1(){
        final SeqGraph g = new SeqGraph(5);
        final SeqVertex v = new SeqVertex("acgt");
        g.addVertex(v);
        final EmptyPathHaplotypeFinderNode finder = new EmptyPathHaplotypeFinderNode(g, v);
        Assert.assertTrue(Double.isNaN(finder.score("acgtgattaca".getBytes(), 2, 4)));
    }

    @Test
    public void testScore2(){
        final SeqGraph g = new SeqGraph(5);
        final SeqVertex v = new SeqVertex("tga");
        g.addVertex(v);
        final EmptyPathHaplotypeFinderNode finder = new EmptyPathHaplotypeFinderNode(g, v);
        Assert.assertEquals(finder.score("acgtgattaca".getBytes(), 3, 3), 0.0);
    }

}
