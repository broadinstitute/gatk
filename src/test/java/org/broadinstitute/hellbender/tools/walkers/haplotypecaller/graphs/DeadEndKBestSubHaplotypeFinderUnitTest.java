package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Collections;

public final class DeadEndKBestSubHaplotypeFinderUnitTest {

    @Test
    public void testCount(){
        final KBestSubHaplotypeFinder finder = DeadEndKBestSubHaplotypeFinder.INSTANCE;
        Assert.assertEquals(finder.getCount(), 0);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testKneg(){
        final KBestSubHaplotypeFinder finder = DeadEndKBestSubHaplotypeFinder.INSTANCE;
        finder.getKBest(-1);
    }
    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testKpos(){
        final KBestSubHaplotypeFinder finder = DeadEndKBestSubHaplotypeFinder.INSTANCE;
        finder.getKBest(1);
    }
    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testKzero(){
        final KBestSubHaplotypeFinder finder = DeadEndKBestSubHaplotypeFinder.INSTANCE;
        finder.getKBest(0);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testScoreNullBases(){
        final KBestSubHaplotypeFinder finder = DeadEndKBestSubHaplotypeFinder.INSTANCE;
        finder.score(null, 1, 1);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testScoreNegOffset(){
        final KBestSubHaplotypeFinder finder = DeadEndKBestSubHaplotypeFinder.INSTANCE;
        finder.score(new byte[]{'a','c','g','t'}, -1, 1);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testScoreNegLength(){
        final KBestSubHaplotypeFinder finder = DeadEndKBestSubHaplotypeFinder.INSTANCE;
        finder.score(new byte[]{'a','c','g','t'}, 1, -1);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testScoreEnd(){
        final KBestSubHaplotypeFinder finder = DeadEndKBestSubHaplotypeFinder.INSTANCE;
        finder.score(new byte[]{'a','c','g','t'}, 2, 3);
    }

    @Test
    public void testScoreEndOK(){
        final KBestSubHaplotypeFinder finder = DeadEndKBestSubHaplotypeFinder.INSTANCE;
        Assert.assertTrue(Double.isNaN(finder.score(new byte[]{'a','c','g','t'}, 2, 1)));
    }

    @Test
    public void testID(){
        final KBestSubHaplotypeFinder finder = DeadEndKBestSubHaplotypeFinder.INSTANCE;
        Assert.assertEquals(finder.id(), "<DEAD>");
    }

    @Test
    public void testLabel(){
        final KBestSubHaplotypeFinder finder = DeadEndKBestSubHaplotypeFinder.INSTANCE;
        Assert.assertEquals(finder.label(), "&lt;DEAD&gt;");
    }

    @Test
    public void testsubFinderLabels(){
        final KBestSubHaplotypeFinder finder = DeadEndKBestSubHaplotypeFinder.INSTANCE;
        Assert.assertEquals(finder.subFinderLabels(), Collections.emptySet());
    }

    @Test
    public void testIsReference(){
        final KBestSubHaplotypeFinder finder = DeadEndKBestSubHaplotypeFinder.INSTANCE;
        Assert.assertEquals(finder.isReference(), false);
    }

}
