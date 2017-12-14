package org.broadinstitute.hellbender.tools.spark.sv.utils;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.util.Iterator;

public class PairedStrandedIntervalTreeTest extends GATKBaseTest {

    @Test(groups = "sv")
    public void testOverlappers() {
        PairedStrandedIntervalTree<Integer> psiTree = new PairedStrandedIntervalTree<>();

        PairedStrandedIntervals v1 =
                new PairedStrandedIntervals(
                        new StrandedInterval(new SVInterval(1, 100, 200),
                        true),
                        new StrandedInterval(new SVInterval(1, 500, 600),
                        false));

        PairedStrandedIntervals v2 =
                new PairedStrandedIntervals(
                        new StrandedInterval(new SVInterval(1, 10, 20),
                        false),
                        new StrandedInterval(new SVInterval(1, 800, 900),
                        true));

        PairedStrandedIntervals v3 =
                new PairedStrandedIntervals(
                        new StrandedInterval(new SVInterval(1, 205, 305),
                        true),
                        new StrandedInterval(new SVInterval(1, 625, 725),
                        false));

        PairedStrandedIntervals query =
                new PairedStrandedIntervals(
                        new StrandedInterval(new SVInterval(1, 110, 210),
                        true),
                        new StrandedInterval(new SVInterval(1, 550, 650),
                        true));

        psiTree.put(v1, 1);
        psiTree.put(v2, 2);
        psiTree.put(v3, 3);

        Iterator<Tuple2<PairedStrandedIntervals, Integer>> overlappers = psiTree.overlappers(query);
        Assert.assertFalse(overlappers.hasNext());

        PairedStrandedIntervals query2 =
                new PairedStrandedIntervals(
                        new StrandedInterval(new SVInterval(1, 110, 210),
                                true),
                        new StrandedInterval(new SVInterval(1, 550, 650),
                                false));

        overlappers = psiTree.overlappers(query2);
        Assert.assertTrue(overlappers.hasNext());

        Tuple2<PairedStrandedIntervals, Integer> next = overlappers.next();
        Assert.assertEquals(next._1(), v1);
        Assert.assertEquals(next._2().intValue(), 1);
        overlappers.remove();

        Assert.assertTrue(overlappers.hasNext());
        next = overlappers.next();
        Assert.assertEquals(next._1(), v3);
        Assert.assertEquals(next._2().intValue(), 3);

        Assert.assertFalse(overlappers.hasNext());

        overlappers = psiTree.overlappers(query2);

        Assert.assertTrue(overlappers.hasNext());
        next = overlappers.next();
        Assert.assertEquals(next._1(), v3);
        Assert.assertEquals(next._2().intValue(), 3);

        Assert.assertFalse(overlappers.hasNext());


    }

    @Test(groups = "sv")
    public void testIterator() {
        PairedStrandedIntervalTree<Integer> psiTree = new PairedStrandedIntervalTree<>();

        PairedStrandedIntervals v1 =
                new PairedStrandedIntervals(
                        new StrandedInterval(new SVInterval(1, 100, 200),
                        true),
                        new StrandedInterval(new SVInterval(1, 500, 600),
                        false));

        PairedStrandedIntervals v2 =
                new PairedStrandedIntervals(
                        new StrandedInterval(new SVInterval(1, 10, 20),
                        false),
                        new StrandedInterval(new SVInterval(1, 800, 900),
                        true));

        psiTree.put(v1, 1);
        psiTree.put(v2, 2);

        Iterator<Tuple2<PairedStrandedIntervals, Integer>> iterator = psiTree.iterator();
        Assert.assertTrue(iterator.hasNext());
        Tuple2<PairedStrandedIntervals, Integer> out1 = iterator.next();
        Assert.assertEquals(v2, out1._1());
        Assert.assertEquals(2, out1._2().intValue());
        Assert.assertTrue(iterator.hasNext());
        Tuple2<PairedStrandedIntervals, Integer> out2 = iterator.next();
        Assert.assertEquals(v1, out2._1());
        Assert.assertEquals(1, out2._2().intValue());

    }

    @Test(groups = "sv")
    public void testIteratorRemove() {
        final PairedStrandedIntervalTree<Integer> psiTree = new PairedStrandedIntervalTree<>();

        final PairedStrandedIntervals v1 =
                new PairedStrandedIntervals(
                        new StrandedInterval(new SVInterval(1, 100, 200),
                        true),
                        new StrandedInterval(new SVInterval(1, 500, 600),
                        false));

        final PairedStrandedIntervals v2 =
                new PairedStrandedIntervals(
                        new StrandedInterval(new SVInterval(1, 10, 20),
                        false),
                        new StrandedInterval(new SVInterval(1, 800, 900),
                        true));

        psiTree.put(v1, 1);
        psiTree.put(v2, 2);

        Iterator<Tuple2<PairedStrandedIntervals, Integer>> iterator = psiTree.iterator();
        Assert.assertTrue(iterator.hasNext());
        Tuple2<PairedStrandedIntervals, Integer> out1 = iterator.next();
        Assert.assertEquals(v2, out1._1());
        Assert.assertEquals(2, out1._2().intValue());

        iterator.remove();

        Assert.assertTrue(iterator.hasNext());
        Tuple2<PairedStrandedIntervals, Integer> out2 = iterator.next();
        Assert.assertEquals(v1, out2._1());
        Assert.assertEquals(1, out2._2().intValue());

        iterator = psiTree.iterator();
        Assert.assertTrue(iterator.hasNext());
        Tuple2<PairedStrandedIntervals, Integer> out3 = iterator.next();
        Assert.assertEquals(v1, out3._1());
        Assert.assertEquals(1, out3._2().intValue());

        iterator.remove();
        Assert.assertFalse(iterator.hasNext());

        iterator = psiTree.iterator();
        Assert.assertFalse(iterator.hasNext());


    }

    @Test(groups = "sv")
    public void testContains() {
        PairedStrandedIntervals query =
                new PairedStrandedIntervals(
                        new StrandedInterval(new SVInterval(1, 100, 200),
                        true),
                        new StrandedInterval(new SVInterval(1, 500, 600),
                        false));
        PairedStrandedIntervalTree<Integer> psiTree = new PairedStrandedIntervalTree<>();
        psiTree.put(query, 1);

        Assert.assertTrue(psiTree.contains(query));

        PairedStrandedIntervals query2 =
                new PairedStrandedIntervals(
                        new StrandedInterval(new SVInterval(1, 100, 200),
                        true),
                        new StrandedInterval(new SVInterval(1, 500, 600),
                        true));

        Assert.assertFalse(psiTree.contains(query2));
    }

}