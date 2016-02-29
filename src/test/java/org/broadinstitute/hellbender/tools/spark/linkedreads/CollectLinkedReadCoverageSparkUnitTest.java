package org.broadinstitute.hellbender.tools.spark.linkedreads;

import htsjdk.samtools.util.IntervalTree;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.HashMap;
import java.util.Map;

public class CollectLinkedReadCoverageSparkUnitTest {

    @Test
    public void testAddReadToIntervals() {

        Map<String, IntervalTree<Integer>> intervals = new HashMap<>();

        final GATKRead read1 = ArtificialReadUtils.createSamBackedRead("Foo", "1", 1750, 100);
        final Map<String, IntervalTree<Integer>> stringIntervalTreeMap = CollectLinkedReadCoverageSpark.addReadToIntervals(intervals, read1);

        Assert.assertNotNull(stringIntervalTreeMap.get("1"));
        Assert.assertEquals(stringIntervalTreeMap.get("1").size(), 1);
        final IntervalTree.Node<Integer> node = stringIntervalTreeMap.get("1").find(1750, 1849);
        Assert.assertNotNull(node);
        Assert.assertEquals((int) node.getValue(), 1);

        final GATKRead read2 = ArtificialReadUtils.createSamBackedRead("Bar", "1", 1900, 100);
        final Map<String, IntervalTree<Integer>> stringIntervalTreeMap2 = CollectLinkedReadCoverageSpark.addReadToIntervals(intervals, read1);
        Assert.assertNotNull(stringIntervalTreeMap2.get("1"));
        Assert.assertEquals(stringIntervalTreeMap2.get("1").size(), 1);
        final IntervalTree.Node<Integer> node2 = stringIntervalTreeMap2.get("1").find(1750, 1999);
        Assert.assertNotNull(node2);
        Assert.assertEquals((int) node2.getValue(), 2);

    }

    @Test
    public void testMergeIntervalTrees() {
        IntervalTree<Integer> tree1 = new IntervalTree<>();
        IntervalTree<Integer> tree2 = new IntervalTree<>();

        tree1.put(1000, 2000, 5);
        tree1.put(3000, 3100, 3);
        tree2.put(1500, 2001, 8);

        final IntervalTree<Integer> mergedTree = CollectLinkedReadCoverageSpark.mergeIntervalTrees(tree1, tree2);
        Assert.assertEquals(mergedTree.size(), 2);
        final IntervalTree.Node<Integer> mergedNode1 = mergedTree.find(1000, 2001);
        Assert.assertNotNull(mergedNode1);
        Assert.assertEquals((int) mergedNode1.getValue(), 13);

        final IntervalTree.Node<Integer> mergedNode2 = mergedTree.find(3000, 3100);
        Assert.assertNotNull(mergedNode2);
        Assert.assertEquals((int) mergedNode2.getValue(), 3);

    }
}