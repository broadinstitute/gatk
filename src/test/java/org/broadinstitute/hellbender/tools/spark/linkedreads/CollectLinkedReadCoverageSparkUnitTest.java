package org.broadinstitute.hellbender.tools.spark.linkedreads;

import htsjdk.samtools.util.IntervalTree;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class CollectLinkedReadCoverageSparkUnitTest {

    @Test
    public void testAddReadToIntervals() {

        final Map<String, IntervalTree<List<CollectLinkedReadCoverageSpark.ReadInfo>>> intervals = new HashMap<>();

        final GATKRead read1 = ArtificialReadUtils.createSamBackedRead("A", "1", 1750, 100);
        final Map<String, IntervalTree<List<CollectLinkedReadCoverageSpark.ReadInfo>>> stringIntervalTreeMap = CollectLinkedReadCoverageSpark.addReadToIntervals(intervals, read1, 1000);

        Assert.assertNotNull(stringIntervalTreeMap.get("1"));
        Assert.assertEquals(stringIntervalTreeMap.get("1").size(), 1);
        final IntervalTree.Node<List<CollectLinkedReadCoverageSpark.ReadInfo>> node = stringIntervalTreeMap.get("1").find(1750, 1849);
        Assert.assertNotNull(node);
        Assert.assertEquals(node.getValue().size(), 1);

        final GATKRead read2 = ArtificialReadUtils.createSamBackedRead("B", "1", 1900, 100);
        final Map<String, IntervalTree<List<CollectLinkedReadCoverageSpark.ReadInfo>>> stringIntervalTreeMap2 = CollectLinkedReadCoverageSpark.addReadToIntervals(intervals, read2, 1000);
        Assert.assertNotNull(stringIntervalTreeMap2.get("1"));
        Assert.assertEquals(stringIntervalTreeMap2.get("1").size(), 1);
        final IntervalTree.Node<List<CollectLinkedReadCoverageSpark.ReadInfo>> node2 = stringIntervalTreeMap2.get("1").find(1750, 1999);
        Assert.assertNotNull(node2);
        Assert.assertEquals(node2.getValue().size(), 2);

        final GATKRead read3 = ArtificialReadUtils.createSamBackedRead("C", "2", 500, 100);
        final Map<String, IntervalTree<List<CollectLinkedReadCoverageSpark.ReadInfo>>> stringIntervalTreeMap3 = CollectLinkedReadCoverageSpark.addReadToIntervals(intervals, read3, 1000);
        Assert.assertNotNull(stringIntervalTreeMap3.get("1"));
        Assert.assertEquals(stringIntervalTreeMap3.get("1").size(), 1);
        final IntervalTree.Node<List<CollectLinkedReadCoverageSpark.ReadInfo>> node3 = stringIntervalTreeMap3.get("1").find(1750, 1999);
        Assert.assertNotNull(node3);
        Assert.assertEquals(node3.getValue().size(), 2);

        Assert.assertNotNull(stringIntervalTreeMap3.get("2"));
        Assert.assertEquals(stringIntervalTreeMap3.get("2").size(), 1);
        final IntervalTree.Node<List<CollectLinkedReadCoverageSpark.ReadInfo>> node4 = stringIntervalTreeMap3.get("2").find(500, 599);
        Assert.assertNotNull(node4);
        Assert.assertEquals(node4.getValue().size(), 1);

    }

    @Test
    public void testMergeIntervalTrees() {
        IntervalTree<List<CollectLinkedReadCoverageSpark.ReadInfo>> tree1 = new IntervalTree<>();
        IntervalTree<List<CollectLinkedReadCoverageSpark.ReadInfo>> tree2 = new IntervalTree<>();

        List<CollectLinkedReadCoverageSpark.ReadInfo> resultList1 = new ArrayList<>();
        resultList1.add(new CollectLinkedReadCoverageSpark.ReadInfo(1000, 1100));
        resultList1.add(new CollectLinkedReadCoverageSpark.ReadInfo(1004, 1104));
        resultList1.add(new CollectLinkedReadCoverageSpark.ReadInfo(1010, 1110));
        resultList1.add(new CollectLinkedReadCoverageSpark.ReadInfo(1100, 1200));
        resultList1.add(new CollectLinkedReadCoverageSpark.ReadInfo(1150, 1250));
        tree1.put(1000, 2000, resultList1);
        List<CollectLinkedReadCoverageSpark.ReadInfo> resultList2 = new ArrayList<>();
        resultList2.add(new CollectLinkedReadCoverageSpark.ReadInfo(3000, 3100));
        resultList2.add(new CollectLinkedReadCoverageSpark.ReadInfo(3004, 3104));
        resultList2.add(new CollectLinkedReadCoverageSpark.ReadInfo(3010, 3110));
        tree1.put(3000, 3100, resultList2);
        List<CollectLinkedReadCoverageSpark.ReadInfo> resultList3 = new ArrayList<>();
        resultList3.add(new CollectLinkedReadCoverageSpark.ReadInfo(1500, 1600));
        resultList3.add(new CollectLinkedReadCoverageSpark.ReadInfo(1504, 1604));
        resultList3.add(new CollectLinkedReadCoverageSpark.ReadInfo(1510, 1610));
        resultList3.add(new CollectLinkedReadCoverageSpark.ReadInfo(1520, 1620));
        resultList3.add(new CollectLinkedReadCoverageSpark.ReadInfo(1524, 1624));
        resultList3.add(new CollectLinkedReadCoverageSpark.ReadInfo(1530, 1630));
        resultList3.add(new CollectLinkedReadCoverageSpark.ReadInfo(1540, 1640));
        resultList3.add(new CollectLinkedReadCoverageSpark.ReadInfo(1544, 1644));
        resultList3.add(new CollectLinkedReadCoverageSpark.ReadInfo(1550, 1650));
        tree2.put(1500, 2001, resultList3);
        List<CollectLinkedReadCoverageSpark.ReadInfo> resultList4 = new ArrayList<>();
        resultList4.add(new CollectLinkedReadCoverageSpark.ReadInfo(3300, 3500));
        tree2.put(3300, 3500, resultList4);

        final IntervalTree<List<CollectLinkedReadCoverageSpark.ReadInfo>> mergedTree = CollectLinkedReadCoverageSpark.mergeIntervalTrees(tree1, tree2, 500);
        Assert.assertEquals(mergedTree.size(), 2);
        final IntervalTree.Node<List<CollectLinkedReadCoverageSpark.ReadInfo>> mergedNode1 = mergedTree.find(1000, 2001);
        Assert.assertNotNull(mergedNode1);
        Assert.assertEquals(mergedNode1.getValue().size(), 13);

        final IntervalTree.Node<List<CollectLinkedReadCoverageSpark.ReadInfo>> mergedNode2 = mergedTree.find(3000, 3500);
        Assert.assertNotNull(mergedNode2);
        Assert.assertEquals(mergedNode2.getValue().size(), 4);

    }
}