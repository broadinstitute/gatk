package org.broadinstitute.hellbender.tools.spark.linkedreads;

import htsjdk.samtools.SAMFileHeader;
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

        final Map<String, IntervalTree<List<GATKRead>>> intervals = new HashMap<>();

        final GATKRead read1 = ArtificialReadUtils.createSamBackedRead("A", "1", 1750, 100);
        final Map<String, IntervalTree<List<GATKRead>>> stringIntervalTreeMap = CollectLinkedReadCoverageSpark.addReadToIntervals(intervals, read1, 1000);

        Assert.assertNotNull(stringIntervalTreeMap.get("1"));
        Assert.assertEquals(stringIntervalTreeMap.get("1").size(), 1);
        final IntervalTree.Node<List<GATKRead>> node = stringIntervalTreeMap.get("1").find(1750, 1849);
        Assert.assertNotNull(node);
        Assert.assertEquals(node.getValue().size(), 1);

        final GATKRead read2 = ArtificialReadUtils.createSamBackedRead("B", "1", 1900, 100);
        final Map<String, IntervalTree<List<GATKRead>>> stringIntervalTreeMap2 = CollectLinkedReadCoverageSpark.addReadToIntervals(intervals, read2, 1000);
        Assert.assertNotNull(stringIntervalTreeMap2.get("1"));
        Assert.assertEquals(stringIntervalTreeMap2.get("1").size(), 1);
        final IntervalTree.Node<List<GATKRead>> node2 = stringIntervalTreeMap2.get("1").find(1750, 1999);
        Assert.assertNotNull(node2);
        Assert.assertEquals(node2.getValue().size(), 2);

        final GATKRead read3 = ArtificialReadUtils.createSamBackedRead("C", "2", 500, 100);
        final Map<String, IntervalTree<List<GATKRead>>> stringIntervalTreeMap3 = CollectLinkedReadCoverageSpark.addReadToIntervals(intervals, read3, 1000);
        Assert.assertNotNull(stringIntervalTreeMap3.get("1"));
        Assert.assertEquals(stringIntervalTreeMap3.get("1").size(), 1);
        final IntervalTree.Node<List<GATKRead>> node3 = stringIntervalTreeMap3.get("1").find(1750, 1999);
        Assert.assertNotNull(node3);
        Assert.assertEquals(node3.getValue().size(), 2);

        Assert.assertNotNull(stringIntervalTreeMap3.get("2"));
        Assert.assertEquals(stringIntervalTreeMap3.get("2").size(), 1);
        final IntervalTree.Node<List<GATKRead>> node4 = stringIntervalTreeMap3.get("2").find(500, 599);
        Assert.assertNotNull(node4);
        Assert.assertEquals(node4.getValue().size(), 1);

    }

    @Test
    public void testMergeIntervalTrees() {
        final IntervalTree<List<GATKRead>> tree1 = new IntervalTree<>();
        final IntervalTree<List<GATKRead>> tree2 = new IntervalTree<>();

        final List<GATKRead> resultList1 = new ArrayList<>();
        final SAMFileHeader artificialSamHeader = ArtificialReadUtils.createArtificialSamHeader(1, 1, 10000);
        resultList1.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1000, 100));
        resultList1.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1004, 100));
        resultList1.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1010, 100));
        resultList1.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1100, 100));
        resultList1.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1150, 100));
        tree1.put(1000, 2000, resultList1);
        final List<GATKRead> resultList2 = new ArrayList<>();
        resultList2.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 3000, 100));
        resultList2.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 3004, 100));
        resultList2.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 3010, 100));
        tree1.put(3000, 3100, resultList2);
        final List<GATKRead> resultList3 = new ArrayList<>();
        resultList3.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1500, 100));
        resultList3.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1504, 100));
        resultList3.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1510, 100));
        resultList3.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1520, 100));
        resultList3.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1524, 100));
        resultList3.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1530, 100));
        resultList3.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1540, 100));
        resultList3.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1544, 100));
        resultList3.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1550, 100));
        tree2.put(1500, 2001, resultList3);
        final List<GATKRead> resultList4 = new ArrayList<>();
        resultList4.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 3300, 200));
        tree2.put(3300, 3500, resultList4);

        final IntervalTree<List<GATKRead>> mergedTree = CollectLinkedReadCoverageSpark.mergeIntervalTrees(tree1, tree2, 500);
        Assert.assertEquals(mergedTree.size(), 2);
        final IntervalTree.Node<List<GATKRead>> mergedNode1 = mergedTree.find(1000, 2001);
        Assert.assertNotNull(mergedNode1);
        Assert.assertEquals(mergedNode1.getValue().size(), 14);

        final IntervalTree.Node<List<GATKRead>> mergedNode2 = mergedTree.find(3000, 3500);
        Assert.assertNotNull(mergedNode2);
        Assert.assertEquals(mergedNode2.getValue().size(), 4);

    }
}