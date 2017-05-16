package org.broadinstitute.hellbender.tools.spark.linkedreads;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.tools.spark.linkedreads.CollectLinkedReadCoverageSpark.ReadInfo;
import org.broadinstitute.hellbender.tools.spark.sv.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.SVIntervalTree;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class CollectLinkedReadCoverageSparkUnitTest {

    @Test
    public void testAddReadToIntervals() {
        
        final SVIntervalTree<List<ReadInfo>> intervals = new SVIntervalTree<>();

        final GATKRead read1 = ArtificialReadUtils.createSamBackedRead("A", "1", 1750, 100);
        final ReadMetadata readMetadata = new ReadMetadata(Collections.emptySet(), ArtificialReadUtils.createArtificialSamHeader(), null, 1, 1, 1, 30);
        final SVIntervalTree<List<ReadInfo>> svIntervalTree = CollectLinkedReadCoverageSpark.addReadToIntervals(intervals, new ReadInfo(readMetadata, read1), 1000);

        final SVIntervalTree.Entry<List<ReadInfo>> node = svIntervalTree.find(new SVInterval(readMetadata.getContigID("1"), 1750, 1849));
        Assert.assertNotNull(node);
        Assert.assertEquals(node.getInterval().getStart(), 1750);
        Assert.assertEquals(node.getInterval().getEnd(), 1849);
        Assert.assertEquals(node.getValue().size(), 1);

        final GATKRead read2 = ArtificialReadUtils.createSamBackedRead("B", "1", 1900, 100);
        final SVIntervalTree<List<ReadInfo>> svIntervalTree2 = CollectLinkedReadCoverageSpark.addReadToIntervals(intervals, new ReadInfo(readMetadata, read2), 1000);
        final SVIntervalTree.Entry<List<ReadInfo>> node2 = svIntervalTree2.find(new SVInterval(readMetadata.getContigID("1"), 1750, 1999));
        Assert.assertNotNull(node2);
        Assert.assertEquals(node2.getInterval().getStart(), 1750);
        Assert.assertEquals(node2.getInterval().getEnd(), 1999);
        Assert.assertEquals(node2.getValue().size(), 2);

        final GATKRead read3 = ArtificialReadUtils.createSamBackedRead("D", "1", 3500, 100);
        final SVIntervalTree<List<ReadInfo>> stringSVIntervalTree3 = CollectLinkedReadCoverageSpark.addReadToIntervals(intervals, new ReadInfo(readMetadata, read3), 1000);
        final SVIntervalTree.Entry<List<ReadInfo>> node3 = stringSVIntervalTree3.find(new SVInterval(readMetadata.getContigID("1"), 1750, 1999));
        Assert.assertNotNull(node3);
        Assert.assertEquals(node3.getInterval().getStart(), 1750);
        Assert.assertEquals(node3.getInterval().getEnd(), 1999);
        Assert.assertEquals(node3.getValue().size(), 2);

        final SVIntervalTree.Entry<List<ReadInfo>> node4 = stringSVIntervalTree3.find(new SVInterval(readMetadata.getContigID("1"), 3500, 3599));
        Assert.assertNotNull(node4);
        Assert.assertEquals(node4.getInterval().getStart(), 3500);
        Assert.assertEquals(node4.getInterval().getEnd(), 3599);
        Assert.assertEquals(node4.getValue().size(), 1);

        final GATKRead read4 = ArtificialReadUtils.createSamBackedRead("C", "2", 500, 100);
        final SVIntervalTree<List<ReadInfo>> svIntervalTree4 = CollectLinkedReadCoverageSpark.addReadToIntervals(intervals, new ReadInfo(readMetadata, read4), 2500);
        final SVIntervalTree.Entry<List<ReadInfo>> node5 = svIntervalTree4.find(new SVInterval(readMetadata.getContigID("1"), 1750, 1999));
        Assert.assertNotNull(node5);
        Assert.assertEquals(node5.getValue().size(), 2);

        final SVIntervalTree.Entry<List<ReadInfo>> node6 = svIntervalTree4.find(new SVInterval(readMetadata.getContigID("2"), 500, 599));
        Assert.assertNotNull(node6);
        Assert.assertEquals(node6.getValue().size(), 1);

    }

    @Test
    public void testMergeSVIntervalTrees() {
        final SAMFileHeader artificialSamHeader = ArtificialReadUtils.createArtificialSamHeader(1, 1, 10000);
        final ReadMetadata readMetadata = new ReadMetadata(Collections.emptySet(), ArtificialReadUtils.createArtificialSamHeader(), null, 1, 1, 1, 30);
        final SVIntervalTree<List<ReadInfo>> tree1 = createTestSVIntervalTree1(artificialSamHeader, readMetadata);
        final SVIntervalTree<List<ReadInfo>> tree2 = new SVIntervalTree<>();

        final List<ReadInfo> resultList2 = new ArrayList<>();
        resultList2.add(new ReadInfo(readMetadata, ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 3000, 100)));
        resultList2.add(new ReadInfo(readMetadata, ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 3004, 100)));
        resultList2.add(new ReadInfo(readMetadata, ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 3010, 100)));
        tree1.put(new SVInterval(readMetadata.getContigID("1"), 3000, 3100), resultList2);

        final List<ReadInfo> resultList3 = new ArrayList<>();
        resultList3.add(new ReadInfo(readMetadata, ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1500, 100)));
        resultList3.add(new ReadInfo(readMetadata, ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1504, 100)));
        resultList3.add(new ReadInfo(readMetadata, ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1510, 100)));
        resultList3.add(new ReadInfo(readMetadata, ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1520, 100)));
        resultList3.add(new ReadInfo(readMetadata, ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1524, 100)));
        resultList3.add(new ReadInfo(readMetadata, ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1530, 100)));
        resultList3.add(new ReadInfo(readMetadata, ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1540, 100)));
        resultList3.add(new ReadInfo(readMetadata, ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1544, 100)));
        resultList3.add(new ReadInfo(readMetadata, ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1550, 100)));
        tree2.put(new SVInterval(readMetadata.getContigID("1"), 1500, 1649), resultList3);
        final List<ReadInfo> resultList4 = new ArrayList<>();
        resultList4.add(new ReadInfo(readMetadata, ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 3300, 200)));
        tree2.put(new SVInterval(readMetadata.getContigID("1"), 3300, 3500), resultList4);

        final SVIntervalTree<List<ReadInfo>> mergedTree = CollectLinkedReadCoverageSpark.mergeIntervalTrees(tree1, tree2, 250);
        Assert.assertEquals(mergedTree.size(), 2);
        final SVIntervalTree.Entry<List<ReadInfo>> mergedNode1 = mergedTree.find(new SVInterval(readMetadata.getContigID("1"), 1000, 1649));
        Assert.assertNotNull(mergedNode1);
        Assert.assertEquals(mergedNode1.getInterval().getStart(), 1000);
        Assert.assertEquals(mergedNode1.getInterval().getEnd(), 1649);
        Assert.assertEquals(mergedNode1.getValue().size(), 14);

        final SVIntervalTree.Entry<List<ReadInfo>> mergedNode2 = mergedTree.find(new SVInterval(readMetadata.getContigID("1"), 3000, 3500));
        Assert.assertNotNull(mergedNode2);
        Assert.assertEquals(mergedNode2.getInterval().getStart(), 3000);
        Assert.assertEquals(mergedNode2.getInterval().getEnd(), 3500);

        Assert.assertEquals(mergedNode2.getValue().size(), 4);

    }

    @Test
    public void testSVIntervalTreeToBed() {
        final SAMFileHeader artificialSamHeader = ArtificialReadUtils.createArtificialSamHeader(1, 1, 10000);
        final ReadMetadata readMetadata = new ReadMetadata(Collections.emptySet(), ArtificialReadUtils.createArtificialSamHeader(), null, 1, 1, 1, 30);
        final SVIntervalTree<List<ReadInfo>> tree = createTestSVIntervalTree1(artificialSamHeader, readMetadata);

        final String barcode = "ACTGACTG";

        final String bedRecord = CollectLinkedReadCoverageSpark.intervalTreeToBedRecord(barcode, tree.iterator().next(), readMetadata);
        Assert.assertEquals(bedRecord, "1\t1000\t1160\tACTGACTG\t5\t+\t1000\t1160\t0,0,255\t5\t10,11,10,9,10\t0,20,25,45,150");

    }

    private SVIntervalTree<List<ReadInfo>> createTestSVIntervalTree1(final SAMFileHeader artificialSamHeader, final ReadMetadata readMetadata) {
        final SVIntervalTree<List<ReadInfo>> tree = new SVIntervalTree<>();
        final List<ReadInfo> resultList1 = new ArrayList<>();
        resultList1.add(new ReadInfo(readMetadata, ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1000, "ACACACACAC".getBytes(), "2222222222".getBytes(), "10M")));
        resultList1.add(new ReadInfo(readMetadata, ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1020, "GTGTGTGTGT".getBytes(), "3333333333".getBytes(), "6M1D4M")));
        resultList1.add(new ReadInfo(readMetadata, ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1025, "AGAGAGAGAG".getBytes(), "4444444444".getBytes(), "10M")));
        resultList1.add(new ReadInfo(readMetadata, ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1045, "GCGCGCGCGC".getBytes(), "5555555555".getBytes(), "5M1I4M")));
        resultList1.add(new ReadInfo(readMetadata, ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1150, "TTTTTTTTTT".getBytes(), "6666666666".getBytes(), "10M")));
        tree.put(new SVInterval(readMetadata.getContigID("1"), 1000, 1160), resultList1);
        return tree;
    }

}