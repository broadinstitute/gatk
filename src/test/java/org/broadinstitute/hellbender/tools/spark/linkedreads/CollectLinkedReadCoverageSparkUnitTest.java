package org.broadinstitute.hellbender.tools.spark.linkedreads;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.IntervalTree;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;
import scala.Tuple2;

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
        final SAMFileHeader artificialSamHeader = ArtificialReadUtils.createArtificialSamHeader(1, 1, 10000);
        final IntervalTree<List<GATKRead>> tree1 = createTestIntervalTree1(artificialSamHeader);
        final IntervalTree<List<GATKRead>> tree2 = new IntervalTree<>();

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

    @Test
    public void testIntervalTreeToBed() {
        final IntervalTree<List<GATKRead>> tree = createTestIntervalTree1(ArtificialReadUtils.createArtificialSamHeader(1, 1, 10000));

        final String barcode = "ACTGACTG";
        final Map<String, IntervalTree<List<GATKRead>>> treeForBarcode = new HashMap<>();
        treeForBarcode.put("1", tree);

        final List<String> bedRecords = CollectLinkedReadCoverageSpark.barcodeLocationsToBed(new Tuple2<>(barcode, treeForBarcode));
        Assert.assertEquals(bedRecords.size(), 1);
        Assert.assertEquals(bedRecords.get(0), "1\t1000\t1160\tACTGACTG\t5\t+\t1000\t1160\t0,0,255\t5\t10,11,10,9,10\t0,20,25,45,150");

    }

    @Test
    public void testIntervalTreeToSAM() {

        final SAMFileHeader samHeader = ArtificialReadUtils.createArtificialSamHeader(1, 1, 10000);


        final IntervalTree<List<GATKRead>> tree2 = createTestIntervalTree2(samHeader);
        final String barcode2 = "GGGGGGG";
        final Map<String, IntervalTree<List<GATKRead>>> treeForBarcode2 = new HashMap<>();
        treeForBarcode2.put("1", tree2);

        final List<String> samRecords2 = CollectLinkedReadCoverageSpark.barcodeLocationsToSam(new Tuple2<>(barcode2, treeForBarcode2), samHeader);
        Assert.assertEquals(samRecords2.size(), 1);
        Assert.assertEquals(samRecords2.get(0), "GGGGGGG\t1\t1\t997\t60\t3M7M18D8M\t*\t0\t0\tACACACACACGTGTGTGT\tSSSSSSSSSSTTTTTTTT\n");

        final IntervalTree<List<GATKRead>> tree = createTestIntervalTree1(samHeader);
        final String barcode = "ACTGACTG";
        final Map<String, IntervalTree<List<GATKRead>>> treeForBarcode = new HashMap<>();
        treeForBarcode.put("1", tree);

        final List<String> samRecords = CollectLinkedReadCoverageSpark.barcodeLocationsToSam(new Tuple2<>(barcode, treeForBarcode), samHeader);
        Assert.assertEquals(samRecords.size(), 1);
        Assert.assertEquals(samRecords.get(0), "ACTGACTG\t1\t1\t1000\t60\t10M10N6M1D4M4M10N5M1I4M96N10M\t*\t0\t0\tACACACACACGTGTGTGTGTAGAGGCGCGCGCGCTTTTTTTTTT\tSSSSSSSSSSTTTTTTTTTTUUUUVVVVVVVVVVWWWWWWWWWW\n");


    }

    private IntervalTree<List<GATKRead>> createTestIntervalTree1(final SAMFileHeader artificialSamHeader) {
        final IntervalTree<List<GATKRead>> tree = new IntervalTree<>();

        final List<GATKRead> resultList1 = new ArrayList<>();
        resultList1.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1000, "ACACACACAC".getBytes(), "2222222222".getBytes(), "10M"));
        resultList1.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1020, "GTGTGTGTGT".getBytes(), "3333333333".getBytes(), "6M1D4M"));
        resultList1.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1025, "AGAGAGAGAG".getBytes(), "4444444444".getBytes(), "10M"));
        resultList1.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1045, "GCGCGCGCGC".getBytes(), "5555555555".getBytes(), "5M1I4M"));
        resultList1.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1150, "TTTTTTTTTT".getBytes(), "6666666666".getBytes(), "10M"));
        tree.put(1000, 1160, resultList1);
        return tree;
    }

    private IntervalTree<List<GATKRead>> createTestIntervalTree2(final SAMFileHeader artificialSamHeader) {
        final IntervalTree<List<GATKRead>> tree = new IntervalTree<>();

        final List<GATKRead> resultList = new ArrayList<>();
        // reads are like:
        //     acaCACACAC
        //           GT--------------------GTGTGTGT
        resultList.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1000, "ACACACACAC".getBytes(), "2222222222".getBytes(), "3S7M"));
        resultList.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1003, "GTGTGTGTGT".getBytes(), "3333333333".getBytes(), "2M20D8M"));
        tree.put(1000, 1033, resultList);
        return tree;
    }

}