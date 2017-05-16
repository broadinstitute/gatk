package org.broadinstitute.hellbender.tools.spark.linkedreads;

import htsjdk.samtools.SAMFileHeader;
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
        
        final SVIntervalTree<List<GATKRead>> intervals = new SVIntervalTree<>();

        final GATKRead read1 = ArtificialReadUtils.createSamBackedRead("A", "1", 1750, 100);
        final ReadMetadata readMetadata = new ReadMetadata(Collections.emptySet(), ArtificialReadUtils.createArtificialSamHeader(), null, 1, 1, 1, 30);
        final SVIntervalTree<List<GATKRead>> svIntervalTree = CollectLinkedReadCoverageSpark.addReadToIntervals(intervals, read1, 1000, readMetadata);

        final SVIntervalTree.Entry<List<GATKRead>> node = svIntervalTree.find(new SVInterval(readMetadata.getContigID("1"), 1750, 1849));
        Assert.assertNotNull(node);
        Assert.assertEquals(node.getInterval().getStart(), 1750);
        Assert.assertEquals(node.getInterval().getEnd(), 1849);
        Assert.assertEquals(node.getValue().size(), 1);

        final GATKRead read2 = ArtificialReadUtils.createSamBackedRead("B", "1", 1900, 100);
        final SVIntervalTree<List<GATKRead>> svIntervalTree2 = CollectLinkedReadCoverageSpark.addReadToIntervals(intervals, read2, 1000, readMetadata);
        final SVIntervalTree.Entry<List<GATKRead>> node2 = svIntervalTree2.find(new SVInterval(readMetadata.getContigID("1"), 1750, 1999));
        Assert.assertNotNull(node2);
        Assert.assertEquals(node2.getInterval().getStart(), 1750);
        Assert.assertEquals(node2.getInterval().getEnd(), 1999);
        Assert.assertEquals(node2.getValue().size(), 2);

        final GATKRead read3 = ArtificialReadUtils.createSamBackedRead("D", "1", 3500, 100);
        final SVIntervalTree<List<GATKRead>> stringSVIntervalTree3 = CollectLinkedReadCoverageSpark.addReadToIntervals(intervals, read3, 1000, readMetadata);
        final SVIntervalTree.Entry<List<GATKRead>> node3 = stringSVIntervalTree3.find(new SVInterval(readMetadata.getContigID("1"), 1750, 1999));
        Assert.assertNotNull(node3);
        Assert.assertEquals(node3.getInterval().getStart(), 1750);
        Assert.assertEquals(node3.getInterval().getEnd(), 1999);
        Assert.assertEquals(node3.getValue().size(), 2);

        final SVIntervalTree.Entry<List<GATKRead>> node4 = stringSVIntervalTree3.find(new SVInterval(readMetadata.getContigID("1"), 3500, 3599));
        Assert.assertNotNull(node4);
        Assert.assertEquals(node4.getInterval().getStart(), 3500);
        Assert.assertEquals(node4.getInterval().getEnd(), 3599);
        Assert.assertEquals(node4.getValue().size(), 1);

        final GATKRead read4 = ArtificialReadUtils.createSamBackedRead("C", "2", 500, 100);
        final SVIntervalTree<List<GATKRead>> svIntervalTree4 = CollectLinkedReadCoverageSpark.addReadToIntervals(intervals, read4, 2500, readMetadata);
        final SVIntervalTree.Entry<List<GATKRead>> node5 = svIntervalTree4.find(new SVInterval(readMetadata.getContigID("1"), 1750, 1999));
        Assert.assertNotNull(node5);
        Assert.assertEquals(node5.getValue().size(), 2);

        final SVIntervalTree.Entry<List<GATKRead>> node6 = svIntervalTree4.find(new SVInterval(readMetadata.getContigID("2"), 500, 599));
        Assert.assertNotNull(node6);
        Assert.assertEquals(node6.getValue().size(), 1);

    }

    @Test
    public void testMergeSVIntervalTrees() {
        final SAMFileHeader artificialSamHeader = ArtificialReadUtils.createArtificialSamHeader(1, 1, 10000);
        final ReadMetadata readMetadata = new ReadMetadata(Collections.emptySet(), ArtificialReadUtils.createArtificialSamHeader(), null, 1, 1, 1, 30);
        final SVIntervalTree<List<GATKRead>> tree1 = createTestSVIntervalTree1(artificialSamHeader, readMetadata);
        final SVIntervalTree<List<GATKRead>> tree2 = new SVIntervalTree<>();

        final List<GATKRead> resultList2 = new ArrayList<>();
        resultList2.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 3000, 100));
        resultList2.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 3004, 100));
        resultList2.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 3010, 100));
        tree1.put(new SVInterval(readMetadata.getContigID("1"), 3000, 3100), resultList2);

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
        tree2.put(new SVInterval(readMetadata.getContigID("1"), 1500, 1649), resultList3);
        final List<GATKRead> resultList4 = new ArrayList<>();
        resultList4.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 3300, 200));
        tree2.put(new SVInterval(readMetadata.getContigID("1"), 3300, 3500), resultList4);

        final SVIntervalTree<List<GATKRead>> mergedTree = CollectLinkedReadCoverageSpark.mergeIntervalTrees(tree1, tree2, 250);
        Assert.assertEquals(mergedTree.size(), 2);
        final SVIntervalTree.Entry<List<GATKRead>> mergedNode1 = mergedTree.find(new SVInterval(readMetadata.getContigID("1"), 1000, 1649));
        Assert.assertNotNull(mergedNode1);
        Assert.assertEquals(mergedNode1.getInterval().getStart(), 1000);
        Assert.assertEquals(mergedNode1.getInterval().getEnd(), 1649);
        Assert.assertEquals(mergedNode1.getValue().size(), 14);

        final SVIntervalTree.Entry<List<GATKRead>> mergedNode2 = mergedTree.find(new SVInterval(readMetadata.getContigID("1"), 3000, 3500));
        Assert.assertNotNull(mergedNode2);
        Assert.assertEquals(mergedNode2.getInterval().getStart(), 3000);
        Assert.assertEquals(mergedNode2.getInterval().getEnd(), 3500);

        Assert.assertEquals(mergedNode2.getValue().size(), 4);

    }

    @Test
    public void testSVIntervalTreeToBed() {
        final SAMFileHeader artificialSamHeader = ArtificialReadUtils.createArtificialSamHeader(1, 1, 10000);
        final ReadMetadata readMetadata = new ReadMetadata(Collections.emptySet(), ArtificialReadUtils.createArtificialSamHeader(), null, 1, 1, 1, 30);
        final SVIntervalTree<List<GATKRead>> tree = createTestSVIntervalTree1(artificialSamHeader, readMetadata);

        final String barcode = "ACTGACTG";

        final String bedRecord = CollectLinkedReadCoverageSpark.intervalTreeToBedRecord(barcode, tree.iterator().next(), readMetadata);
        Assert.assertEquals(bedRecord, "1\t1000\t1160\tACTGACTG\t5\t+\t1000\t1160\t0,0,255\t5\t10,11,10,9,10\t0,20,25,45,150");

    }


    @Test
    public void testSVIntervalTreeToSAM() {

        final SAMFileHeader samHeader = ArtificialReadUtils.createArtificialSamHeader(1, 1, 10000);
        final ReadMetadata readMetadata = new ReadMetadata(Collections.emptySet(), ArtificialReadUtils.createArtificialSamHeader(), null, 1, 1, 1, 30);

        final SVIntervalTree<List<GATKRead>> tree3 = createTestSVIntervalTree3(samHeader, readMetadata);
        final String barcode3 = "AAAAAAA";

        final GATKRead samRecord3 = CollectLinkedReadCoverageSpark.intervalTreeToGATKRead(barcode3, samHeader, tree3.iterator().next(), readMetadata);
        Assert.assertEquals(samRecord3.convertToSAMRecord(samHeader).getSAMString(), "AAAAAAA\t1\t1\t997\t60\t3M7M18N10M\t*\t0\t0\tACACACACACGTGTGTGTGT\tSSSSSSSSSSTTTTTTTTTT\tBX:Z:AAAAAAA\n");

        final SVIntervalTree<List<GATKRead>> tree4 = createTestSVIntervalTree4(samHeader, readMetadata);
        final String barcode4 = "CCCCCCC";

        final GATKRead samRecord4 = CollectLinkedReadCoverageSpark.intervalTreeToGATKRead(barcode4, samHeader, tree4.iterator().next(), readMetadata);
        Assert.assertEquals(samRecord4.convertToSAMRecord(samHeader).getSAMString(), "CCCCCCC\t1\t1\t997\t60\t3M7M2M4M\t*\t0\t0\tACACACACACGTGTGT\tSSSSSSSSSSTTTTTT\tBX:Z:CCCCCCC\n");

        final SVIntervalTree<List<GATKRead>> tree2 = createTestSVIntervalTree2(samHeader, readMetadata);
        final String barcode2 = "GGGGGGG";

        final GATKRead samRecord2 = CollectLinkedReadCoverageSpark.intervalTreeToGATKRead(barcode2, samHeader, tree2.iterator().next(), readMetadata);
        Assert.assertEquals(samRecord2.convertToSAMRecord(samHeader).getSAMString(), "GGGGGGG\t1\t1\t997\t60\t3M7M18D8M\t*\t0\t0\tACACACACACGTGTGTGT\tSSSSSSSSSSTTTTTTTT\tBX:Z:GGGGGGG\n");

        final SVIntervalTree<List<GATKRead>> tree = createTestSVIntervalTree1(samHeader, readMetadata);
        final String barcode = "ACTGACTG";

        final GATKRead samRecord = CollectLinkedReadCoverageSpark.intervalTreeToGATKRead(barcode, samHeader, tree.iterator().next(), readMetadata);
        Assert.assertEquals(samRecord.convertToSAMRecord(samHeader).getSAMString(), "ACTGACTG\t1\t1\t1000\t60\t10M10N6M1D4M4M10N5M1I4M96N10M\t*\t0\t0\tACACACACACGTGTGTGTGTAGAGGCGCGCGCGCTTTTTTTTTT\tSSSSSSSSSSTTTTTTTTTTUUUUVVVVVVVVVVWWWWWWWWWW\tBX:Z:ACTGACTG\n");

        final SVIntervalTree<List<GATKRead>> tree5 = createTestSVIntervalTree5(samHeader, readMetadata);
        final String barcode5 = "TTTTTTTT";

        final GATKRead samRecord5 = CollectLinkedReadCoverageSpark.intervalTreeToGATKRead(barcode5, samHeader, tree5.iterator().next(), readMetadata);
        Assert.assertEquals(samRecord5.convertToSAMRecord(samHeader).getSAMString(), "TTTTTTTT\t1\t1\t997\t60\t3M7M2M20D2M2M\t*\t0\t0\tACACACACACGTGTGT\tSSSSSSSSSSTTTTTT\tBX:Z:TTTTTTTT\n");

        final SVIntervalTree<List<GATKRead>> tree6 = createTestSVIntervalTree6(samHeader, readMetadata);
        final String barcode6 = "GGGGGGG";

        final GATKRead samRecord6 = CollectLinkedReadCoverageSpark.intervalTreeToGATKRead(barcode6, samHeader, tree6.iterator().next(), readMetadata);
        Assert.assertEquals(samRecord6.convertToSAMRecord(samHeader).getSAMString(),
                "GGGGGGG\t1\t1\t1000\t60\t150M1M6M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTGTGTG\tSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSTTTTTTT\tBX:Z:GGGGGGG\n");

    }

    private SVIntervalTree<List<GATKRead>> createTestSVIntervalTree1(final SAMFileHeader artificialSamHeader, final ReadMetadata readMetadata) {
        final SVIntervalTree<List<GATKRead>> tree = new SVIntervalTree<>();
        final List<GATKRead> resultList1 = new ArrayList<>();
        resultList1.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1000, "ACACACACAC".getBytes(), "2222222222".getBytes(), "10M"));
        resultList1.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1020, "GTGTGTGTGT".getBytes(), "3333333333".getBytes(), "6M1D4M"));
        resultList1.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1025, "AGAGAGAGAG".getBytes(), "4444444444".getBytes(), "10M"));
        resultList1.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1045, "GCGCGCGCGC".getBytes(), "5555555555".getBytes(), "5M1I4M"));
        resultList1.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1150, "TTTTTTTTTT".getBytes(), "6666666666".getBytes(), "10M"));
        tree.put(new SVInterval(readMetadata.getContigID("1"), 1000, 1160), resultList1);
        return tree;
    }

    private SVIntervalTree<List<GATKRead>> createTestSVIntervalTree4(final SAMFileHeader artificialSamHeader, final ReadMetadata readMetadata) {
        final SVIntervalTree<List<GATKRead>> tree = new SVIntervalTree<>();

        final List<GATKRead> resultList = new ArrayList<>();
        // reads are like:
        //     acaCACACAC
        //           gtgtgtGTGT
        resultList.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1000, "ACACACACAC".getBytes(), "2222222222".getBytes(), "3S7M"));
        resultList.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1009, "GTGTGTGTGT".getBytes(), "3333333333".getBytes(), "6S4M"));
        tree.put(new SVInterval(readMetadata.getContigID("1"),1000, 1013), resultList);
        return tree;
    }


    private SVIntervalTree<List<GATKRead>> createTestSVIntervalTree2(final SAMFileHeader artificialSamHeader, final ReadMetadata readMetadata) {
        final SVIntervalTree<List<GATKRead>> tree = new SVIntervalTree<>();

        final List<GATKRead> resultList = new ArrayList<>();
        // reads are like:
        //     acaCACACAC
        //           GT--------------------GTGTGTGT
        resultList.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1000, "ACACACACAC".getBytes(), "2222222222".getBytes(), "3S7M"));
        resultList.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1003, "GTGTGTGTGT".getBytes(), "3333333333".getBytes(), "2M20D8M"));
        tree.put(new SVInterval(readMetadata.getContigID("1"), 1000, 1033), resultList);
        return tree;
    }


    private SVIntervalTree<List<GATKRead>> createTestSVIntervalTree3(final SAMFileHeader artificialSamHeader, final ReadMetadata readMetadata) {
        final SVIntervalTree<List<GATKRead>> tree = new SVIntervalTree<>();

        final List<GATKRead> resultList = new ArrayList<>();
        // Simple non overlapping reads:
        //     acaCACACAC
        //                        GTGTGTGTGT
        resultList.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1000, "ACACACACAC".getBytes(), "2222222222".getBytes(), "3S7M"));
        resultList.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1025, "GTGTGTGTGT".getBytes(), "3333333333".getBytes(), "10M"));
        tree.put(new SVInterval(readMetadata.getContigID("1"), 1000, 1035), resultList);
        return tree;
    }

    private SVIntervalTree<List<GATKRead>> createTestSVIntervalTree5(final SAMFileHeader artificialSamHeader, final ReadMetadata readMetadata) {
        final SVIntervalTree<List<GATKRead>> tree = new SVIntervalTree<>();

        final List<GATKRead> resultList = new ArrayList<>();
        // reads are like:
        //     acaCACACAC
        //           GTGTGT--------------------GTgt
        resultList.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1000, "ACACACACAC".getBytes(), "2222222222".getBytes(), "3S7M"));
        resultList.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1003, "GTGTGTGTGT".getBytes(), "3333333333".getBytes(), "6M20D2M2S"));
        tree.put(new SVInterval(readMetadata.getContigID("1"), 1000, 1033), resultList);
        return tree;
    }

    private SVIntervalTree<List<GATKRead>> createTestSVIntervalTree6(final SAMFileHeader artificialSamHeader, final ReadMetadata readMetadata) {
        final SVIntervalTree<List<GATKRead>> tree = new SVIntervalTree<>();

        final List<GATKRead> resultList = new ArrayList<>();
        resultList.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1000,
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".getBytes(),
                "222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222".getBytes(),
                "150M"));
        resultList.add(ArtificialReadUtils.createArtificialRead(artificialSamHeader, "1", 0, 1002,
                "GTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG".getBytes(),
                "3333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333".getBytes(),
                "75M4D70M6S"));
        tree.put(new SVInterval(readMetadata.getContigID("1"), 1000, 1151), resultList);
        return tree;
    }

}