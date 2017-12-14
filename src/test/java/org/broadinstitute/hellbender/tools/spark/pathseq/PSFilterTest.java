package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndexCache;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.File;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

public class PSFilterTest extends CommandLineProgramTest {

    private static final String BWA_IMAGE_PATH = publicTestDir + "hg19mini.fasta.img";

    private static final GATKRead fastqRecordToGATKRead(final FastqRecord rec) {
        final GATKRead read = new SAMRecordToGATKReadAdapter(new SAMRecord(null));
        read.setName(rec.getReadName());
        read.setBaseQualities(rec.getBaseQualities());
        read.setBases(rec.getReadBases());
        return read;
    }

    private static final List<GATKRead> getReadsFromFastq(final File fastqFile) {
        try (final FastqReader reader = new FastqReader(fastqFile)) {
            final List<GATKRead> readList = new ArrayList<>();
            reader.forEachRemaining(read -> readList.add(fastqRecordToGATKRead(read)));
            return readList;
        }
    }

    private static final List<GATKRead> getPairedReadsFromFastq(final File fastqFile1, final File fastqFile2) {
        final List<GATKRead> readList1 = getReadsFromFastq(fastqFile1);
        final List<GATKRead> readList2 = getReadsFromFastq(fastqFile2);
        for (final GATKRead gatkRead : readList1) {
            gatkRead.setIsPaired(true);
            gatkRead.setIsFirstOfPair();
        }
        for (final GATKRead gatkRead : readList2) {
            gatkRead.setIsPaired(true);
            gatkRead.setIsSecondOfPair();
        }
        readList1.addAll(readList2);
        readList1.sort(Comparator.comparing(GATKRead::getName));
        return readList1;
    }

    private static final List<GATKRead> makeReadSet(final SAMFileHeader header) {
        final List<GATKRead> readList = new ArrayList<>(3);

        final List<GATKRead> readPair1 = ArtificialReadUtils.createPair(header, "paired_1", 101, 1, 102, false, false);
        readList.addAll(readPair1);

        final List<GATKRead> readPair2 = ArtificialReadUtils.createPair(header, "paired_2", 101, 1, 102, false, false);
        readList.addAll(readPair2);

        final GATKRead unpairedRead = ArtificialReadUtils.createRandomRead(101);
        unpairedRead.setName("unpaired_1");
        readList.add(unpairedRead);

        final List<GATKRead> readPair3 = ArtificialReadUtils.createPair(header, "paired_3", 101, 1, 102, false, false);
        final GATKRead formerlyPairedRead = readPair3.get(0);
        readList.add(formerlyPairedRead);

        return readList;
    }

    @Test(groups = "spark")
    public void testDoSetPairFlags() {

        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final SAMSequenceDictionary seq = new SAMSequenceDictionary();
        seq.addSequence(new SAMSequenceRecord("test_seq", 1000));
        final SAMFileHeader header = new SAMFileHeader(seq);

        final List<GATKRead> readList = makeReadSet(header);
        final JavaRDD<GATKRead> reads = ctx.parallelize(readList);
        ;

        final List<GATKRead> result = PSFilter.setPairFlags(reads, 100).collect();

        Assert.assertEquals(result.size(), 6);
        for (final GATKRead read : result) {
            if (read.getName().equals("paired_1") || read.getName().equals("paired_2")) {
                Assert.assertTrue(read.isPaired());
            } else {
                Assert.assertFalse(read.isPaired());
            }
        }

    }

    @Test
    public void testCanonicalizeRead() {

        final GATKRead read_1 = ArtificialReadUtils.createRandomRead(100);
        read_1.setBases(ArtificialReadUtils.createRandomReadBases(100, true));

        final GATKRead read_2 = read_1.copy();
        final byte[] bases_2 = read_1.getBases();
        SequenceUtil.reverseComplement(bases_2);
        read_2.setBases(bases_2);

        final GATKRead read_3 = read_1.copy();

        final GATKRead read_4 = ArtificialReadUtils.createRandomRead(100);
        read_4.setBases(ArtificialReadUtils.createRandomReadBases(100, true));

        final GATKRead read_5 = ArtificialReadUtils.createRandomRead(100);
        final byte[] bases_3 = read_4.getBases();
        switch (bases_3[bases_3.length - 1]) {
            case 'A':
                bases_3[bases_3.length - 1] = 'G';
                break;
            default:
                bases_3[bases_3.length - 1] = 'A';
                break;
        }
        read_5.setBases(bases_3);

        Tuple2<Long, GATKRead> result_1 = PSFilter.canonicalizeRead(read_1);
        Tuple2<Long, GATKRead> result_2 = PSFilter.canonicalizeRead(read_2);
        Tuple2<Long, GATKRead> result_3 = PSFilter.canonicalizeRead(read_3);
        Tuple2<Long, GATKRead> result_4 = PSFilter.canonicalizeRead(read_4);
        Tuple2<Long, GATKRead> result_5 = PSFilter.canonicalizeRead(read_5);
        Assert.assertEquals(result_1._1, result_2._1);
        Assert.assertEquals(result_1._1, result_3._1);
        Assert.assertNotEquals(result_1._1, result_4._1);
        Assert.assertNotEquals(result_4._1, result_5._1);

        Assert.assertEquals(read_1, result_1._2);
        Assert.assertEquals(read_2, result_2._2);
        Assert.assertEquals(read_3, result_3._2);
        Assert.assertEquals(read_4, result_4._2);
        Assert.assertEquals(read_5, result_5._2);
    }

    @DataProvider(name = "fastq")
    public Object[][] loadReads() {
        return new Object[][]{
                {"hg19mini.0.fq", null, 3},
                {"hg19mini.1.fq", "hg19mini.2.fq", 2},
                {"hg19mini.negative.0.fq", null, 45},
                {"hg19mini.negative.1.fq", "hg19mini.negative.2.fq", 60}
        };
    }

    @Test(dataProvider = "fastq", groups = "spark")
    public void testBwaFilter(final String reads1, final String reads2, final int expectedNum) {

        final List<GATKRead> readList;
        if (reads2 == null) {
            final File fastqFile = getTestFile(reads1);
            readList = getReadsFromFastq(fastqFile);
        } else {
            final File fastqFile1 = getTestFile(reads1);
            final File fastqFile2 = getTestFile(reads2);
            readList = getPairedReadsFromFastq(fastqFile1, fastqFile2);
        }

        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final JavaRDD<GATKRead> reads = ctx.parallelize(readList);

        final List<GATKRead> result = PSFilter.doBwaFilter(reads, BWA_IMAGE_PATH, 19, 1, 65).collect();
        BwaMemIndexCache.closeAllDistributedInstances(ctx);
        Assert.assertEquals(result.size(), expectedNum);
    }

    @Test(groups = "spark")
    public void testFilterDuplicateSequences() {

        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        final List<GATKRead> readList = new ArrayList<>(100000);
        final int numReads = 1000;
        for (int i = 0; i < numReads; i++) {
            final GATKRead read = ArtificialReadUtils.createRandomRead(100);
            read.setBases(ArtificialReadUtils.createRandomReadBases(100, true));
            read.setReadGroup("A");
            readList.add(read);
        }
        readList.get(2).setBases(readList.get(4).getBases());
        readList.get(2).setName("dupRead_1a");
        readList.get(4).setName("dupRead_1b");

        final byte[] bases_3 = readList.get(3).getBases();
        readList.get(9).setBases(bases_3);
        SequenceUtil.reverseComplement(bases_3);
        readList.get(1).setBases(bases_3);
        readList.get(1).setName("dupRead_2a");
        readList.get(3).setName("dupRead_2b");
        readList.get(9).setName("dupRead_2c");

        final JavaRDD<GATKRead> reads = ctx.parallelize(readList);
        final List<GATKRead> result = PSFilter.filterDuplicateSequences(reads).collect();

        Assert.assertEquals(result.size(), numReads - 3);
        int numDup_1 = 0;
        int numDup_2 = 0;
        for (final GATKRead read : result) {
            if (read.getName().startsWith("dupRead_1")) numDup_1++;
            else if (read.getName().startsWith("dupRead_2")) numDup_2++;
        }
        Assert.assertEquals(numDup_1, 1);
        Assert.assertEquals(numDup_2, 1);
    }

}