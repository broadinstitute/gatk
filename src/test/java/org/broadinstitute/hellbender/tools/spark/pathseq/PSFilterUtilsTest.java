package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.commons.io.IOUtils;
import org.apache.commons.io.LineIterator;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.TestException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

public class PSFilterUtilsTest extends CommandLineProgramTest {

    final static List<List<String>> readFastq(final File fastqFile) {
        try {
            final LineIterator itr = IOUtils.lineIterator(new FileReader(fastqFile));
            final List<List<String>> seqs = new ArrayList<>(100);
            while (itr.hasNext()) {
                final List<String> read = new ArrayList<>(3);
                String readName = itr.next();
                if (readName.endsWith("/1") || readName.endsWith("/2")) {
                    readName = readName.substring(0, readName.length() - 2);
                }
                read.add(readName);
                read.add(itr.next());
                itr.next();
                read.add(itr.next());
                seqs.add(read);
            }
            return seqs;
        } catch (IOException e) {
            logger.error("Could not open test file " + fastqFile.getAbsolutePath(), e);
        }
        throw new TestException("Invalid fastq test file " + fastqFile);
    }

    final static List<GATKRead> getReadsFromFastq(final File fastqFile) {
        final List<List<String>> fastqList = readFastq(fastqFile);
        return fastqList.stream().map(list -> {
            final String name = list.get(0);
            final byte[] bases = list.get(1).getBytes();
            final byte[] qual = list.get(2).getBytes();
            final GATKRead read = ArtificialReadUtils.createArtificialRead(bases, qual, "*");
            read.setName(name);
            return read;
        }).collect(Collectors.toList());
    }

    final static List<GATKRead> getPairedReadsFromFastq(final File fastqFile1, final File fastqFile2) {
        final List<GATKRead> readList1 = getReadsFromFastq(fastqFile1);
        final List<GATKRead> readList2 = getReadsFromFastq(fastqFile2);
        readList1.stream().forEach(read -> {
            read.setIsPaired(true);
            read.setIsFirstOfPair();
        });
        readList2.stream().forEach(read -> {
            read.setIsPaired(true);
            read.setIsSecondOfPair();
        });
        readList1.addAll(readList2);
        readList1.sort(Comparator.comparing(GATKRead::getName));
        return readList1;
    }

    final static List<GATKRead> makeReadSet(final SAMFileHeader header) {
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
        final JavaRDD<GATKRead> reads = ctx.parallelize(readList);;

        final List<GATKRead> result = PSFilterUtils.doSetPairFlags(reads).collect();

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

        Tuple2<Long, GATKRead> result_1 = PSFilterUtils.canonicalizeRead(read_1);
        Tuple2<Long, GATKRead> result_2 = PSFilterUtils.canonicalizeRead(read_2);
        Tuple2<Long, GATKRead> result_3 = PSFilterUtils.canonicalizeRead(read_3);
        Tuple2<Long, GATKRead> result_4 = PSFilterUtils.canonicalizeRead(read_4);
        Assert.assertEquals(result_1._1,result_2._1);
        Assert.assertEquals(result_1._1,result_3._1);
        Assert.assertNotEquals(result_1._1,result_4._1);
    }

    @DataProvider(name = "fastq")
    public Object[][] loadReads() {
        return new Object[][]{
                {"hg19mini.0.fq", null, 3},
                {"hg19mini.1.fq", "hg19mini.2.fq", 4},
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
            readList = getPairedReadsFromFastq(fastqFile1,fastqFile2);
        }

        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final JavaRDD<GATKRead> reads = ctx.parallelize(readList);

        final String indexPath = getTestFile("hg19mini.fasta.bwa_image").getPath();

        final List<GATKRead> result = PSFilterUtils.doBwaFilter(ctx, reads, hg19Header, indexPath, 19, 1, 70, 65).collect();
        Assert.assertEquals(result.size(), expectedNum);
    }
    @Test(groups = "spark")
    public void testFilterDuplicateSequences() {

        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        final List<GATKRead> readList = new ArrayList<>(100);
        for (int i = 0; i < 10; i++) {
            final GATKRead read = ArtificialReadUtils.createRandomRead(100);
            read.setBases(ArtificialReadUtils.createRandomReadBases(100, true));
            read.setReadGroup("A");
            readList.add(read);
        }
        readList.get(2).setBases(readList.get(4).getBases());
        readList.get(2).setName("dupRead_1a");
        readList.get(4).setName("dupRead_1b");

        readList.get(9).setBases(readList.get(3).getBases());
        readList.get(1).setBases(readList.get(3).getBases());
        readList.get(1).setName("dupRead_2a");
        readList.get(3).setName("dupRead_2b");
        readList.get(9).setName("dupRead_2c");

        final JavaRDD<GATKRead> reads = ctx.parallelize(readList);
        final List<GATKRead> result = PSFilterUtils.filterDuplicateSequences(reads).collect();

        Assert.assertEquals(result.size(), 7);
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