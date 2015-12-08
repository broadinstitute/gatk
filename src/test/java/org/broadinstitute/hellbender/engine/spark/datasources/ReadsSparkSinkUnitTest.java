package org.broadinstitute.hellbender.engine.spark.datasources;


import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

public class ReadsSparkSinkUnitTest extends BaseTest {

    private static String testDataDir = publicTestDir + "org/broadinstitute/hellbender/";

    @DataProvider(name = "loadReadsBAM")
    public Object[][] loadReadsBAM() {
        return new Object[][]{
                {testDataDir + "tools/BQSR/HiSeq.1mb.1RG.2k_lines.bam", "ReadsSparkSinkUnitTest1", ".bam"},
                {testDataDir + "tools/BQSR/expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.bam", "ReadsSparkSinkUnitTest2", ".bam"},

                // This file has unmapped reads that are set to the position of their mates -- the ordering check
                // in the tests below will fail if our ordering of these reads relative to the mapped reads
                // is not consistent with the definition of coordinate sorting as defined in
                // htsjdk.samtools.SAMRecordCoordinateComparator
                {testDataDir + "tools/BQSR/CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam", "ReadsSparkSinkUnitTest3", ".bam"},

                {NA12878_20_21_WGS_bam , "ReadsSparkSinkUnitTest4", ".bam"}
        };
    }

    @DataProvider(name = "loadReadsADAM")
    public Object[][] loadReadsADAM() {
        return new Object[][]{
                {testDataDir + "tools/BQSR/HiSeq.1mb.1RG.2k_lines.bam", "ReadsSparkSinkUnitTest1_ADAM"},
                {testDataDir + "tools/BQSR/expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.bam", "ReadsSparkSinkUnitTest2_ADAM"},

                // This file has unmapped reads that are set to the position of their mates -- the ordering check
                // in the tests below will fail if our ordering of these reads relative to the mapped reads
                // is not consistent with the definition of coordinate sorting as defined in
                // htsjdk.samtools.SAMRecordCoordinateComparator
                //
                // This test is currently disabled, as this test case doesn't pass for ADAM (we have an open ticket for this:
                // https://github.com/broadinstitute/gatk/issues/1254)
                // {testDataDir + "tools/BQSR/CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam", "ReadsSparkSinkUnitTest3_ADAM"},

                //This test is disabled because it fails on travis (passes locally though)
                //https://github.com/broadinstitute/gatk/issues/1254
                // {NA12878_chr17_1k_BAM, "ReadsSparkSinkUnitTest4_ADAM"}
        };
    }

    @Test(dataProvider = "loadReadsBAM", groups = "spark")
    public void readsSinkTest(String inputBam, String outputFileName, String outputFileExtension) throws IOException {
        final File outputFile = createTempFile(outputFileName, outputFileExtension);
        outputFile.deleteOnExit();
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        JavaRDD<GATKRead> rddParallelReads = readSource.getParallelReads(inputBam, null);
        SAMFileHeader header = ReadsSparkSource.getHeader(ctx, inputBam, null);

        ReadsSparkSink.writeReads(ctx, outputFile.getAbsolutePath(), rddParallelReads, header, ReadsWriteFormat.SINGLE);

        JavaRDD<GATKRead> rddParallelReads2 = readSource.getParallelReads(outputFile.getAbsolutePath(), null);
        final List<GATKRead> writtenReads = rddParallelReads2.collect();

        final SAMRecordCoordinateComparator comparator = new SAMRecordCoordinateComparator();
        // Assert that the reads are sorted.
        final int size = writtenReads.size();
        for (int i = 0; i < size-1; ++i) {
            final SAMRecord smaller = writtenReads.get(i).convertToSAMRecord(header);
            final SAMRecord larger = writtenReads.get(i + 1).convertToSAMRecord(header);
            final int compare = comparator.compare(smaller, larger);
            Assert.assertTrue(compare < 0, "Reads are out of order (compare=" + compare+"): " + smaller.getSAMString() + " and " + larger.getSAMString());
        }
        Assert.assertEquals(rddParallelReads.count(), rddParallelReads2.count());
    }

    @Test(dataProvider = "loadReadsBAM", groups = "spark")
    public void readsSinkShardedTest(String inputBam, String outputFileName, String outputFileExtension) throws IOException {
        final File outputFile = createTempFile(outputFileName, outputFileExtension);
        outputFile.deleteOnExit();
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        JavaRDD<GATKRead> rddParallelReads = readSource.getParallelReads(inputBam, null);
        rddParallelReads = rddParallelReads.repartition(2); // ensure that the output is in two shards
        SAMFileHeader header = ReadsSparkSource.getHeader(ctx, inputBam, null);

        ReadsSparkSink.writeReads(ctx, outputFile.getAbsolutePath(), rddParallelReads, header, ReadsWriteFormat.SHARDED);
        int shards = outputFile.listFiles((dir, name) -> !name.startsWith(".") && !name.startsWith("_")).length;
        Assert.assertEquals(shards, 2);

        JavaRDD<GATKRead> rddParallelReads2 = readSource.getParallelReads(outputFile.getAbsolutePath(), null);
        // reads are not globally sorted, so don't test that
        Assert.assertEquals(rddParallelReads.count(), rddParallelReads2.count());
    }

    @Test(dataProvider = "loadReadsADAM", groups = "spark")
    public void readsSinkADAMTest(String inputBam, String outputDirectoryName) throws IOException {
        // Since the test requires that we not create the actual output directory in advance,
        // we instead create its parent directory and mark it for deletion on exit. This protects
        // us from naming collisions across multiple instances of the test suite.
        final File outputParentDirectory = Files.createTempDirectory(outputDirectoryName + "_parent").toFile();
        IOUtils.deleteRecursivelyOnExit(outputParentDirectory);
        final File outputDirectory = new File(outputParentDirectory, outputDirectoryName);

        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        JavaRDD<GATKRead> rddParallelReads = readSource.getParallelReads(inputBam, null);
        SAMFileHeader header = ReadsSparkSource.getHeader(ctx, inputBam, null);

        ReadsSparkSink.writeReads(ctx, outputDirectory.getAbsolutePath(), rddParallelReads, header, ReadsWriteFormat.ADAM);

        JavaRDD<GATKRead> rddParallelReads2 = readSource.getADAMReads(outputDirectory.getAbsolutePath(), null, header);
        Assert.assertEquals(rddParallelReads.count(), rddParallelReads2.count());

        // Test the round trip
        List<GATKRead> samList = rddParallelReads.collect();
        List<GATKRead> adamList = rddParallelReads2.collect();
        Comparator<GATKRead> comparator = new ReadCoordinateComparator(header);
        samList.sort(comparator);
        adamList.sort(comparator);
        for (int i = 0; i < samList.size(); i++) {
            SAMRecord expected = samList.get(i).convertToSAMRecord(header);
            SAMRecord observed = adamList.get(i).convertToSAMRecord(header);
            // manually test equality of some fields, as there are issues with roundtrip BAM -> ADAM -> BAM
            // see https://github.com/bigdatagenomics/adam/issues/823
            Assert.assertEquals(observed.getReadName(), expected.getReadName(), "readname");
            Assert.assertEquals(observed.getAlignmentStart(), expected.getAlignmentStart(), "getAlignmentStart");
            Assert.assertEquals(observed.getAlignmentEnd(), expected.getAlignmentEnd(), "getAlignmentEnd");
            Assert.assertEquals(observed.getFlags(), expected.getFlags(), "getFlags");
            Assert.assertEquals(observed.getMappingQuality(), expected.getMappingQuality(), "getMappingQuality");
            Assert.assertEquals(observed.getMateAlignmentStart(), expected.getMateAlignmentStart(), "getMateAlignmentStart");
            Assert.assertEquals(observed.getCigar(), expected.getCigar(), "getCigar");
        }
    }

    @Test
    public void testGetBamFragments() throws IOException {
        final Path fragmentDir = new Path(publicTestDir + "org/broadinstitute/hellbender/engine/spark/ReadSparkSink_fragments_test");
        final FileSystem fs = fragmentDir.getFileSystem(new Configuration());

        final FileStatus[] bamFragments = ReadsSparkSink.getBamFragments(fragmentDir, fs);
        final List<String> expectedFragmentNames = Arrays.asList("part-r-00000", "part-r-00001", "part-r-00002", "part-r-00003");

        Assert.assertEquals(bamFragments.length, expectedFragmentNames.size(), "Wrong number of fragments returned by ReadsSparkSink.getBamFragments()");
        for ( int i = 0; i < bamFragments.length; ++i ) {
            Assert.assertEquals(bamFragments[i].getPath().getName(), expectedFragmentNames.get(i), "Fragments are not in correct order");
        }
    }
}
