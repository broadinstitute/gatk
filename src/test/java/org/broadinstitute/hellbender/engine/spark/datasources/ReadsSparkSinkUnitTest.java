package org.broadinstitute.hellbender.engine.spark.datasources;


import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SBIIndex;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.hdfs.MiniDFSCluster;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.MiniClusterUtils;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

public class ReadsSparkSinkUnitTest extends GATKBaseTest {
    private MiniDFSCluster cluster;

    private static String testDataDir = publicTestDir + "org/broadinstitute/hellbender/";

    @BeforeClass(alwaysRun = true)
    private void setupMiniCluster() throws IOException {
        cluster = MiniClusterUtils.getMiniCluster();
    }

    @AfterClass(alwaysRun = true)
    private void shutdownMiniCluster() {
        MiniClusterUtils.stopCluster(cluster);
    }

    @DataProvider(name = "loadReadsBAM")
    public Object[][] loadReadsBAM() {
        return new Object[][]{
                {testDataDir + "tools/BQSR/HiSeq.1mb.1RG.2k_lines.bam", "ReadsSparkSinkUnitTest1", null, ".bam"},
                {testDataDir + "tools/BQSR/expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.bam", "ReadsSparkSinkUnitTest2", null, ".bam"},

                // This file has unmapped reads that are set to the position of their mates -- the ordering check
                // in the tests below will fail if our ordering of these reads relative to the mapped reads
                // is not consistent with the definition of coordinate sorting as defined in
                // htsjdk.samtools.SAMRecordCoordinateComparator
                {testDataDir + "tools/BQSR/CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam", "ReadsSparkSinkUnitTest3", null, ".bam"},
                {testDataDir + "tools/BQSR/NA12878.chr17_69k_70k.dictFix.cram", "ReadsSparkSinkUnitTest5",
                                                publicTestDir + "human_g1k_v37.chr17_1Mb.fasta", ".cram"},

                {testDataDir + "tools/BQSR/HiSeq.1mb.1RG.2k_lines.bam", "ReadsSparkSinkUnitTest6", null, ".sam"},
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
    public void readsSinkTest(String inputBam, String outputFileName, String referenceFile, String outputFileExtension) throws IOException {
        final File outputFile = createTempFile(outputFileName, outputFileExtension);
        assertSingleShardedWritingWorks(inputBam, referenceFile, outputFile.getAbsolutePath(), null);
    }

    @Test(dataProvider = "loadReadsBAM", groups = "spark")
    public void testSpecifyPartsDir(String inputBam, String outputFileName, String referenceFile, String outputFileExtension) throws IOException {
        final File outputFile = createTempFile(outputFileName, outputFileExtension);
        final File nonDefaultShardsDir = createTempDir(outputFileName + ".someOtherPlace");
        nonDefaultShardsDir.delete();

        final java.nio.file.Path defaultPartsDir = IOUtils.getPath(ReadsSparkSink.getDefaultPartsDirectory(outputFile.getAbsolutePath()));
        final java.nio.file.Path subpath = defaultPartsDir.resolve("subpath");

        try {
            // Make a directory with unusable permissions in place of where the default file will live
            Files.createDirectory(defaultPartsDir);
            Files.createFile(subpath);
            Runtime.getRuntime().exec("chmod a-w -R " + defaultPartsDir + "/");

            //show this succeeds when specifying a different path for the parts directory
            assertSingleShardedWritingWorks(inputBam, referenceFile, outputFile.getAbsolutePath(), nonDefaultShardsDir.getAbsolutePath());

            // Test that the file wasn't deleted when spark cleared its temp directory
            Assert.assertTrue(Files.exists(defaultPartsDir));

        } finally {
            // Remove the file this time
            Runtime.getRuntime().exec("rm -r " + defaultPartsDir );
        }
    }

    @Test(dataProvider = "loadReadsBAM", groups = "spark")
    public void readsSinkHDFSTest(String inputBam, String outputFileName, String referenceFileName, String outputFileExtension) throws IOException {
        final String outputHDFSPath = MiniClusterUtils.getTempPath(cluster, outputFileName, outputFileExtension).toString();
        Assert.assertTrue(BucketUtils.isHadoopUrl(outputHDFSPath));
        assertSingleShardedWritingWorks(inputBam, referenceFileName, outputHDFSPath, null);
    }

    @Test(dataProvider = "loadReadsBAM", groups = "spark")
    public void testWritingToAnExistingFileHDFS(String inputBam, String outputFileName, String referenceFileName, String outputFileExtension) throws IOException {
        final Path outputPath = MiniClusterUtils.getTempPath(cluster, outputFileName, outputFileExtension);
        final FileSystem fs = outputPath.getFileSystem(new Configuration());
        Assert.assertTrue(fs.createNewFile(outputPath));
        Assert.assertTrue(fs.exists(outputPath));
        assertSingleShardedWritingWorks(inputBam, referenceFileName, outputPath.toString(), null);
    }

    @Test(groups = "spark")
    public void testWritingToFileURL() throws IOException {
        String inputBam = testDataDir + "tools/BQSR/HiSeq.1mb.1RG.2k_lines.bam";
        String outputUrl = "file:///" + createTempFile("ReadsSparkSinkUnitTest1", ".bam").getAbsolutePath();
        assertSingleShardedWritingWorks(inputBam, null, outputUrl, null);
    }

    private void assertSingleShardedWritingWorks(String inputBam, String referenceFile, String outputPath, String outputPartsPath) throws IOException {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        JavaRDD<GATKRead> rddParallelReads = readSource.getParallelReads(inputBam, referenceFile);
        SAMFileHeader header = readSource.getHeader(inputBam, referenceFile);

        ReadsSparkSink.writeReads(ctx, outputPath, referenceFile, rddParallelReads, header, ReadsWriteFormat.SINGLE, 0, outputPartsPath);

        // check that a splitting bai file is created
        if (IOUtils.isBamFileName(outputPath)) {
            // Reinstate once https://github.com/disq-bio/disq/issues/45 is done
            //Assert.assertTrue(Files.exists(IOUtils.getPath(outputPath + SBIIndex.FILE_EXTENSION)));
        }

        JavaRDD<GATKRead> rddParallelReads2 = readSource.getParallelReads(outputPath, referenceFile);
        final List<GATKRead> writtenReads = rddParallelReads2.collect();

        assertReadsAreSorted(header, writtenReads);
        Assert.assertEquals(rddParallelReads.count(), rddParallelReads2.count());
    }

    private static void assertReadsAreSorted(SAMFileHeader header, List<GATKRead> writtenReads) {
        final SAMRecordCoordinateComparator comparator = new SAMRecordCoordinateComparator();
        // Assert that the reads are sorted.
        final int size = writtenReads.size();
        for (int i = 0; i < size-1; ++i) {
            final SAMRecord smaller = writtenReads.get(i).convertToSAMRecord(header);
            final SAMRecord larger = writtenReads.get(i + 1).convertToSAMRecord(header);
            final int compare = comparator.compare(smaller, larger);
            Assert.assertTrue(compare < 0, "Reads are out of order (compare=" + compare+"): " + smaller.getSAMString() + " and " + larger.getSAMString());
        }
    }

    @Test(dataProvider = "loadReadsBAM", groups = "spark")
    public void readsSinkShardedTest(String inputBam, String outputFileName, String referenceFile, String outputFileExtension) throws IOException {
        final File outputFile = createTempFile(outputFileName, outputFileExtension);
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        JavaRDD<GATKRead> rddParallelReads = readSource.getParallelReads(inputBam, referenceFile);
        rddParallelReads = rddParallelReads.repartition(2); // ensure that the output is in two shards
        SAMFileHeader header = readSource.getHeader(inputBam, referenceFile);

        ReadsSparkSink.writeReads(ctx, outputFile.getAbsolutePath(), referenceFile, rddParallelReads, header, ReadsWriteFormat.SHARDED);
        int shards = outputFile.listFiles((dir, name) -> !name.startsWith(".") && !name.startsWith("_")).length;
        Assert.assertEquals(shards, 2);
        // check that no local .crc files are created
        int crcs = outputFile.listFiles((dir, name) -> name.startsWith(".") && name.endsWith(".crc")).length;
        Assert.assertEquals(crcs, 0);

        JavaRDD<GATKRead> rddParallelReads2 = readSource.getParallelReads(outputFile.getAbsolutePath(), referenceFile);
        // reads are not globally sorted, so don't test that
        Assert.assertEquals(rddParallelReads.count(), rddParallelReads2.count());
    }

    @Test(enabled = false, dataProvider = "loadReadsADAM", groups = "spark")
    public void readsSinkADAMTest(String inputBam, String outputDirectoryName) throws IOException {
        // Since the test requires that we not create the actual output directory in advance,
        // we instead create its parent directory and mark it for deletion on exit. This protects
        // us from naming collisions across multiple instances of the test suite.
        final File outputParentDirectory = createTempDir(outputDirectoryName + "_parent");
        final File outputDirectory = new File(outputParentDirectory, outputDirectoryName);

        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        JavaRDD<GATKRead> rddParallelReads = readSource.getParallelReads(inputBam, null)
                .filter(r -> !r.isUnmapped()); // filter out unmapped reads (see comment below)
        SAMFileHeader header = readSource.getHeader(inputBam, null);

        ReadsSparkSink.writeReads(ctx, outputDirectory.getAbsolutePath(), null, rddParallelReads, header, ReadsWriteFormat.ADAM);

        JavaRDD<GATKRead> rddParallelReads2 = readSource.getADAMReads(outputDirectory.getAbsolutePath(), null, header);
        Assert.assertEquals(rddParallelReads.count(), rddParallelReads2.count());

        // Test the round trip
        List<GATKRead> samList = new ArrayList<>(rddParallelReads.collect());//make a mutable copy for sort
        List<GATKRead> adamList = new ArrayList<>(rddParallelReads2.collect());//make a mutable copy for sort
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
}
