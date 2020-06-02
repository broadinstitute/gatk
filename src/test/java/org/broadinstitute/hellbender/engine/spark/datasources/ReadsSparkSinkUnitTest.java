package org.broadinstitute.hellbender.engine.spark.datasources;


import htsjdk.samtools.*;
import htsjdk.samtools.util.FileExtensions;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.hdfs.MiniDFSCluster;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.GATKPath;
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
                {testDataDir + "tools/BQSR/HiSeq.1mb.1RG.2k_lines.bam", "ReadsSparkSinkUnitTest1", null, ".bam", true, true, 100L},
                {testDataDir + "tools/BQSR/HiSeq.1mb.1RG.2k_lines.bam", "ReadsSparkSinkUnitTest1", null, ".bam", true, true, 1L}, // check SBI granularity setting
                {testDataDir + "tools/BQSR/HiSeq.1mb.1RG.2k_lines.bam", "ReadsSparkSinkUnitTest1", null, ".bam", true, false, 100L}, // write BAI, don't write SBI
                {testDataDir + "tools/BQSR/HiSeq.1mb.1RG.2k_lines.bam", "ReadsSparkSinkUnitTest1", null, ".bam", false, true, 100L}, // don't write BAI, write SBI
                {testDataDir + "tools/BQSR/HiSeq.1mb.1RG.2k_lines.bam", "ReadsSparkSinkUnitTest1", null, ".bam", false, false, 100L}, // don't write BAI, don't write SBI
                {testDataDir + "tools/BQSR/expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.bam", "ReadsSparkSinkUnitTest2", null, ".bam", true, true, 100L},

                // This file has unmapped reads that are set to the position of their mates -- the ordering check
                // in the tests below will fail if our ordering of these reads relative to the mapped reads
                // is not consistent with the definition of coordinate sorting as defined in
                // htsjdk.samtools.SAMRecordCoordinateComparator
                {testDataDir + "tools/BQSR/CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam", "ReadsSparkSinkUnitTest3", null, ".bam", true, true, 100L},
                {testDataDir + "tools/BQSR/NA12878.chr17_69k_70k.dictFix.cram", "ReadsSparkSinkUnitTest5",
                                                publicTestDir + "human_g1k_v37.chr17_1Mb.fasta", ".cram", true, true, 100L},

                {testDataDir + "tools/BQSR/HiSeq.1mb.1RG.2k_lines.bam", "ReadsSparkSinkUnitTest6", null, ".sam", true, true, 100L},
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

    // This bam was samtools sorted queryname bam, we expect if this were sorted to match the header that this would no longer match read-for-read due to differences in queryname-sort definitions compared to htsjdk
    @Test
    public void testReadsSparkSinkNotSortingReadsToHeader() throws IOException {
        final GATKPath inputBam = new GATKPath(testDataDir + "engine/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.10000000-10000020.with.unmapped.queryname.samtools.sam");
        final File outputFile = createTempFile("ReadsSparkSinkNotSorting", ".bam");
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        JavaRDD<GATKRead> rddParallelReads = readSource.getParallelReads(inputBam, null);
        SAMFileHeader header = readSource.getHeader(inputBam, null);

        ReadsSparkSink.writeReads(ctx, outputFile.getPath(), null, rddParallelReads, header, ReadsWriteFormat.SINGLE, 0, null, true, true, false, SBIIndexWriter.DEFAULT_GRANULARITY);

        JavaRDD<GATKRead> rddParallelReads2 = readSource.getParallelReads(new GATKPath(outputFile.getPath()), null);
        final List<GATKRead> writtenReads = rddParallelReads2.collect();

        JavaRDD<GATKRead> rddParallelReads3 = readSource.getParallelReads(inputBam, null);
        final List<GATKRead> inputReads = rddParallelReads3.collect();

        Assert.assertEquals(writtenReads.size(), inputReads.size());
        for (int i = 0; i < writtenReads.size(); i++) {
            Assert.assertEquals(writtenReads.get(i), inputReads.get(i), "These bams were likely out of order to eachother, which may have been caused by automatic sorting of the output bam");
        }
    }

    @Test(dataProvider = "loadReadsBAM", groups = "spark")
    public void readsSinkTest(String inputBam, String outputFileName, String referenceFile, String outputFileExtension, boolean writeBai, boolean writeSbi, long sbiGranularity) throws IOException {
        final File outputFile = createTempFile(outputFileName, outputFileExtension);
        assertSingleShardedWritingWorks(new GATKPath(inputBam), referenceFile, outputFile.getAbsolutePath(), null, writeBai, writeSbi, sbiGranularity);
    }

    @Test(dataProvider = "loadReadsBAM", groups = "spark")
    public void testSpecifyPartsDir(String inputBam, String outputFileName, String referenceFile, String outputFileExtension, boolean writeBai, boolean writeSbi, long sbiGranularity) throws IOException {
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
            assertSingleShardedWritingWorks(new GATKPath(inputBam), referenceFile, outputFile.getAbsolutePath(), nonDefaultShardsDir.getAbsolutePath(), writeBai, writeSbi, sbiGranularity);

            // Test that the file wasn't deleted when spark cleared its temp directory
            Assert.assertTrue(Files.exists(defaultPartsDir));

        } finally {
            // Remove the file this time
            Runtime.getRuntime().exec("rm -r " + defaultPartsDir );
        }
    }

    @Test(dataProvider = "loadReadsBAM", groups = "spark")
    public void readsSinkHDFSTest(String inputBam, String outputFileName, String referenceFileName, String outputFileExtension, boolean writeBai, boolean writeSbi, long sbiGranularity) throws IOException {
        final String outputHDFSPath = MiniClusterUtils.getTempPath(cluster, outputFileName, outputFileExtension).toString();
        Assert.assertTrue(BucketUtils.isHadoopUrl(outputHDFSPath));
        assertSingleShardedWritingWorks(new GATKPath(inputBam), referenceFileName, outputHDFSPath, null, writeBai, writeSbi, sbiGranularity);
    }

    @Test(dataProvider = "loadReadsBAM", groups = "spark")
    public void testWritingToAnExistingFileHDFS(String inputBam, String outputFileName, String referenceFileName, String outputFileExtension, boolean writeBai, boolean writeSbi, long sbiGranularity) throws IOException {
        final Path outputPath = MiniClusterUtils.getTempPath(cluster, outputFileName, outputFileExtension);
        final FileSystem fs = outputPath.getFileSystem(new Configuration());
        Assert.assertTrue(fs.createNewFile(outputPath));
        Assert.assertTrue(fs.exists(outputPath));
        assertSingleShardedWritingWorks(new GATKPath(inputBam), referenceFileName, outputPath.toString(), null, writeBai, writeSbi, sbiGranularity);
    }

    @Test(groups = "spark")
    public void testWritingToFileURL() throws IOException {
        GATKPath inputBam = new GATKPath(testDataDir + "tools/BQSR/HiSeq.1mb.1RG.2k_lines.bam");
        String outputUrl = "file:///" + createTempFile("ReadsSparkSinkUnitTest1", ".bam").getAbsolutePath();
        assertSingleShardedWritingWorks(inputBam, null, outputUrl, null, true, true, 100L);
    }

    private void assertSingleShardedWritingWorks(GATKPath inputBam, String referenceFile, String outputPath, String outputPartsPath, boolean writeBai, boolean writeSbi, long sbiGranularity) throws IOException {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        final GATKPath referencePath = referenceFile == null ? null : new GATKPath(referenceFile);

        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        JavaRDD<GATKRead> rddParallelReads = readSource.getParallelReads(inputBam, referencePath);
        SAMFileHeader header = readSource.getHeader(inputBam, referencePath);

        ReadsSparkSink.writeReads(ctx, outputPath, referencePath, rddParallelReads, header, ReadsWriteFormat.SINGLE, 0, outputPartsPath, writeBai, writeSbi, true, sbiGranularity);

        // check that a bai file is created
        if (new GATKPath(outputPath).isBam() && writeBai) {
            Assert.assertTrue(Files.exists(IOUtils.getPath(outputPath + FileExtensions.BAI_INDEX)));
        }
        // check that a splitting bai file is created with correct granularity
        if (new GATKPath(outputPath).isBam() && writeSbi) {
            final java.nio.file.Path sbiPath = IOUtils.getPath(outputPath + FileExtensions.SBI);
            Assert.assertTrue(Files.exists(sbiPath));
            final SBIIndex sbi = SBIIndex.load(sbiPath);
            Assert.assertEquals(sbi.getGranularity(), sbiGranularity);
        }

        JavaRDD<GATKRead> rddParallelReads2 = readSource.getParallelReads(new GATKPath(outputPath), referencePath);
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
    public void readsSinkShardedTest(String inputBam, String outputFileName, String referenceFile, String outputFileExtension, boolean writeBai, boolean writeSbi, long sbiGranularity) throws IOException {
        final GATKPath inputBamSpecifier = new GATKPath(inputBam);
        final File outputFile = createTempFile(outputFileName, outputFileExtension);
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        final GATKPath referencePath = referenceFile == null ? null : new GATKPath(referenceFile);

        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        JavaRDD<GATKRead> rddParallelReads = readSource.getParallelReads(inputBamSpecifier, referencePath);
        rddParallelReads = rddParallelReads.repartition(2); // ensure that the output is in two shards
        SAMFileHeader header = readSource.getHeader(inputBamSpecifier, referencePath);

        ReadsSparkSink.writeReads(ctx, outputFile.getAbsolutePath(), referencePath, rddParallelReads, header, ReadsWriteFormat.SHARDED, 0, null, false, sbiGranularity);
        int shards = outputFile.listFiles((dir, name) -> !name.startsWith(".") && !name.startsWith("_")).length;
        Assert.assertEquals(shards, 2);
        // check that no local .crc files are created
        int crcs = outputFile.listFiles((dir, name) -> name.startsWith(".") && name.endsWith(".crc")).length;
        Assert.assertEquals(crcs, 0);

        JavaRDD<GATKRead> rddParallelReads2 = readSource.getParallelReads(new GATKPath(outputFile.getAbsolutePath()), referencePath);
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
        final GATKPath inputBamSpecifier = new GATKPath(inputBam);

        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        JavaRDD<GATKRead> rddParallelReads = readSource.getParallelReads(inputBamSpecifier, null)
                .filter(r -> !r.isUnmapped()); // filter out unmapped reads (see comment below)
        SAMFileHeader header = readSource.getHeader(inputBamSpecifier, null);

        ReadsSparkSink.writeReads(ctx, outputDirectory.getAbsolutePath(), null, rddParallelReads, header, ReadsWriteFormat.ADAM, 0, null, true, SBIIndexWriter.DEFAULT_GRANULARITY);

        JavaRDD<GATKRead> rddParallelReads2 = readSource.getADAMReads(new GATKPath(outputDirectory.getAbsolutePath()), null, header);
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
