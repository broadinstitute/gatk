package org.broadinstitute.hellbender.engine.spark.datasources;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import htsjdk.samtools.*;
import org.apache.hadoop.fs.Path;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.GATKPathSpecifier;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.TraversalParameters;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadConstants;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.MiniClusterUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;

import static org.testng.Assert.assertEquals;

public class ReadsSparkSourceUnitTest extends GATKBaseTest {
    private static final String dirBQSR = toolsTestDir + "BQSR/";

    @DataProvider(name = "loadReads")
    public Object[][] loadReads() {
        return new Object[][]{
                {NA12878_chr17_1k_BAM, null},
                {dirBQSR + "HiSeq.1mb.1RG.2k_lines.alternate.bam", null},
                {dirBQSR + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.bam", null},
                {NA12878_chr17_1k_CRAM, v37_chr17_1Mb_Reference},
                {toolsTestDir + "valid.cram", toolsTestDir + "valid.fasta"}
        };
    }

    @DataProvider(name = "loadShardedReads")
    public Object[][] loadShardedReads() {
        String dir = dirBQSR;
        return new Object[][]{
                {dir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.bam", dir + "HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.sharded.bam", null},
        };
    }

    @DataProvider(name = "loadReadsValidation") // these require validation lenient or silent
    public Object[][] loadReadsValidation() {
        return new Object[][]{
                {NA12878_chr17_1k_BAM, null},
                {NA12878_chr17_1k_CRAM, v37_chr17_1Mb_Reference},
        };
    }

    @Test(dataProvider = "loadReads")
    public void testGetHeader(final String bam, final String referencePath) {
        doTestGetHeader(bam, referencePath);
    }

    @Test(groups="cloud")
    public void testGetHeaderOnGCSInputs() {
        doTestGetHeader(
                GCS_GATK_TEST_RESOURCES + "org/broadinstitute/hellbender/tools/count_reads_sorted.cram",
                GCS_GATK_TEST_RESOURCES + "org/broadinstitute/hellbender/tools/count_reads.fasta");
    }

    private void doTestGetHeader(final String bam, final String referencePath) {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final ReadsSparkSource readSource = new ReadsSparkSource(ctx, ReadConstants.DEFAULT_READ_VALIDATION_STRINGENCY);
        final SAMFileHeader samHeader = readSource.getHeader(new GATKPathSpecifier(bam), referencePath == null ? null : new GATKPathSpecifier(referencePath));
        Assert.assertNotNull(samHeader);
    }

    @Test
    public void testCheckCRAMReference() {
        // the "hdfs://" case is covered in testCRAMReferenceFromHDFS
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        ReadsSparkSource.checkCramReference(ctx, new GATKPathSpecifier(NA12878_chr17_1k_CRAM), new GATKPathSpecifier(v37_chr17_1Mb_Reference));
    }

    // this tests handling the case where a reference is specified but the file doesn't exist
    @Test(expectedExceptions = UserException.MissingReference.class)
    public void loadReadsNonExistentReference() {
        doLoadReads(toolsTestDir + "valid.cram",
                GATKBaseTest.getSafeNonExistentFile("NonExistentReference.fasta").getAbsolutePath(),
                ReadConstants.DEFAULT_READ_VALIDATION_STRINGENCY);
    }

    private void doLoadReadsTest(String bam, String referencePath) {
        doLoadReads(bam, referencePath, ReadConstants.DEFAULT_READ_VALIDATION_STRINGENCY);
    }

    private void doLoadReads(String bam, String referencePathName, ValidationStringency validationStringency) {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        final GATKPathSpecifier inputBamSpecifier = new GATKPathSpecifier(bam);
        final GATKPathSpecifier referenceInputPath = referencePathName == null ? null : new GATKPathSpecifier(referencePathName);
        ReadsSparkSource readSource = new ReadsSparkSource(ctx, validationStringency);
        JavaRDD<GATKRead> rddSerialReads = getSerialReads(ctx, bam, referenceInputPath, validationStringency);
        JavaRDD<GATKRead> rddParallelReads = readSource.getParallelReads(inputBamSpecifier, referenceInputPath);

        List<GATKRead> serialReads = rddSerialReads.collect();
        List<GATKRead> parallelReads = rddParallelReads.collect();
        Assert.assertEquals(serialReads.size(), parallelReads.size());
    }

    //TODO: this test is probably no longer needed since this is tested at the GATKPathSpecifier level
    @Test(expectedExceptions = UserException.class, expectedExceptionsMessageRegExp = ".*java.net.UnknownHostException: bogus.*")
    public void readsSparkSourceUnknownHostTest() {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        ReadsSparkSource readSource = new ReadsSparkSource(ctx, ReadConstants.DEFAULT_READ_VALIDATION_STRINGENCY);
        readSource.getParallelReads(new GATKPathSpecifier("hdfs://bogus/path.bam"), null);
    }

    @Test(dataProvider = "loadReads", groups = "spark")
    public void readsSparkSourceTest(String bam, String referencePath) {
        doLoadReadsTest(bam, referencePath);
    }

    @Test(dataProvider = "loadReadsValidation", groups = "spark", expectedExceptions = SAMFormatException.class)
    public void readsSparkSourceTestStrict(String bam, String referencePath) {
        doLoadReads(bam, referencePath, ValidationStringency.STRICT);
    }

    @Test(dataProvider = "loadReadsValidation", groups = "spark")
    public void readsSparkSourceTestLenient(String bam, String referencePath) {
        doLoadReads(bam, referencePath, ValidationStringency.LENIENT);
    }

    @Test(dataProvider = "loadReadsValidation", groups = "spark")
    public void readsSparkSourceTestSilent(String bam, String referencePath) {
        doLoadReads(bam, referencePath, ValidationStringency.SILENT);
    }

    @Test(dataProvider = "loadShardedReads", groups = "spark")
    public void shardedReadsSparkSourceTest(String expectedBam, String shardedBam, String referencePath) {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final GATKPathSpecifier referenceInputPath = referencePath == null ? null : new GATKPathSpecifier(referencePath);
        final GATKPathSpecifier shardedBamSpecifier = new GATKPathSpecifier(shardedBam);

        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        JavaRDD<GATKRead> rddSerialReads = getSerialReads(ctx, expectedBam, referenceInputPath, ReadConstants.DEFAULT_READ_VALIDATION_STRINGENCY);
        JavaRDD<GATKRead> rddParallelReads = readSource.getParallelReads(shardedBamSpecifier, referenceInputPath);

        List<GATKRead> serialReads = rddSerialReads.collect();
        List<GATKRead> parallelReads = rddParallelReads.collect();
        Assert.assertEquals(parallelReads.size(), serialReads.size());
    }

    @Test(groups = "spark")
    public void testHeadersAreStripped() {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        final List<GATKRead> reads = readSource.getParallelReads(new GATKPathSpecifier(dirBQSR + "HiSeq.1mb.1RG.2k_lines.alternate.bam"), null).collect();

        for ( final GATKRead read : reads ) {
            Assert.assertNull(((SAMRecordToGATKReadAdapter)read).getEncapsulatedSamRecord().getHeader(), "ReadSparkSource failed to null out header for read");
        }
    }

    @Test(groups = "spark", expectedExceptions=UserException.class)
    public void testReject2BitCRAMReference() {
        doLoadReadsTest(NA12878_chr17_1k_CRAM, b37_2bit_reference_20_21);
    }

    @Test(groups = "spark", expectedExceptions=UserException.class)
    public void testCRAMReferenceRequired() {
        doLoadReadsTest(NA12878_chr17_1k_CRAM, null);
    }

    @Test(groups = "spark")
    public void testPartitionSizing(){

        final GATKPathSpecifier bam = new GATKPathSpecifier(dirBQSR + "HiSeq.1mb.1RG.2k_lines.alternate.bam"); //file is ~220 kB
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        JavaRDD<GATKRead> allInOnePartition = readSource.getParallelReads(bam, null);
        JavaRDD<GATKRead> smallPartitions = readSource.getParallelReads(bam, null,  100 * 1024); // 100 kB
        Assert.assertEquals(allInOnePartition.partitions().size(), 1);
        Assert.assertEquals(smallPartitions.partitions().size(), 3);
    }

    @Test(groups = "spark")
    public void testReadFromFileAndHDFS() throws Exception {
        //TODO: add or change getTestFile overload with GATKPathSpecifier
        final GATKPathSpecifier bam = new GATKPathSpecifier(getTestFile("hdfs_file_test.bam").getAbsolutePath());
        final File bai = getTestFile("hdfs_file_test.bai");
        MiniClusterUtils.runOnIsolatedMiniCluster( cluster -> {
            final Path workingDirectory = MiniClusterUtils.getWorkingDir(cluster);
            final Path bamPath = new Path(workingDirectory,"hdfs.bam");
            final Path baiPath = new Path(workingDirectory, "hdfs.bai");
            cluster.getFileSystem().copyFromLocalFile(new Path(bam.getURI()), bamPath);
            cluster.getFileSystem().copyFromLocalFile(new Path(bai.toURI()), baiPath);

            final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
            final ReadsSparkSource readsSparkSource = new ReadsSparkSource(ctx);
            final List<GATKRead> localReads = readsSparkSource.getParallelReads(bam, null).collect();
            final List<GATKRead> hdfsReads = readsSparkSource.getParallelReads(new GATKPathSpecifier(bamPath.toUri().toString()), null).collect();

            Assert.assertFalse(localReads.isEmpty());
            Assert.assertEquals(localReads, hdfsReads);
        });
    }

    @Test(groups = "spark")
    public void testCRAMReferenceFromHDFS() throws Exception {
        final GATKPathSpecifier cramSpecifier = new GATKPathSpecifier(NA12878_chr17_1k_CRAM);
        final GATKPathSpecifier reference = new GATKPathSpecifier(v37_chr17_1Mb_Reference);
        final GATKPathSpecifier referenceIndex = new GATKPathSpecifier(v37_chr17_1Mb_Reference + ".fai");

        MiniClusterUtils.runOnIsolatedMiniCluster( cluster -> {
            final Path workingDirectory = MiniClusterUtils.getWorkingDir(cluster);
            final Path cramHDFSPath = new Path(workingDirectory, "hdfs.cram");
            final Path refHDFSPath = new Path(workingDirectory, "hdfs.fasta");
            final Path refIndexHDFSPath = new Path(workingDirectory, "hdfs.fasta.fai");
            cluster.getFileSystem().copyFromLocalFile(new Path(cramSpecifier.getURI()), cramHDFSPath);
            cluster.getFileSystem().copyFromLocalFile(new Path(reference.getURI()), refHDFSPath);
            cluster.getFileSystem().copyFromLocalFile(new Path(referenceIndex.getURI()), refIndexHDFSPath);

            final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
            final ReadsSparkSource readsSparkSource = new ReadsSparkSource(ctx);
            final SAMFileHeader samHeader = readsSparkSource.getHeader(new GATKPathSpecifier(cramHDFSPath.toString()), new GATKPathSpecifier(refHDFSPath.toString()));
            Assert.assertNotNull(samHeader);
            final List<GATKRead> localReads = readsSparkSource.getParallelReads(cramSpecifier, reference).collect();
            final List<GATKRead> hdfsReads = readsSparkSource.getParallelReads(new GATKPathSpecifier(cramHDFSPath.toUri().toString()), new GATKPathSpecifier(refHDFSPath.toUri().toString())).collect();

            Assert.assertFalse(localReads.isEmpty());
            Assert.assertEquals(localReads, hdfsReads);
        });
    }

    @Test(groups = "spark")
    public void testIntervals() throws IOException {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        List<SimpleInterval> intervals =
                ImmutableList.of(new SimpleInterval("17", 69010, 69040), new SimpleInterval("17", 69910, 69920));
        TraversalParameters traversalParameters = new TraversalParameters(intervals, false);
        JavaRDD<GATKRead> reads = readSource.getParallelReads(new GATKPathSpecifier(NA12878_chr17_1k_BAM), null, traversalParameters);

        SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
        try (SamReader samReader = samReaderFactory.open(new File(NA12878_chr17_1k_BAM))) {
            int seqIndex = samReader.getFileHeader().getSequenceIndex("17");
            SAMRecordIterator query = samReader.query(new QueryInterval[]{new QueryInterval(seqIndex, 69010, 69040), new QueryInterval(seqIndex, 69910, 69920)}, false);
            Assert.assertEquals(reads.count(), Iterators.size(query));
        }
    }

    @Test(groups = "spark")
    public void testIntervalsWithUnmapped() throws IOException {
        final GATKPathSpecifier bamPathSpecifier = new GATKPathSpecifier(publicTestDir + "org/broadinstitute/hellbender/engine/CEUTrio.HiSeq.WGS.b37.NA12878.snippet_with_unmapped.bam");
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        List<SimpleInterval> intervals = ImmutableList.of(new SimpleInterval("20", 10000009, 10000011));
        TraversalParameters traversalParameters = new TraversalParameters(intervals, true);
        JavaRDD<GATKRead> reads = readSource.getParallelReads(bamPathSpecifier, null, traversalParameters);

        SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
        try (SamReader samReader = samReaderFactory.open(bamPathSpecifier.toPath())) {
            int seqIndex = samReader.getFileHeader().getSequenceIndex("20");
            SAMRecordIterator query = samReader.query(new QueryInterval[]{new QueryInterval(seqIndex, 10000009, 10000011)}, false);
            int queryReads = Iterators.size(query);
            query.close();
            SAMRecordIterator queryUnmapped = samReader.queryUnmapped();
            queryReads += Iterators.size(queryUnmapped);
            Assert.assertEquals(reads.count(), queryReads);
        }
    }

    /**
     * Loads Reads using samReaderFactory, then calling ctx.parallelize.
     * @param bam file to load
     * @return RDD of (SAMRecord-backed) GATKReads from the file.
     */
    public JavaRDD<GATKRead> getSerialReads(final JavaSparkContext ctx, final String bam, final GATKPathSpecifier referencePath, final ValidationStringency validationStringency) {
        final SAMFileHeader readsHeader = new ReadsSparkSource(ctx, validationStringency).getHeader(new GATKPathSpecifier(bam), referencePath);

        final SamReaderFactory samReaderFactory;
        if (referencePath != null) {
            samReaderFactory = SamReaderFactory.makeDefault().validationStringency(validationStringency).referenceSequence(referencePath.toPath());
        } else {
            samReaderFactory = SamReaderFactory.makeDefault().validationStringency(validationStringency);
        }

        ReadsDataSource bam2 = new ReadsDataSource(IOUtils.getPath(bam), samReaderFactory);
        List<GATKRead> records = Lists.newArrayList();
        for ( GATKRead read : bam2 ) {
            records.add(read);
        }
        return ctx.parallelize(records);
    }
}
