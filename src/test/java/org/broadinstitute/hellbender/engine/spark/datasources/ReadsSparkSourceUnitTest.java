package org.broadinstitute.hellbender.engine.spark.datasources;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import htsjdk.samtools.*;
import org.apache.hadoop.fs.Path;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadConstants;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.MiniClusterUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;

public class ReadsSparkSourceUnitTest extends BaseTest {

    private static final String dir = "src/test/resources/org/broadinstitute/hellbender/tools/";
    private static final String dirBQSR = dir + "BQSR/";


    @DataProvider(name = "loadReads")
    public Object[][] loadReads() {
        return new Object[][]{
                {NA12878_chr17_1k_BAM, null},
                {dirBQSR + "HiSeq.1mb.1RG.2k_lines.alternate.bam", null},
                {dirBQSR + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.bam", null},
                {NA12878_chr17_1k_CRAM, v37_chr17_1Mb_Reference},
                {dir + "valid.cram", dir + "valid.fasta"}
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

    private void doLoadReadsTest(String bam, String referencePath) {
        doLoadReads(bam, referencePath, ReadConstants.DEFAULT_READ_VALIDATION_STRINGENCY);
    }

    private void doLoadReads(String bam, String referencePath, ValidationStringency validationStringency) {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        ReadsSparkSource readSource = new ReadsSparkSource(ctx, validationStringency);
        JavaRDD<GATKRead> rddSerialReads = getSerialReads(ctx, bam, referencePath, validationStringency);
        JavaRDD<GATKRead> rddParallelReads = readSource.getParallelReads(bam, referencePath);

        List<GATKRead> serialReads = rddSerialReads.collect();
        List<GATKRead> parallelReads = rddParallelReads.collect();
        Assert.assertEquals(serialReads.size(), parallelReads.size());
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

        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        JavaRDD<GATKRead> rddSerialReads = getSerialReads(ctx, expectedBam, referencePath, ReadConstants.DEFAULT_READ_VALIDATION_STRINGENCY);
        JavaRDD<GATKRead> rddParallelReads = readSource.getParallelReads(shardedBam, referencePath);

        List<GATKRead> serialReads = rddSerialReads.collect();
        List<GATKRead> parallelReads = rddParallelReads.collect();
        Assert.assertEquals(serialReads.size(), parallelReads.size());
    }

    @Test(groups = "spark")
    public void testHeadersAreStripped() {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        final List<GATKRead> reads = readSource.getParallelReads(dirBQSR + "HiSeq.1mb.1RG.2k_lines.alternate.bam", null).collect();

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

    @Test
    public void testPartitionSizing(){

        String bam = dirBQSR + "HiSeq.1mb.1RG.2k_lines.alternate.bam"; //file is ~220 kB
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        JavaRDD<GATKRead> allInOnePartition = readSource.getParallelReads(bam, null);
        JavaRDD<GATKRead> smallPartitions = readSource.getParallelReads(bam, null,  100 * 1024); // 100 kB
        Assert.assertEquals(allInOnePartition.partitions().size(), 1);
        Assert.assertEquals(smallPartitions.partitions().size(), 2);
    }

    @Test
    public void testReadFromFileAndHDFS() throws Exception {
        final File bam = getTestFile("hdfs_file_test.bam");
        final File bai = getTestFile("hdfs_file_test.bai");
        MiniClusterUtils.runOnIsolatedMiniCluster( cluster -> {
            final Path workingDirectory = MiniClusterUtils.getWorkingDir(cluster);
            final Path bamPath = new Path(workingDirectory,"hdfs.bam");
            final Path baiPath = new Path(workingDirectory, "hdfs.bai");
            cluster.getFileSystem().copyFromLocalFile(new Path(bam.toURI()), bamPath);
            cluster.getFileSystem().copyFromLocalFile(new Path(bai.toURI()), baiPath);

            final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
            final ReadsSparkSource readsSparkSource = new ReadsSparkSource(ctx);
            final List<GATKRead> localReads = readsSparkSource.getParallelReads(bam.toURI().toString(), null).collect();
            final List<GATKRead> hdfsReads = readsSparkSource.getParallelReads(bamPath.toUri().toString(), null).collect();

            Assert.assertFalse(localReads.isEmpty());
            Assert.assertEquals(localReads, hdfsReads);
        });
    }

    @Test
    public void testCRAMReferenceFromHDFS() throws Exception {
        final File cram = new File(NA12878_chr17_1k_CRAM);
        final File reference = new File(v37_chr17_1Mb_Reference);
        final File referenceIndex = new File(v37_chr17_1Mb_Reference + ".fai");

        MiniClusterUtils.runOnIsolatedMiniCluster( cluster -> {
            final Path workingDirectory = MiniClusterUtils.getWorkingDir(cluster);
            final Path cramHDFSPath = new Path(workingDirectory, "hdfs.cram");
            final Path refHDFSPath = new Path(workingDirectory, "hdfs.fasta");
            final Path refIndexHDFSPath = new Path(workingDirectory, "hdfs.fasta.fai");
            cluster.getFileSystem().copyFromLocalFile(new Path(cram.toURI()), cramHDFSPath);
            cluster.getFileSystem().copyFromLocalFile(new Path(reference.toURI()), refHDFSPath);
            cluster.getFileSystem().copyFromLocalFile(new Path(referenceIndex.toURI()), refIndexHDFSPath);

            final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
            final ReadsSparkSource readsSparkSource = new ReadsSparkSource(ctx);
            final List<GATKRead> localReads = readsSparkSource.getParallelReads(cram.toURI().toString(), reference.toURI().toString()).collect();
            final List<GATKRead> hdfsReads = readsSparkSource.getParallelReads(cramHDFSPath.toUri().toString(), refHDFSPath.toUri().toString()).collect();

            Assert.assertFalse(localReads.isEmpty());
            Assert.assertEquals(localReads, hdfsReads);
        });
    }

    @Test
    public void testIntervals() throws IOException {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        List<SimpleInterval> intervals =
                ImmutableList.of(new SimpleInterval("17", 69010, 69040), new SimpleInterval("17", 69910, 69920));
        JavaRDD<GATKRead> reads = readSource.getParallelReads(NA12878_chr17_1k_BAM, null, intervals);

        SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
        try (SamReader samReader = samReaderFactory.open(new File(NA12878_chr17_1k_BAM))) {
            int seqIndex = samReader.getFileHeader().getSequenceIndex("17");
            SAMRecordIterator query = samReader.query(new QueryInterval[]{new QueryInterval(seqIndex, 69010, 69040), new QueryInterval(seqIndex, 69910, 69920)}, false);
            Assert.assertEquals(reads.count(), Iterators.size(query));
        }
    }

    /**
     * Loads Reads using samReaderFactory, then calling ctx.parallelize.
     * @param bam file to load
     * @return RDD of (SAMRecord-backed) GATKReads from the file.
     */
    public JavaRDD<GATKRead> getSerialReads(final JavaSparkContext ctx, final String bam, final String referencePath, final ValidationStringency validationStringency) {
        final SAMFileHeader readsHeader = new ReadsSparkSource(ctx, validationStringency).getHeader(bam, referencePath, null);
        List<SimpleInterval> intervals = IntervalUtils.getAllIntervalsForReference(readsHeader.getSequenceDictionary());

        final SamReaderFactory samReaderFactory;
        if (referencePath != null) {
            final File reference = new File(referencePath);
            samReaderFactory = SamReaderFactory.makeDefault().validationStringency(validationStringency).referenceSequence(reference);
        } else {
            samReaderFactory = SamReaderFactory.makeDefault().validationStringency(validationStringency);
        }

        ReadsDataSource bam2 = new ReadsDataSource(new File(bam), samReaderFactory);
        bam2.setTraversalBounds(intervals);
        List<GATKRead> records = Lists.newArrayList();
        for ( GATKRead read : bam2 ) {
            records.add(read);
        }
        return ctx.parallelize(records);
    }
}
