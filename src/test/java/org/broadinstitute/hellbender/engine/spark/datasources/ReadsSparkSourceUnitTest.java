package org.broadinstitute.hellbender.engine.spark.datasources;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.hdfs.MiniDFSCluster;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.test.BaseTest;
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

                {dirBQSR + "HiSeq.1mb.1RG.2k_lines.alternate.bam",  null},
                {dirBQSR + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.bam", null},

                //reading CRAM on Spark needs a way to relax validation stringency https://github.com/broadinstitute/gatk/issues/1261
//                {NA12878_chr17_1k_CRAM, new File(v37_chr17_1Mb_Reference)},
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

    @Test(dataProvider = "loadReads", groups = "spark")
    public void readsSparkSourceTest(String bam, String referencePath) {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        JavaRDD<GATKRead> rddSerialReads = getSerialReads(ctx, bam, referencePath);
        JavaRDD<GATKRead> rddParallelReads = readSource.getParallelReads(bam, referencePath);

        List<GATKRead> serialReads = rddSerialReads.collect();
        List<GATKRead> parallelReads = rddParallelReads.collect();
        Assert.assertEquals(serialReads.size(), parallelReads.size());
    }

    @Test(dataProvider = "loadShardedReads", groups = "spark")
    public void shardedReadsSparkSourceTest(String expectedBam, String shardedBam, String referencePath) {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        JavaRDD<GATKRead> rddSerialReads = getSerialReads(ctx, expectedBam, referencePath);
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
    public void testReadFromFileAndHDFS() throws IOException {
        final File bam = new File( getToolTestDataDir(), "hdfs_file_test.bam");
        final File bai = new File( getToolTestDataDir(), "hdfs_file_test.bai");
        MiniDFSCluster cluster = null;
        try {
            cluster = new MiniDFSCluster.Builder(new Configuration()).build();
            final Path workingDirectory = cluster.getFileSystem().getWorkingDirectory();
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
        } finally {

            if (cluster != null) {
                cluster.shutdown();
            }
        }
    }

    /**
     * Loads Reads using samReaderFactory, then calling ctx.parallelize.
     * @param bam file to load
     * @return RDD of (SAMRecord-backed) GATKReads from the file.
     */
    public JavaRDD<GATKRead> getSerialReads(final JavaSparkContext ctx, final String bam, final String referencePath) {
        final SAMFileHeader readsHeader = ReadsSparkSource.getHeader(ctx, bam, null);
        List<SimpleInterval> intervals = IntervalUtils.getAllIntervalsForReference(readsHeader.getSequenceDictionary());

        final SamReaderFactory samReaderFactory;
        if (referencePath != null) {
            final File reference = new File(referencePath);
            samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).referenceSequence(reference);
        } else {
            samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
        }

        ReadsDataSource bam2 = new ReadsDataSource(new File(bam), samReaderFactory);
        bam2.setIntervalsForTraversal(intervals);
        List<GATKRead> records = Lists.newArrayList();
        for ( GATKRead read : bam2 ) {
            records.add(read);
        }
        return ctx.parallelize(records);
    }
}
