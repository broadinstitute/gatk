package org.broadinstitute.hellbender.engine.spark.datasources;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import htsjdk.samtools.*;
import org.apache.hadoop.fs.Path;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.TraversalParameters;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaSparkEngine;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import static org.testng.Assert.assertEquals;

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

    // this tests handling the case where a reference is specified but the file doesn't exist
    @Test(expectedExceptions = UserException.MissingReference.class)
    public void loadReadsNonExistentReference() {
        doLoadReads(dir + "valid.cram",
                BaseTest.getSafeNonExistentFile("NonExistentReference.fasta").getAbsolutePath(),
                ReadConstants.DEFAULT_READ_VALIDATION_STRINGENCY);
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

    @Test(expectedExceptions = UserException.class, expectedExceptionsMessageRegExp = ".*Failed to read bam header from hdfs://bogus/path.bam.*")
    public void readsSparkSourceUnknownHostTest() {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        ReadsSparkSource readSource = new ReadsSparkSource(ctx, ReadConstants.DEFAULT_READ_VALIDATION_STRINGENCY);
        readSource.getParallelReads("hdfs://bogus/path.bam", null);
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
        Assert.assertEquals(parallelReads.size(), serialReads.size());
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

    @Test(groups = "spark")
    public void testPartitionSizing(){

        String bam = dirBQSR + "HiSeq.1mb.1RG.2k_lines.alternate.bam"; //file is ~220 kB
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        JavaRDD<GATKRead> allInOnePartition = readSource.getParallelReads(bam, null);
        JavaRDD<GATKRead> smallPartitions = readSource.getParallelReads(bam, null,  100 * 1024); // 100 kB
        Assert.assertEquals(allInOnePartition.partitions().size(), 1);
        Assert.assertEquals(smallPartitions.partitions().size(), 2);
    }

    @Test(groups = "spark")
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

    @Test(groups = "spark")
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

    @Test(groups = "spark")
    public void testIntervals() throws IOException {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        List<SimpleInterval> intervals =
                ImmutableList.of(new SimpleInterval("17", 69010, 69040), new SimpleInterval("17", 69910, 69920));
        TraversalParameters traversalParameters = new TraversalParameters(intervals, false);
        JavaRDD<GATKRead> reads = readSource.getParallelReads(NA12878_chr17_1k_BAM, null, traversalParameters);

        SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
        try (SamReader samReader = samReaderFactory.open(new File(NA12878_chr17_1k_BAM))) {
            int seqIndex = samReader.getFileHeader().getSequenceIndex("17");
            SAMRecordIterator query = samReader.query(new QueryInterval[]{new QueryInterval(seqIndex, 69010, 69040), new QueryInterval(seqIndex, 69910, 69920)}, false);
            Assert.assertEquals(reads.count(), Iterators.size(query));
        }
    }

    @Test(groups = "spark")
    public void testIntervalsWithUnmapped() throws IOException {
        String bam = publicTestDir + "org/broadinstitute/hellbender/engine/CEUTrio.HiSeq.WGS.b37.NA12878.snippet_with_unmapped.bam";
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        List<SimpleInterval> intervals = ImmutableList.of(new SimpleInterval("20", 10000009, 10000011));
        TraversalParameters traversalParameters = new TraversalParameters(intervals, true);
        JavaRDD<GATKRead> reads = readSource.getParallelReads(bam, null, traversalParameters);

        SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
        try (SamReader samReader = samReaderFactory.open(new File(bam))) {
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
    public JavaRDD<GATKRead> getSerialReads(final JavaSparkContext ctx, final String bam, final String referencePath, final ValidationStringency validationStringency) {
        final SAMFileHeader readsHeader = new ReadsSparkSource(ctx, validationStringency).getHeader(bam, referencePath);

        final SamReaderFactory samReaderFactory;
        if (referencePath != null) {
            final File reference = new File(referencePath);
            samReaderFactory = SamReaderFactory.makeDefault().validationStringency(validationStringency).referenceSequence(reference);
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

    @DataProvider(name="readPairsAndPartitions")
    public Object[][] readPairsAndPartitions() {
        return new Object[][] {
                // number of pairs, number of partitions, expected reads per partition
                { 1, 1, new int[] {2} },
                { 2, 2, new int[] {2, 2} },
                { 3, 2, new int[] {4, 2} },
                { 3, 3, new int[] {2, 2, 2} },
                { 6, 2, new int[] {6, 6} },
                { 6, 3, new int[] {4, 4, 4} },
                { 6, 4, new int[] {4, 2, 4, 2} },
        };
    }

    @Test(dataProvider = "readPairsAndPartitions")
    public void testPutPairsInSamePartition(int numPairs, int numPartitions, int[] expectedReadsPerPartition) throws IOException {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        JavaRDD<GATKRead> reads = createPairedReads(ctx, header, numPairs, numPartitions);
        ReadsSparkSource readsSparkSource = new ReadsSparkSource(ctx);
        JavaRDD<GATKRead> pairedReads = readsSparkSource.putPairsInSamePartition(header, reads);
        List<List<GATKRead>> partitions = pairedReads.mapPartitions((FlatMapFunction<Iterator<GATKRead>, List<GATKRead>>) it ->
                Iterators.singletonIterator(Lists.newArrayList(it))).collect();
        assertEquals(partitions.size(), numPartitions);
        for (int i = 0; i < numPartitions; i++) {
            assertEquals(partitions.get(i).size(), expectedReadsPerPartition[i]);
        }
        assertEquals(Arrays.stream(expectedReadsPerPartition).sum(), numPairs * 2);
    }

    private JavaRDD<GATKRead> createPairedReads(JavaSparkContext ctx, SAMFileHeader header, int numPairs, int numPartitions) {
        final int readSize = 151;
        final int fragmentLen = 400;
        final String templateName = "readpair";
        int leftStart = 10000;
        List<GATKRead> reads = new ArrayList<>();
        for (int i = 0; i < numPairs;i++) {
            leftStart += readSize * 2;
            int rightStart = leftStart + fragmentLen - readSize;
            reads.addAll(ArtificialReadUtils.createPair(header, templateName + i, readSize, leftStart, rightStart, true, false));
        }
        return ctx.parallelize(reads, numPartitions);
    }
}
