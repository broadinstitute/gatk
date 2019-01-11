package org.broadinstitute.hellbender.utils.spark;

import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import htsjdk.samtools.*;
import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.ReadQueryNameComparator;
import org.broadinstitute.hellbender.testutils.MiniClusterUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.testng.Assert.assertEquals;

public class SparkUtilsUnitTest extends GATKBaseTest {

    @Test
    public void testConvertHeaderlessHadoopBamShardToBam() {
        final File bamShard = new File(publicTestDir + "org/broadinstitute/hellbender/utils/spark/reads_data_source_test1.bam.headerless.part-r-00000");
        final File output = createTempFile("testConvertHadoopBamShardToBam", ".bam");
        final File headerSource = new File(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam");
        final int expectedReadCount = 11;

        boolean shardIsNotValidBam = false;
        try ( final ReadsDataSource readsSource = new ReadsDataSource(bamShard.toPath()) ) {
            for ( final GATKRead read : readsSource ) {}
        }
        catch ( SAMFormatException e ) {
            shardIsNotValidBam = true;
        }

        Assert.assertTrue(shardIsNotValidBam, "Input shard should not be a valid BAM");

        SAMFileHeader header = null;
        try ( final SamReader headerReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(headerSource) ) {
            header = headerReader.getFileHeader();
        }
        catch ( IOException e ) {
            throw new UserException("Error reading header from " + headerSource.getAbsolutePath(), e);
        }

        SparkUtils.convertHeaderlessHadoopBamShardToBam(bamShard, header, output);

        int actualCount = 0;
        try ( final ReadsDataSource readsSource = new ReadsDataSource(output.toPath()) ) {
            for ( final GATKRead read : readsSource ) { ++actualCount; }
        }

        Assert.assertEquals(actualCount, expectedReadCount, "Wrong number of reads in final BAM file");
    }

    @Test
    public void testPathExists() throws Exception {
        MiniClusterUtils.runOnIsolatedMiniCluster( cluster -> {
            //use the HDFS on the mini cluster
            final Path workingDirectory = MiniClusterUtils.getWorkingDir(cluster);
            final Path tempPath = new Path(workingDirectory, "testFileExists.txt");
            final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

            Assert.assertFalse(SparkUtils.pathExists(ctx, tempPath));
            final FileSystem fs = tempPath.getFileSystem(ctx.hadoopConfiguration());
            final FSDataOutputStream fsOutStream = fs.create(tempPath);
            fsOutStream.close();
            fs.deleteOnExit(tempPath);
            Assert.assertTrue(SparkUtils.pathExists(ctx, tempPath));
        });
    }

    @DataProvider(name="readPairsAndPartitions")
    public Object[][] readPairsAndPartitions() {
        return new Object[][] {
                // number of pairs, number of partitions, number of reads per pair,  expected reads per partition
                { 1, 1, 2, new int[] {2} },
                { 2, 2, 2, new int[] {4, 0} },
                { 3, 2, 2, new int[] {4, 2} },
                { 3, 3, 2, new int[] {4, 2, 0} },
                { 6, 2, 2, new int[] {8, 4} },
                { 6, 3, 2, new int[] {6, 4, 2} },
                { 6, 4, 2, new int[] {4, 4, 2, 2} },
                { 2, 2, 3, new int[] {6, 0} },
                { 3, 2, 10, new int[] {20, 10} },
                { 6, 4, 3, new int[] {6, 6, 3, 3} },
                { 20, 7, 5, new int[] {15, 15, 15, 15, 15, 15, 10} },
        };
    }

    @Test(dataProvider = "readPairsAndPartitions")
    public void testPutReadsWithSameNameInSamePartition(int numPairs, int numPartitions, int numReadsInPair, int[] expectedReadsPerPartition) {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        JavaRDD<GATKRead> reads =  ctx.parallelize(createPairedReads(header, numPairs, numReadsInPair), numPartitions);
        JavaRDD<GATKRead> pairedReads = SparkUtils.putReadsWithTheSameNameInTheSamePartition(header, reads, ctx);
        List<List<GATKRead>> partitions = pairedReads.mapPartitions((FlatMapFunction<Iterator<GATKRead>, List<GATKRead>>) it ->
                Iterators.singletonIterator(Lists.newArrayList(it))).collect();
        assertEquals(partitions.size(), numPartitions);
        for (int i = 0; i < numPartitions; i++) {
            assertEquals(partitions.get(i).size(), expectedReadsPerPartition[i]);
        }
        assertEquals(Arrays.stream(expectedReadsPerPartition).sum(), numPairs * numReadsInPair);
    }

    private static List<GATKRead> createPairedReads(SAMFileHeader header, int numPairs, int numReadsInPair) {
        final int readSize = 151;
        final int fragmentLen = 400;
        final String templateName = "readpair";
        int leftStart = 10000;
        List<GATKRead> reads = new ArrayList<>();
        for (int i = 0; i < numPairs;i++) {
            leftStart += readSize * 2;
            int rightStart = leftStart + fragmentLen - readSize;
            reads.addAll(ArtificialReadUtils.createPair(header, templateName + i, readSize, leftStart, rightStart, true, false));
            // Copying a secondary alignment for the second read to fill out the read group
            GATKRead readToCopy = reads.get(reads.size()-1).copy();
            readToCopy.setIsSecondaryAlignment(true);
            for (int j = 2; j < numReadsInPair; j++) {
                reads.add(readToCopy.copy());
            }
        }
        return reads;
    }

    @Test(expectedExceptions = GATKException.class)
    public void testReadsPairsSpanningMultiplePartitionsCrash() {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        List<GATKRead> reads = createPairedReads(header, 40, 2);
        // Creating one group in the middle that should cause problems
        reads.addAll(40, createPairedReads(header, 1, 30));

        JavaRDD<GATKRead> problemReads = ctx.parallelize(reads,5 );
        SparkUtils.putReadsWithTheSameNameInTheSamePartition(header, problemReads, ctx);
    }

    @Test
    public void testReadsMustBeQueryGroupedToFixPartitions(){
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        List<GATKRead> reads = createPairedReads(header, 40, 2);
        final JavaRDD<GATKRead> readsRDD = ctx.parallelize(reads, 5);
        Assert.assertThrows(IllegalArgumentException.class, () -> SparkUtils.putReadsWithTheSameNameInTheSamePartition(header, readsRDD, ctx));
    }

    @Test(expectedExceptions = GATKException.class)
    public void testInvalidSort(){
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        List<GATKRead> reads = new ArrayList<>();
        for(int i = 0; i < 10; i++){
            //create reads with alternating contigs and  decreasing start position
            reads.add(ArtificialReadUtils.createArtificialRead(header, "READ"+i, i % header.getSequenceDictionary().size() , 100, 100));
        }
        final JavaRDD<GATKRead> readsRDD = ctx.parallelize(reads);
        SparkUtils.sortReadsAccordingToHeader(readsRDD, header, 0);
    }

    @Test
    public void testSortCoordinateSortMatchesHtsjdk() {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        List<GATKRead> reads = new ArrayList<>();
        for(int i = 0; i < 2000; i++){
            //create reads with alternating contigs and  decreasing start position
            reads.add(ArtificialReadUtils.createArtificialRead(header, "READ"+i, i % header.getSequenceDictionary().size() , 3000 - i, 100));
        }
        final JavaRDD<GATKRead> readsRDD = ctx.parallelize(reads);
        final List<GATKRead> coordinateSorted = SparkUtils.sortReadsAccordingToHeader(readsRDD, header, 0).collect();
        assertSorted(coordinateSorted, new ReadCoordinateComparator(header));
        assertSorted(coordinateSorted.stream().map(read -> read.convertToSAMRecord(header)).collect(Collectors.toList()), new SAMRecordCoordinateComparator());
    }

    @Test
    public void testSortQuerynameSortMatchesHtsjdk() {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        List<GATKRead> reads = new ArrayList<>();
        final int numReads = 2000;
        for(int i = 0; i < numReads; i++) {

            //create reads with non-lexicographically ordered names
            //names are created in lexicographically decreasing order, with 2 repetitions to create "pairs"
            reads.add(ArtificialReadUtils.createArtificialRead(header, "READ" + (numReads - i) % (numReads / 2),
                                                               i % header.getSequenceDictionary().size(),
                                                               3000 - i,
                                                               100));
        }
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        final JavaRDD<GATKRead> readsRDD = ctx.parallelize(reads);
        final List<GATKRead> querynameSorted = SparkUtils.sortReadsAccordingToHeader(readsRDD, header, 31).collect();
        assertSorted(querynameSorted, new ReadQueryNameComparator());
        assertSorted(querynameSorted.stream().map(read -> read.convertToSAMRecord(header)).collect(Collectors.toList()), new SAMRecordQueryNameComparator());
    }

    @Test
    public void testSortQuerynameFixesPartitionBoundaries(){
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        final int numReadsWithSameName = 4;
        final List<GATKRead> pairedReads = createPairedReads(header, 100, numReadsWithSameName);
        final int numPartitions = 7;
        final JavaRDD<GATKRead> reads = ctx.parallelize(pairedReads, numPartitions);

        //assert that the grouping is not correct before sorting
        final List<GATKRead>[] partitions = reads.collectPartitions(IntStream.range(0, reads.getNumPartitions()).toArray());
        Assert.assertTrue(
                Arrays.stream(partitions)
                        //look through each partition and count the number of each read name seen
                        .flatMap( readsInPartition -> readsInPartition.stream()
                            .collect(Collectors.groupingBy(GATKRead::getName))
                            .values()
                            .stream()
                            .map(List::size)
                        )
                        //check that at least one partition was not correctly distributed
                .anyMatch(size -> size != numReadsWithSameName), "The partitioning was correct before sorting so the test is meaningless.");

        final JavaRDD<GATKRead> sorted = SparkUtils.sortReadsAccordingToHeader(reads, header, numPartitions);

        //assert that the grouping is fixed after sorting
        final List<GATKRead>[] sortedPartitions = sorted.collectPartitions(IntStream.range(0, sorted.getNumPartitions()).toArray());
        Assert.assertTrue(Arrays.stream(sortedPartitions)
                .flatMap( readsInPartition -> readsInPartition.stream()
                        .collect(Collectors.groupingBy(GATKRead::getName))
                        .values()
                        .stream()
                        .map(List::size)
                )
                .allMatch(size -> size == numReadsWithSameName), "Some reads names were split between multiple partitions");
    }

    @Test
    public void testSortUsingElementsAsKeys(){
        final List<Integer> unsorted = Arrays.asList(4, 2, 6, 0, 8);
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final JavaRDD<Integer> unsortedRDD = ctx.parallelize(unsorted);
        final JavaRDD<Integer> sorted = SparkUtils.sortUsingElementsAsKeys(unsortedRDD, Comparator.naturalOrder(), 2);
        assertSorted(sorted.collect(), Comparator.naturalOrder());
    }
}
