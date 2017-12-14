package org.broadinstitute.hellbender.engine.spark;

import com.google.common.collect.*;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.PairFlatMapFunction;
import org.broadinstitute.hellbender.engine.Shard;
import org.broadinstitute.hellbender.engine.ShardBoundary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.IOException;
import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;

import static org.testng.Assert.*;

public class SparkSharderUnitTest extends GATKBaseTest implements Serializable {

    private static final long serialVersionUID = 1L;

    private static final int STANDARD_READ_LENGTH = 3;

    private SAMSequenceDictionary sequenceDictionary = new SAMSequenceDictionary(
            ImmutableList.of(new SAMSequenceRecord("1", 100), new SAMSequenceRecord("2", 50)));


    @Test
    public void testLocatablesPerShard() throws IOException {

        // Consider the following reads (divided into four partitions), and intervals.
        // This test checks that iterating over the reads and intervals at the same time produces the overlapping
        // reads for each interval.

        //                      1                   2
        //    1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7
        // ---------------------------------------------------------
        // Reads in partition 0
        //   [-----]
        //           [-----]
        //               [-----]
        // ---------------------------------------------------------
        // Reads in partition 1
        //               [-----]
        //               [-----]
        //               [-----]
        // ---------------------------------------------------------
        // Reads in partition 2
        //               [-----]
        //                       [-----]
        //                         [-----]
        // ---------------------------------------------------------
        // Reads in partition 3
        //                                   [-----]
        //                                           [-----]
        //                                                   [-----]
        // ---------------------------------------------------------
        // Per-partition read extents
        //   [-----------------]
        //               [-----]
        //               [---------------]
        //                                   [---------------------]
        // ---------------------------------------------------------
        // Intervals
        //     [-----]
        //                 [---------]
        //                       [-----------------------]
        //
        //                      1                   2
        // ---------------------------------------------------------
        //    1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7

        List<TestRead> reads = ImmutableList.of(
                new TestRead(1, 3), new TestRead(5, 7), new TestRead(7, 9),
                new TestRead(7, 9), new TestRead(7, 9), new TestRead(7, 9),
                new TestRead(7, 9), new TestRead(11, 13), new TestRead(12, 14),
                new TestRead(17, 19), new TestRead(21, 23), new TestRead(25, 27)
        );

        List<SimpleInterval> intervals = ImmutableList.of(
                new SimpleInterval("1", 2, 4),
                new SimpleInterval("1", 8, 12),
                new SimpleInterval("1", 11, 22));

        Iterator<Tuple2<SimpleInterval, Iterable<TestRead>>> it = SparkSharder.locatablesPerShard(reads.iterator(), intervals.iterator(), sequenceDictionary, STANDARD_READ_LENGTH);
        assertTrue(it.hasNext());
        Tuple2<SimpleInterval, Iterable<TestRead>> next = it.next();
        assertEquals(next._1(), intervals.get(0));
        assertEquals(next._2(), ImmutableList.of(reads.get(0)));

        assertTrue(it.hasNext());
        next = it.next();
        assertEquals(next._1(), intervals.get(1));
        assertEquals(next._2(), ImmutableList.of(reads.get(2), reads.get(3), reads.get(4), reads.get(5), reads.get(6), reads.get(7), reads.get(8)));

        assertTrue(it.hasNext());
        next = it.next();
        assertEquals(next._1(), intervals.get(2));
        assertEquals(next._2(), ImmutableList.of(reads.get(7), reads.get(8), reads.get(9), reads.get(10)));

        assertFalse(it.hasNext());
    }

    @Test
    public void testSingleContig() throws IOException {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        // Consider the following reads (divided into four partitions), and intervals.
        // This test counts the number of reads that overlap each interval.

        //                      1                   2
        //    1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7
        // ---------------------------------------------------------
        // Reads in partition 0
        //   [-----]
        //           [-----]
        //               [-----]
        // ---------------------------------------------------------
        // Reads in partition 1
        //               [-----]
        //               [-----]
        //               [-----]
        // ---------------------------------------------------------
        // Reads in partition 2
        //               [-----]
        //                       [-----]
        //                         [-----]
        // ---------------------------------------------------------
        // Reads in partition 3
        //                                   [-----]
        //                                           [-----]
        //                                                   [-----]
        // ---------------------------------------------------------
        // Per-partition read extents
        //   [-----------------]
        //               [-----]
        //               [---------------]
        //                                   [---------------------]
        // ---------------------------------------------------------
        // Intervals
        //     [-----]
        //                 [---------]
        //                       [-----------------------]
        //
        //                      1                   2
        // ---------------------------------------------------------
        //    1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7

        JavaRDD<TestRead> reads = ctx.parallelize(ImmutableList.of(
                new TestRead(1, 3), new TestRead(5, 7), new TestRead(7, 9),
                new TestRead(7, 9), new TestRead(7, 9), new TestRead(7, 9),
                new TestRead(7, 9), new TestRead(11, 13), new TestRead(12, 14),
                new TestRead(17, 19), new TestRead(21, 23), new TestRead(25, 27)
        ), 4);

        List<SimpleInterval> intervals = ImmutableList.of(
                new SimpleInterval("1", 2, 4),
                new SimpleInterval("1", 8, 12),
                new SimpleInterval("1", 11, 22));

        List<ShardBoundary> shardBoundaries = intervals.stream().map(si -> new ShardBoundary(si, si)).collect(Collectors.toList());

        ImmutableMap<SimpleInterval, Integer> expectedReadsPerInterval = ImmutableMap.of(intervals.get(0), 1, intervals.get(1), 7, intervals.get(2), 4);

        JavaPairRDD<Locatable, Integer> readsPerInterval =
                SparkSharder.shard(ctx, reads, TestRead.class, sequenceDictionary, shardBoundaries, STANDARD_READ_LENGTH, false)
                .flatMapToPair(new CountOverlappingReadsFunction());
        assertEquals(readsPerInterval.collectAsMap(), expectedReadsPerInterval);

        JavaPairRDD<Locatable, Integer> readsPerIntervalShuffle =
                SparkSharder.shard(ctx, reads, TestRead.class, sequenceDictionary, shardBoundaries, STANDARD_READ_LENGTH, true)
                .flatMapToPair(new CountOverlappingReadsFunction());
        assertEquals(readsPerIntervalShuffle.collectAsMap(), expectedReadsPerInterval);

        try {
            int maxReadLength = STANDARD_READ_LENGTH - 1; // max read length less than actual causes exception
            SparkSharder.shard(ctx, reads, TestRead.class, sequenceDictionary, shardBoundaries, maxReadLength, true)
                    .flatMapToPair(new CountOverlappingReadsFunction()).collect();
        } catch (Exception e) {
            assertEquals(e.getCause().getClass(), UserException.class);
        }
    }

    @Test
    public void testContigBoundary() throws IOException {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        // Consider the following reads (divided into four partitions), and intervals.
        // This test counts the number of reads that overlap each interval.

        //                      1                   2
        //    1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7
        // ---------------------------------------------------------
        // Reads in partition 0
        //   [-----] chr 1
        //           [-----] chr 1
        //               [-----] chr 1
        //   [-----] chr 2
        //     [-----] chr 2
        // ---------------------------------------------------------
        // Per-partition read extents
        //   [-----------------] chr 1
        //   [-------] chr 2
        // ---------------------------------------------------------
        // Intervals
        //     [-----] chr 1
        //                 [---------] chr 1
        //   [-----------------------] chr 2
        // ---------------------------------------------------------
        //    1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7

        JavaRDD<TestRead> reads = ctx.parallelize(ImmutableList.of(
                new TestRead("1", 1, 3), new TestRead("1", 5, 7), new TestRead("1", 7, 9),
                new TestRead("2", 1, 3), new TestRead("2", 2, 4)
        ), 1);

        List<SimpleInterval> intervals = ImmutableList.of(
                new SimpleInterval("1", 2, 4),
                new SimpleInterval("1", 8, 12),
                new SimpleInterval("2", 1, 12));

        List<ShardBoundary> shardBoundaries = intervals.stream().map(si -> new ShardBoundary(si, si)).collect(Collectors.toList());

        ImmutableMap<SimpleInterval, Integer> expectedReadsPerInterval = ImmutableMap.of(intervals.get(0), 1, intervals.get(1), 1, intervals.get(2), 2);

        JavaPairRDD<Locatable, Integer> readsPerInterval =
                SparkSharder.shard(ctx, reads, TestRead.class, sequenceDictionary, shardBoundaries, STANDARD_READ_LENGTH, false)
                        .flatMapToPair(new CountOverlappingReadsFunction());
        assertEquals(readsPerInterval.collectAsMap(), expectedReadsPerInterval);

        JavaPairRDD<Locatable, Integer> readsPerIntervalShuffle =
                SparkSharder.shard(ctx, reads, TestRead.class, sequenceDictionary, shardBoundaries, STANDARD_READ_LENGTH, true)
                        .flatMapToPair(new CountOverlappingReadsFunction());
        assertEquals(readsPerIntervalShuffle.collectAsMap(), expectedReadsPerInterval);

    }

    @Test
    public void testPartitionReadExtents() throws IOException {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        // Consider the following reads.
        // This test checks the partition read extents when the reads are divided into
        // different numbers of partitions (1, 2, or 3), and with different sequence dictionaries.

        //                      1                   2
        //    1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7
        // ---------------------------------------------------------
        // Reads
        //   [-----] chr 1
        //           [-----] chr 1
        //               [-----] chr 1
        //                       [-----] chr 2
        //                         [-----] chr 2
        //                           [-----] chr 2
        // ---------------------------------------------------------
        //    1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7

        ImmutableList<TestRead> reads = ImmutableList.of(
                new TestRead("1", 1, 3), new TestRead("1", 5, 7), new TestRead("1", 7, 9),
                new TestRead("2", 11, 13), new TestRead("2", 12, 14), new TestRead("2", 13, 15)
        );

        assertEquals(SparkSharder.computePartitionReadExtents(ctx.parallelize(reads, 1), sequenceDictionary, STANDARD_READ_LENGTH),
                ImmutableList.of(
                        new SparkSharder.PartitionLocatable<>(0, new SimpleInterval("1", 1, 100)),
                        new SparkSharder.PartitionLocatable<>(0, new SimpleInterval("2", 1, 50))
                ));

        assertEquals(SparkSharder.computePartitionReadExtents(ctx.parallelize(reads, 2), sequenceDictionary, STANDARD_READ_LENGTH),
                ImmutableList.of(
                        new SparkSharder.PartitionLocatable<>(0, new SimpleInterval("1", 1, 100)),
                        new SparkSharder.PartitionLocatable<>(0, new SimpleInterval("2", 1, 14)), // since last read of partition 0 _could_ end at start of first read in partition 1 + max read length
                        new SparkSharder.PartitionLocatable<>(1, new SimpleInterval("2", 11, 50))
                ));

        assertEquals(SparkSharder.computePartitionReadExtents(ctx.parallelize(reads, 3), sequenceDictionary, STANDARD_READ_LENGTH),
                ImmutableList.of(
                        new SparkSharder.PartitionLocatable<>(0, new SimpleInterval("1", 1, 10)),
                        new SparkSharder.PartitionLocatable<>(1, new SimpleInterval("1", 7, 100)),
                        new SparkSharder.PartitionLocatable<>(1, new SimpleInterval("2", 1, 15)),
                        new SparkSharder.PartitionLocatable<>(2, new SimpleInterval("2", 12, 50))
                ));

        // Use a different dictionary with contig 1 missing
        SAMSequenceDictionary sequenceDictionaryMissing1 = new SAMSequenceDictionary(
                ImmutableList.of(new SAMSequenceRecord("2", 50)));
        try {
            SparkSharder.computePartitionReadExtents(ctx.parallelize(reads, 1), sequenceDictionaryMissing1, STANDARD_READ_LENGTH);
            fail("Should throw IllegalStateException");
        } catch (IllegalStateException e) {
            // expected
        }

        // Use a different dictionary with contig 2 missing
        SAMSequenceDictionary sequenceDictionaryMissing2 = new SAMSequenceDictionary(
                ImmutableList.of(new SAMSequenceRecord("1", 50)));
        try {
            SparkSharder.computePartitionReadExtents(ctx.parallelize(reads, 2), sequenceDictionaryMissing2, STANDARD_READ_LENGTH);
            fail("Should throw IllegalStateException");
        } catch (IllegalStateException e) {
            // expected
        }

        // Use a different dictionary with contig 3 at the end
        SAMSequenceDictionary sequenceDictionary123 = new SAMSequenceDictionary(
                ImmutableList.of(new SAMSequenceRecord("1", 100), new SAMSequenceRecord("2", 50), new SAMSequenceRecord("3", 25)));

        assertEquals(SparkSharder.computePartitionReadExtents(ctx.parallelize(reads, 1), sequenceDictionary123, STANDARD_READ_LENGTH),
                ImmutableList.of(
                        new SparkSharder.PartitionLocatable<>(0, new SimpleInterval("1", 1, 100)),
                        new SparkSharder.PartitionLocatable<>(0, new SimpleInterval("2", 1, 50)),
                        new SparkSharder.PartitionLocatable<>(0, new SimpleInterval("3", 1, 25)) // partition could contain contig 3 reads
                ));

        assertEquals(SparkSharder.computePartitionReadExtents(ctx.parallelize(reads, 2), sequenceDictionary123, STANDARD_READ_LENGTH),
                ImmutableList.of(
                        new SparkSharder.PartitionLocatable<>(0, new SimpleInterval("1", 1, 100)),
                        new SparkSharder.PartitionLocatable<>(0, new SimpleInterval("2", 1, 14)), // since last read of partition 0 _could_ end at start of first read in partition 1 + max read length
                        new SparkSharder.PartitionLocatable<>(1, new SimpleInterval("2", 11, 50)),
                        new SparkSharder.PartitionLocatable<>(1, new SimpleInterval("3", 1, 25)) // partition could contain contig 3 reads
                ));

        assertEquals(SparkSharder.computePartitionReadExtents(ctx.parallelize(reads, 3), sequenceDictionary123, STANDARD_READ_LENGTH),
                ImmutableList.of(
                        new SparkSharder.PartitionLocatable<>(0, new SimpleInterval("1", 1, 10)),
                        new SparkSharder.PartitionLocatable<>(1, new SimpleInterval("1", 7, 100)),
                        new SparkSharder.PartitionLocatable<>(1, new SimpleInterval("2", 1, 15)),
                        new SparkSharder.PartitionLocatable<>(2, new SimpleInterval("2", 12, 50)),
                        new SparkSharder.PartitionLocatable<>(2, new SimpleInterval("3", 1, 25)) // partition could contain contig 3 reads
                ));

        // Use a different dictionary with contig X between contigs 1 and 2
        SAMSequenceDictionary sequenceDictionary1X2 = new SAMSequenceDictionary(
                ImmutableList.of(new SAMSequenceRecord("1", 100), new SAMSequenceRecord("X", 75), new SAMSequenceRecord("2", 50)));

        assertEquals(SparkSharder.computePartitionReadExtents(ctx.parallelize(reads, 1), sequenceDictionary1X2, STANDARD_READ_LENGTH),
                ImmutableList.of(
                        new SparkSharder.PartitionLocatable<>(0, new SimpleInterval("1", 1, 100)),
                        new SparkSharder.PartitionLocatable<>(0, new SimpleInterval("X", 1, 75)), // partition could contain contig X reads
                        new SparkSharder.PartitionLocatable<>(0, new SimpleInterval("2", 1, 50))
                ));

        assertEquals(SparkSharder.computePartitionReadExtents(ctx.parallelize(reads, 2), sequenceDictionary1X2, STANDARD_READ_LENGTH),
                ImmutableList.of(
                        new SparkSharder.PartitionLocatable<>(0, new SimpleInterval("1", 1, 100)),
                        new SparkSharder.PartitionLocatable<>(0, new SimpleInterval("X", 1, 75)), // partition could contain contig X reads
                        new SparkSharder.PartitionLocatable<>(0, new SimpleInterval("2", 1, 14)), // since last read of partition 0 _could_ end at start of first read in partition 1 + max read length
                        new SparkSharder.PartitionLocatable<>(1, new SimpleInterval("2", 11, 50))
                ));

        assertEquals(SparkSharder.computePartitionReadExtents(ctx.parallelize(reads, 3), sequenceDictionary1X2, STANDARD_READ_LENGTH),
                ImmutableList.of(
                        new SparkSharder.PartitionLocatable<>(0, new SimpleInterval("1", 1, 10)),
                        new SparkSharder.PartitionLocatable<>(1, new SimpleInterval("1", 7, 100)),
                        new SparkSharder.PartitionLocatable<>(1, new SimpleInterval("X", 1, 75)), // partition could contain contig X reads
                        new SparkSharder.PartitionLocatable<>(1, new SimpleInterval("2", 1, 15)),
                        new SparkSharder.PartitionLocatable<>(2, new SimpleInterval("2", 12, 50))
                ));
    }

    private static class TestRead implements Locatable {
        private static final long serialVersionUID = 1L;
        private final String contig;
        private final int start;
        private final int end;

        public TestRead(int start, int end) {
            this("1", start, end);
        }

        public TestRead(String contig, int start, int end) {
            this.contig = contig;
            this.start = start;
            this.end = end;
        }

        @Override
        public String getContig() {
            return contig;
        }

        @Override
        public int getStart() {
            return start;
        }

        @Override
        public int getEnd() {
            return end;
        }

        @Override
        public String toString() {
            return "TestRead{" +
                    "contig='" + contig + '\'' +
                    ", start=" + start +
                    ", end=" + end +
                    '}';
        }
    }

    private static class CountOverlappingReadsFunction implements PairFlatMapFunction<Shard<TestRead>, Locatable, Integer> {
        private static final long serialVersionUID = 1L;

        @Override
        public Iterator<Tuple2<Locatable, Integer>> call(Shard<TestRead> s) throws Exception {
            Locatable interval = s.getInterval();
            Iterator<TestRead> iterator = s.iterator();
            int count = 0;
            OverlapDetector<Locatable> overlapDetector = OverlapDetector.create(ImmutableList.of(interval));
            while (iterator.hasNext()) {
                count += overlapDetector.getOverlaps(iterator.next()).size();
            }
            return ImmutableList.of(new Tuple2<>(interval, count)).iterator();
        }
    }
}
