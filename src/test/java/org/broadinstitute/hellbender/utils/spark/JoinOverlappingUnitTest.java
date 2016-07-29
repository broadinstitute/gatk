package org.broadinstitute.hellbender.utils.spark;

import com.google.common.collect.*;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.api.java.function.FlatMapFunction2;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.IOException;
import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;

import static org.testng.Assert.assertEquals;

public class JoinOverlappingUnitTest extends BaseTest implements Serializable {

    private static final long serialVersionUID = 1L;

    @Test
    public void testSingleContig() throws IOException {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        //                      1                   2
        //    1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7
        //   [-----]
        //           [-----]
        //               [-----]
        // ---------------------------------------------------------
        //               [-----]
        //               [-----]
        //               [-----]
        // ---------------------------------------------------------
        //               [-----]
        //                       [-----]
        //                         [-----]
        // ---------------------------------------------------------
        //                                   [-----]
        //                                           [-----]
        //                                                   [-----]
        // ---------------------------------------------------------
        //   [-----------------]                                     <- per-partition read extents
        //               [-----]
        //               [---------------]
        //                                   [---------------------]
        // ---------------------------------------------------------
        //
        //     [-----]                                                 <- intervals
        //                 [---------]
        //                       [-----------------------]
        //
        //                      1                   2
        //    1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7

        // maxEndPartitionIndexes: [2, 1, 3, 3]

        JavaRDD<TestRead> reads = ctx.parallelize(ImmutableList.of(
                new TestRead(1, 3), new TestRead(5, 7), new TestRead(7, 9),
                new TestRead(7, 9), new TestRead(7, 9), new TestRead(7, 9),
                new TestRead(7, 9), new TestRead(11, 13), new TestRead(12, 14),
                new TestRead(17, 19), new TestRead(21, 23), new TestRead(25, 27)
        ), 4);

        List<Locatable> intervals = ImmutableList.of(
                new SimpleInterval("1", 2, 4),
                new SimpleInterval("1", 8, 12),
                new SimpleInterval("1", 11, 22));

        JavaPairRDD<Locatable, Integer> readsPerInterval = SparkUtils.joinOverlapping(ctx, reads, TestRead.class, intervals,
                new CountOverlappingReadsFunction()).mapToPair(t -> t);
        assertEquals(readsPerInterval.collectAsMap(), ImmutableMap.of(intervals.get(0), 1, intervals.get(1), 7, intervals.get(2), 4));

        JavaPairRDD<Locatable, Integer> readsPerIntervalNaive = SparkUtils.joinOverlappingShuffle(ctx, reads, TestRead.class, intervals,
                new CountOverlappingReadsFlatFunction()).mapToPair(t -> t);
        assertEquals(readsPerIntervalNaive.collectAsMap(), ImmutableMap.of(intervals.get(0), 1, intervals.get(1), 7, intervals.get(2), 4));
    }

    @Test
    public void testContigBoundary() throws IOException {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        //                      1                   2
        //    1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7
        //   [-----] chr 1
        //           [-----] chr 1
        //               [-----] chr 1
        //                       [-----] chr 2
        //                         [-----] chr 2
        // ---------------------------------------------------------
        //   [-----------------] [-------]                             <- per-partition read extents
        // ---------------------------------------------------------
        //
        //     [-----] chr 1                                           <- intervals
        //                 [---------] chr 1
        //                       [-----------------------] chr 2

        JavaRDD<TestRead> reads = ctx.parallelize(ImmutableList.of(
                new TestRead("1", 1, 3), new TestRead("1", 5, 7), new TestRead("1", 7, 9),
                new TestRead("2", 11, 13), new TestRead("2", 12, 14)
        ), 1);

        List<Locatable> intervals = ImmutableList.of(
                new SimpleInterval("1", 2, 4),
                new SimpleInterval("1", 8, 12),
                new SimpleInterval("2", 11, 22));

        JavaPairRDD<Locatable, Integer> readsPerInterval = SparkUtils.joinOverlapping(ctx, reads, TestRead.class, intervals,
                new CountOverlappingReadsFunction()).mapToPair(t -> t);
        assertEquals(readsPerInterval.collectAsMap(), ImmutableMap.of(intervals.get(0), 1, intervals.get(1), 1, intervals.get(2), 2));

        JavaPairRDD<Locatable, Integer> readsPerIntervalNaive = SparkUtils.joinOverlappingShuffle(ctx, reads, TestRead.class, intervals,
                new CountOverlappingReadsFlatFunction()).mapToPair(t -> t);
        assertEquals(readsPerIntervalNaive.collectAsMap(), ImmutableMap.of(intervals.get(0), 1, intervals.get(1), 1, intervals.get(2), 2));

    }

    static class TestRead implements Locatable {
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
    }

    @SuppressWarnings("unchecked")
    private static class CountOverlappingReadsFunction implements FlatMapFunction2<Iterator<TestRead>, Iterator<Locatable>, Tuple2<Locatable, Integer>> {
        private static final long serialVersionUID = 1L;

        @Override
        public Iterable<Tuple2<Locatable, Integer>> call(Iterator<TestRead> readIterator,
                                           Iterator<Locatable> intervalIterator) {
            List<Locatable> intervals = Lists.newArrayList(intervalIterator);
            List<OverlapDetector<Locatable>> overlapDetectors = intervals.stream()
                    .map(interval -> OverlapDetector.create(ImmutableList.of(interval))).collect(Collectors.toList());

            Map<Locatable, Integer> counts = Maps.newLinkedHashMap();
            while (readIterator.hasNext()) {
                TestRead read = readIterator.next();
                for (OverlapDetector<Locatable> overlapDetector : overlapDetectors) {
                    Set<Locatable> overlaps = overlapDetector.getOverlaps(read);
                    for (Locatable overlap : overlaps) {
                        Integer count = counts.get(overlap);
                        counts.put(overlap, count == null ? 1 : count + 1);
                    }
                }
            }
            return counts.entrySet().stream()
                    .map(entry -> new Tuple2<>(entry.getKey(), entry.getValue()))
                    .collect(Collectors.toList());
        }
    }

    private static class CountOverlappingReadsFlatFunction implements FlatMapFunction<Tuple2<Locatable, Iterable<TestRead>>, Tuple2<Locatable, Integer>> {
        private static final long serialVersionUID = 1L;

        @Override
        public Iterable<Tuple2<Locatable, Integer>> call(Tuple2<Locatable, Iterable<TestRead>> t) throws Exception {
            Locatable interval = t._1;
            int count = 0;
            OverlapDetector<Locatable> overlapDetector = OverlapDetector.create(ImmutableList.of(interval));
            for (TestRead read : t._2) {
                count += overlapDetector.getOverlaps(read).size();
            }
            return ImmutableList.of(new Tuple2<>(interval, count));
        }
    }
}
