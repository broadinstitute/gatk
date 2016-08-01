package org.broadinstitute.hellbender.tools.spark.pipelines;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.collect.Iterators;
import com.google.cloud.dataflow.sdk.repackaged.com.google.common.collect.PeekingIterator;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.Function2;
import org.apache.spark.api.java.function.PairFlatMapFunction;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.*;
import java.util.stream.Collectors;

@CommandLineProgramProperties(summary = "Calculates the read coverage (depth) in the input SAM/BAM",
        oneLineSummary = "Calculate coverage on Spark",
        programGroup = SparkProgramGroup.class)
public final class CalculateCoverageSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc="the output file path", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String outputFile;

    @Argument(doc = "whether to use the shuffle implementation or not", shortName = "shuffle", fullName = "shuffle", optional = true)
    public boolean shuffle = false;

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final JavaRDD<GATKRead> reads = getReads();
        JavaPairRDD<Integer, Integer> coverage = shuffle ?
                calculateCoverage(reads) : calculateCoverageBatched(reads);
        coverage.saveAsTextFile(outputFile);
    }

    private static JavaPairRDD<Integer, Integer> calculateCoverage(final JavaRDD<GATKRead> reads) {
        return reads.flatMapToPair((PairFlatMapFunction<GATKRead, Integer, Integer>) read -> {
            List<Tuple2<Integer, Integer>> depths = new ArrayList<>();
            for (int pos = read.getStart(); pos <= read.getEnd(); pos++) {
                depths.add(new Tuple2<>(pos, 1));
            }
            return depths;
        }).foldByKey(0, (Function2<Integer, Integer, Integer>) (v1, v2) -> v1 + v2);
    }

    private static JavaPairRDD<Integer, Integer> calculateCoverageBatched(final JavaRDD<GATKRead> reads) {
        // Strategy:
        // 1. work out depths for all reads in partition (keep track of partition range while doing this)
        // 2. emit coverage ranges (ranges of N positions, e.g. N=10^5) using SortedMap#submap
        // 3. aggregate by coverage range (to account for ranges that span more than one partition)
        // 4. write ranges out in an appropriate format (e.g. parquet contig/pos/depth)

        JavaPairRDD<Integer, SortedMap<Integer, Integer>> rangesRdd = reads.mapPartitions((FlatMapFunction<Iterator<GATKRead>, Tuple2<Integer, SortedMap<Integer, Integer>>>) it -> {
            int rangeSize = 100000;
            // TODO: loop overs steps 1. and 2. per contig. Use Iterators.peekingIterator
            // 1.
            int minStart = 0;
            int maxEnd = 0;
            SortedMap<Integer, Integer> partitionDepths = new TreeMap<>(); // TODO: use a more efficient sorted map for ints (e.g. fastutil's Int2IntSortedMap)
            while (it.hasNext()) {
                GATKRead read = it.next();
                for (int pos = read.getStart(); pos <= read.getEnd(); pos++) {
                    Integer count = partitionDepths.get(pos);
                    partitionDepths.put(pos, count == null ? 1 : count + 1);
                }
                minStart = Math.min(minStart, read.getStart());
                maxEnd = Math.max(maxEnd, read.getEnd());
            }
            // 2.
            List<Tuple2<Integer, SortedMap<Integer, Integer>>> ranges = new ArrayList<>();
            int firstRangeStart = rangeSize * (minStart / rangeSize);
            int lastRangeStart = rangeSize * (maxEnd / rangeSize);
            for (int rangeStart = firstRangeStart; rangeStart <= lastRangeStart; rangeStart += rangeSize) {
                SortedMap<Integer, Integer> range = partitionDepths.subMap(rangeStart, rangeStart + rangeSize);
                ranges.add(new Tuple2<>(rangeStart, range));
            }
            return ranges;
        }).mapToPair(p -> p);

        // 3.
        JavaPairRDD<Integer, TreeMap<Integer, Integer>> rdd = rangesRdd.aggregateByKey(new TreeMap<>(), (Function2<TreeMap<Integer, Integer>, SortedMap<Integer, Integer>, TreeMap<Integer, Integer>>) (v1, v2) -> {
            v1.putAll(v2);
            return v1;
        }, (Function2<TreeMap<Integer, Integer>, TreeMap<Integer, Integer>, TreeMap<Integer, Integer>>) (v1, v2) -> {
            v1.putAll(v2);
            return v1;
        });

        // 4.
        // TODO: rather than creating Tuple2 objects, it would be better to write out directly from the previous step, perhaps with custom output format
        return rdd.flatMapToPair((PairFlatMapFunction<Tuple2<Integer, TreeMap<Integer, Integer>>, Integer, Integer>) pair -> pair._2.entrySet().stream().map(entry -> new Tuple2<>(entry.getKey(), entry.getValue())).collect(Collectors.toList()));
    }

}
