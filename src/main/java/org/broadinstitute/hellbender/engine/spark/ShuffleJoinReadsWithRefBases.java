package org.broadinstitute.hellbender.engine.spark;

import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import com.google.common.collect.Lists;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.ReferenceShard;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import scala.Tuple2;

import java.util.List;
import java.util.stream.Collectors;

/**
 * RefBasesForReads queries the Google Genomics API for reference bases overlapping all of the reads.
 *
 * step 1: key reads by reference shards
 *
 * |--------- shard 0 ----------|---------- shard 1 ----------|--------- shard 2 ----------|--------- shard 3 ----------|
 *           |------ read a -----|                             |-- read b --|   |---- read c ----|
 *
 * step 2: group reads by the shard they start in
 *
 *  |--------- shard 0 ----------|
 *            |------ read a -----|
 *
 *
 *  |--------- shard 2 ----------|
 *   |-- read b --|   |---- read c ----|
 *
 *  step 3: query the Google Genomics API for all bases needed for each shard
 *

 * |--- ref bases 1 ---|        |--------- ref bases 2 -----------|
 * |------ read a -----|        |-- read b --|   |---- read c ----|
 *
 *  step 4: pair the ref bases needed for each read with the read
 *
 * |------ read a -----|        |-- read b --|   |---- read c ----|
 * |-- ref bases 1a ---|        |ref bases 2b|   |- ref bases 2c -|
 *
 * or in code,
 *  KV<read a, ref bases 1a>
 *  KV<read b, ref bases 2b>
 *  KV<read c, ref bases 2c>
 *
 * The reference bases paired with each read can be customized by passing in a reference window function
 * inside the {@link ReferenceMultiSparkSource} argument to {@link #addBases}. See {@link ReferenceWindowFunctions} for examples.
 */
public final class ShuffleJoinReadsWithRefBases {

    /**
     * Joins each read of an RDD<GATKRead> with that read's corresponding reference sequence.
     *
     * @param referenceDataflowSource The source of the reference sequence information
     * @param reads The reads for which to extract reference sequence information
     * @return The JavaPairRDD that contains each read along with the corresponding ReferenceBases object
     */
    public static JavaPairRDD<GATKRead, ReferenceBases> addBases(final ReferenceMultiSparkSource referenceDataflowSource,
                                                                 final JavaRDD<GATKRead> reads) {
        // TODO: reimpl this method by calling out to the more complex version?
        SerializableFunction<GATKRead, SimpleInterval> windowFunction = referenceDataflowSource.getReferenceWindowFunction();

        JavaPairRDD<ReferenceShard, GATKRead> shardRead = reads.mapToPair(gatkRead -> {
            ReferenceShard shard = ReferenceShard.getShardNumberFromInterval(windowFunction.apply(gatkRead));
            return new Tuple2<>(shard, gatkRead);
        });

        JavaPairRDD<ReferenceShard, Iterable<GATKRead>> shardiRead = shardRead.groupByKey();

        return shardiRead.flatMapToPair(in -> {
            List<Tuple2<GATKRead, ReferenceBases>> out = Lists.newArrayList();
            Iterable<GATKRead> iReads = in._2();

            // Apply the reference window function to each read to produce a set of intervals representing
            // the desired reference bases for each read.
            final List<SimpleInterval> readWindows = Utils.stream(iReads).map(read -> windowFunction.apply(read)).collect(Collectors.toList());

            SimpleInterval interval = IntervalUtils.getSpanningInterval(readWindows);
            ReferenceBases bases = referenceDataflowSource.getReferenceBases(interval);
            for (GATKRead r : iReads) {
                final ReferenceBases subset = bases.getSubset(windowFunction.apply(r));
                out.add(new Tuple2<>(r, subset));
            }
            return out.iterator();
        });
    }

    /**
     * Joins each read of an RDD<GATKRead, T> with key's corresponding reference sequence.
     *
     * @param referenceDataflowSource The source of the reference sequence information
     * @param keyedByRead The read-keyed RDD for which to extract reference sequence information
     * @return The JavaPairRDD that contains each read along with the corresponding ReferenceBases object and the value
     */
    public static <T> JavaPairRDD<GATKRead, Tuple2<T, ReferenceBases>> addBases(final ReferenceMultiSparkSource referenceDataflowSource,
                                                                                final JavaPairRDD<GATKRead, T> keyedByRead) {
        SerializableFunction<GATKRead, SimpleInterval> windowFunction = referenceDataflowSource.getReferenceWindowFunction();

        JavaPairRDD<ReferenceShard, Tuple2<GATKRead, T>> shardRead = keyedByRead.mapToPair(pair -> {
            ReferenceShard shard = ReferenceShard.getShardNumberFromInterval(windowFunction.apply(pair._1()));
            return new Tuple2<>(shard, pair);
        });

        JavaPairRDD<ReferenceShard, Iterable<Tuple2<GATKRead, T>>> shardiRead = shardRead.groupByKey();

        return shardiRead.flatMapToPair(in -> {
            List<Tuple2<GATKRead, Tuple2<T, ReferenceBases>>> out = Lists.newArrayList();
            Iterable<Tuple2<GATKRead, T>> iReads = in._2();

            // Apply the reference window function to each read to produce a set of intervals representing
            // the desired reference bases for each read.
            final List<SimpleInterval> readWindows = Utils.stream(iReads).map(pair -> windowFunction.apply(pair._1())).collect(Collectors.toList());

            SimpleInterval interval = IntervalUtils.getSpanningInterval(readWindows);
            // TODO: don't we need to support GCS PipelineOptions?
            ReferenceBases bases = referenceDataflowSource.getReferenceBases(interval);
            for (Tuple2<GATKRead, T> p : iReads) {
                final ReferenceBases subset = bases.getSubset(windowFunction.apply(p._1()));
                out.add(new Tuple2<>(p._1(), new Tuple2<>(p._2(), subset)));
            }
            return out.iterator();
        });
    }
}
