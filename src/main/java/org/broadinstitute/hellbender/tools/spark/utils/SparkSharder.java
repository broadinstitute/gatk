package org.broadinstitute.hellbender.tools.spark.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.Shard;
import org.broadinstitute.hellbender.engine.ShardBoundary;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipList;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by valentin on 8/24/17.
 */
public final class SparkSharder {

    private final JavaSparkContext ctx;
    private final int shardSize;
    private final int padding;
    private final SAMSequenceDictionary dictionary;
    private final IntervalsSkipList<SimpleInterval> intervals;
    private IntervalsSkipList<ShardBoundary> shards;
    private Broadcast<IntervalsSkipList<ShardBoundary>> shardsBroadcast;

    public SparkSharder(final JavaSparkContext ctx, final SAMSequenceDictionary referenceDictionary, final Collection<SimpleInterval> intervals, final int shardSize, final int shardPadding) {
        Utils.nonNull(ctx, "the input context must not be null");
        Utils.nonNull(referenceDictionary, "the input reference dictionary must not be null");
        ParamUtils.isPositive(shardSize, "the input shard-size must be positive");
        this.intervals = new IntervalsSkipList<>(intervals == null ? IntervalUtils.getAllIntervalsForReference(referenceDictionary) : IntervalUtils.sortAndMergeIntervals(intervals, referenceDictionary, IntervalMergingRule.ALL));
        this.dictionary = new SAMSequenceDictionary(Utils.nonNull(referenceDictionary).getSequences());
        this.ctx = Utils.nonNull(ctx, "the input context cannot be null");
        this.shardSize = shardSize;
        this.padding = shardPadding;
    }

    public IntervalsSkipList<ShardBoundary> shards() {
        if (shards == null) {
            shards = new IntervalsSkipList<>(Utils.stream(intervals)
                    .flatMap(s -> Shard.divideIntervalIntoShards(s, shardSize, padding, dictionary).stream())
                    .collect(Collectors.toList()));
        }
        return shards;
    }

    private Broadcast<IntervalsSkipList<ShardBoundary>> getShardsBroadcast() {
        if (shardsBroadcast == null) {
            shardsBroadcast = ctx.broadcast(shards());
        }
        return shardsBroadcast;
    }

    public <T> ShardedRDD<T> shard(final JavaRDD<T> rdd, final Function<? super T, ? extends Iterable<SimpleInterval>> toLoc) {
        return new ShardedRDD<>(rdd.flatMapToPair(t -> {
            final IntervalsSkipList<ShardBoundary> shards = getShardsBroadcast().getValue();
            return Utils.stream(toLoc.call(t))
                    .flatMap(l -> shards.getOverlapping(l).stream())
                    .map(x -> new Tuple2<>(x, t)).iterator();
        }), toLoc);
    }

    public <T, U> ShardedPairRDD<T, U> join(final ShardedRDD<T> t, final ShardedRDD<U> u) {
        final JavaPairRDD<ShardBoundary, List<T>> reducedTs = reduceByKey(t);
        final JavaPairRDD<ShardBoundary, List<U>> reducedUs = reduceByKey(u);
        return new ShardedPairRDD<T, U>(reducedTs.join(reducedUs), t.toLoc, u.toLoc);
    }

    private <T> JavaPairRDD<ShardBoundary, List<T>> reduceByKey(ShardedRDD<T> t) {
        return t.rdd.mapToPair(tt -> {
            final List<T> v = new ArrayList<>();
            v.add(tt._2());
            return new Tuple2<>(tt._1(), v);
        }).reduceByKey((a, b) -> { a.addAll(b); return a; });
    }
}
