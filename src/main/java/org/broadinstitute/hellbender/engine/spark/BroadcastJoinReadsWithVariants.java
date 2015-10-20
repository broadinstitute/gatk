package org.broadinstitute.hellbender.engine.spark;

import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipList;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.Variant;
import scala.Tuple2;

/**
 * Joins an RDD of GATKReads to variant data using a broadcast strategy.
 *
 * The variants RDD is materialized as a List then broadcast using Spark's Broadcast variable mechanism.  The reads are
 * then mapped over and overlapping variants are added for each read.
 */
public class BroadcastJoinReadsWithVariants {
    public static JavaPairRDD<GATKRead, Iterable<Variant>> join(
            final JavaRDD<GATKRead> reads, final JavaRDD<Variant> variants) {
        JavaSparkContext ctx = new JavaSparkContext(reads.context());
        final IntervalsSkipList<Variant> variantSkipList = new IntervalsSkipList<>(variants.collect());
        Broadcast<IntervalsSkipList<Variant>> variantsBroadcast = ctx.broadcast(variantSkipList);

        return reads.mapToPair(r -> {
            IntervalsSkipList<Variant> intervalsSkipList = variantsBroadcast.getValue();
            return new Tuple2<>(r, intervalsSkipList.getOverlapping(new SimpleInterval(r)));
        });
    }
}