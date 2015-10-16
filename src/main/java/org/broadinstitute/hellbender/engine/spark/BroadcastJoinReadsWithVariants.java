package org.broadinstitute.hellbender.engine.spark;

import com.google.common.collect.Lists;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.Variant;
import scala.Tuple2;

import java.util.List;

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
        Broadcast<List<Variant>> variantsBroadcast = ctx.broadcast(variants.collect());
        return reads.mapToPair(r -> {
            List<Variant> overlappingVariants = Lists.newArrayList();
            SimpleInterval interval = new SimpleInterval(r);
            for (Variant v : variantsBroadcast.getValue()) {
                if (interval.overlaps(v)) {
                    overlappingVariants.add(v);
                }
            }
            return new Tuple2<>(r, overlappingVariants);
        });
    }
}