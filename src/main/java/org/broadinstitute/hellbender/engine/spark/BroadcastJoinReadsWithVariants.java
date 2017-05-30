package org.broadinstitute.hellbender.engine.spark;

import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipList;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import scala.Tuple2;

import java.util.Collections;
import java.util.List;

/**
 * Joins an RDD of GATKReads to variant data using a broadcast strategy.
 *
 * The variants RDD is materialized as a List then broadcast using Spark's Broadcast variable mechanism.  The reads are
 * then mapped over and overlapping variants are added for each read.
 */
public final class BroadcastJoinReadsWithVariants {
    private BroadcastJoinReadsWithVariants(){}

    public static JavaPairRDD<GATKRead, Iterable<GATKVariant>> join(final JavaRDD<GATKRead> reads, final JavaRDD<GATKVariant> variants, final List<String> variantsPaths ) {
        final JavaSparkContext ctx = new JavaSparkContext(reads.context());
        final Broadcast<IntervalsSkipList<GATKVariant>> variantsBroadcast = variantsPaths == null ? ctx.broadcast(new IntervalsSkipList<>(variants.collect())) : null;

        return reads.mapToPair(r -> {
            final IntervalsSkipList<GATKVariant> intervalsSkipList = variantsPaths == null ? variantsBroadcast.getValue() :
                    KnownSitesCache.getVariants(variantsPaths);
            if (SimpleInterval.isValid(r.getContig(), r.getStart(), r.getEnd())) {
                return new Tuple2<>(r, intervalsSkipList.getOverlapping(new SimpleInterval(r)));
            } else {
                //Sometimes we have reads that do not form valid intervals (reads that do not consume any ref bases, eg CIGAR 61S90I
                //In those cases, we'll just say that nothing overlaps the read
                return new Tuple2<>(r, Collections.emptyList());
            }
        });
    }
}