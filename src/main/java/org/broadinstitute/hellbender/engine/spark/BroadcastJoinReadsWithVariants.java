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

    /**
     * Joins each read of an RDD<GATKRead> with overlapping variants from an RDD of GATKVariants. Broadcasts the
     * variants, so this is only suitable for collections of variants that are <2GB (due to Spark's broadcast limitation).
     *
     * @param reads the RDD of reads, in coordinate-sorted order
     * @param variants the RDD of variants
     * @return an RDD that contains each read along with the overlapping variants
     */
    public static JavaPairRDD<GATKRead, Iterable<GATKVariant>> join(final JavaRDD<GATKRead> reads, final JavaRDD<GATKVariant> variants) {
        final JavaSparkContext ctx = new JavaSparkContext(reads.context());
        final Broadcast<IntervalsSkipList<GATKVariant>> variantsBroadcast = ctx.broadcast(new IntervalsSkipList<>(variants.collect()));
        return reads.mapToPair(r -> getOverlapping(r, variantsBroadcast.getValue()));
    }

    /**
     * Joins each read of an RDD<GATKRead> with overlapping variants from an RDD of GATKVariants. Can be used for any size of
     * variants (although they are still read into memory) since Spark broadcast is not used.
     *
     * @param reads the RDD of reads, in coordinate-sorted order
     * @param variantsPaths the path to the variants file
     * @return an RDD that contains each read along with the overlapping variants
     */
    public static JavaPairRDD<GATKRead, Iterable<GATKVariant>> join(final JavaRDD<GATKRead> reads, final List<String> variantsPaths) {
        return reads.mapToPair(r -> getOverlapping(r, KnownSitesCache.getVariants(variantsPaths)));
    }

    private static Tuple2<GATKRead, Iterable<GATKVariant>> getOverlapping(final GATKRead read, final IntervalsSkipList<GATKVariant> intervalsSkipList) {
        if (SimpleInterval.isValid(read.getContig(), read.getStart(), read.getEnd())) {
            return new Tuple2<>(read, intervalsSkipList.getOverlapping(new SimpleInterval(read)));
        } else {
            //Sometimes we have reads that do not form valid intervals (reads that do not consume any ref bases, eg CIGAR 61S90I
            //In those cases, we'll just say that nothing overlaps the read
            return new Tuple2<>(read, Collections.emptyList());
        }
    }
}