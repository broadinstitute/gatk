package org.broadinstitute.hellbender.engine.spark;

import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import scala.Tuple2;

/**
 * Joins an RDD of GATKReads to reference data using a broadcast strategy.
 *
 * The ReferenceDataflowSource is broadcast using Spark's Broadcast variable mechanism.  The reads are then mapped
 * over and a reference query is executed on each read.  This makes sense for ReferenceDataflowSource implementations
 * that contain the reference data in memory (e.g., ReferenceTwoBitSource), but will likely be much slower for
 * implementations that have to query other resources for the reference sequences.
 */
public class BroadcastJoinReadsWithRefBases {

    /**
     * Joins each read of an RDD<GATKRead> with that read's corresponding reference sequence.
     *
     * @param referenceDataflowSource The source of the reference sequence information
     * @param reads The reads for which to extract reference sequence information
     * @return The JavaPairRDD that contains each read along with the corresponding ReferenceBases object
     */
    public static JavaPairRDD<GATKRead, ReferenceBases> addBases(final ReferenceMultiSparkSource referenceDataflowSource,
                                                                 final JavaRDD<GATKRead> reads) {
        JavaSparkContext ctx = new JavaSparkContext(reads.context());
        Broadcast<ReferenceMultiSparkSource> bReferenceSource = ctx.broadcast(referenceDataflowSource);
        return reads.mapToPair(read -> {
            SimpleInterval interval = bReferenceSource.getValue().getReferenceWindowFunction().apply(read);
            return new Tuple2<>(read, bReferenceSource.getValue().getReferenceBases(interval));
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
        JavaSparkContext ctx = new JavaSparkContext(keyedByRead.context());
        Broadcast<ReferenceMultiSparkSource> bReferenceSource = ctx.broadcast(referenceDataflowSource);
        return keyedByRead.mapToPair(pair -> {
            SimpleInterval interval = bReferenceSource.getValue().getReferenceWindowFunction().apply(pair._1());
            return new Tuple2<>(pair._1(), new Tuple2<>(pair._2(), bReferenceSource.getValue().getReferenceBases(
                    interval)));
        });
    }
}
