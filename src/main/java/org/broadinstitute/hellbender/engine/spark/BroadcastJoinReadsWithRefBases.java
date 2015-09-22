package org.broadinstitute.hellbender.engine.spark;

import com.google.cloud.dataflow.sdk.transforms.SerializableFunction;
import org.apache.spark.SparkContext;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceDataflowSource;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceTwoBitSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import scala.Tuple2;

/**
 * Created by laserson on 9/18/15.
 */
public class BroadcastJoinReadsWithRefBases {
    public static JavaPairRDD<GATKRead, ReferenceBases> addBases(final ReferenceDataflowSource referenceDataflowSource,
                                                                 final JavaRDD<GATKRead> reads) {
        JavaSparkContext ctx = new JavaSparkContext(reads.context());
        Broadcast<ReferenceDataflowSource> bReferenceSource = ctx.broadcast(referenceDataflowSource);
        return reads.mapToPair(read -> {
            SimpleInterval interval = bReferenceSource.getValue().getReferenceWindowFunction().apply(read);
            return new Tuple2<>(read, bReferenceSource.getValue().getReferenceBases(null, interval));
        });
    }

    public static <T> JavaPairRDD<GATKRead, Tuple2<T, ReferenceBases>> addBases(final ReferenceDataflowSource referenceDataflowSource,
                                                                                final JavaPairRDD<GATKRead, T> keyedByRead) {
        JavaSparkContext ctx = new JavaSparkContext(keyedByRead.context());
        Broadcast<ReferenceDataflowSource> bReferenceSource = ctx.broadcast(referenceDataflowSource);
        return keyedByRead.mapToPair(pair -> {
            SimpleInterval interval = bReferenceSource.getValue().getReferenceWindowFunction().apply(pair._1());
            return new Tuple2<>(pair._1(), new Tuple2<>(pair._2(), bReferenceSource.getValue().getReferenceBases(null, interval)));
        });
    }
}
