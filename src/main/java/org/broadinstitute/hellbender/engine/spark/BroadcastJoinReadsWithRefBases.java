package org.broadinstitute.hellbender.engine.spark;

import com.google.cloud.dataflow.sdk.transforms.SerializableFunction;
import org.apache.spark.SparkContext;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceTwoBitSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import scala.Tuple2;

/**
 * Created by laserson on 9/18/15.
 */
public class BroadcastJoinReadsWithRefBases {
    public static JavaPairRDD<GATKRead, ReferenceBases> addBases(final ReferenceTwoBitSource referenceTwoBitSource,
                                                                 final SerializableFunction<GATKRead, SimpleInterval> windowFunction,
                                                                 final JavaRDD<GATKRead> reads) {
        JavaSparkContext ctx = new JavaSparkContext(reads.context());
        Broadcast<ReferenceTwoBitSource> bReferenceSource = ctx.broadcast(referenceTwoBitSource);
        return reads.mapToPair(read -> {
            SimpleInterval interval = windowFunction.apply(read);
            return new Tuple2<>(read, bReferenceSource.getValue().getReferenceBases(null, interval));
        });
    }
}
