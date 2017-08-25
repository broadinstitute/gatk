package org.broadinstitute.hellbender.tools.spark.utils;

import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.function.Function;
import org.broadinstitute.hellbender.engine.ShardBoundary;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

/**
 * Created by valentin on 8/24/17.
 */
public class ShardedRDD<T> {

    final JavaPairRDD<ShardBoundary, T> rdd;
    final Function<? super T, ? extends Iterable<SimpleInterval>> toLoc;

    ShardedRDD(final JavaPairRDD<ShardBoundary, T> rdd, final Function<? super T, ? extends Iterable<SimpleInterval>> toLoc) {
        this.rdd = Utils.nonNull(rdd);
        this.toLoc = Utils.nonNull(toLoc);
    }
}
