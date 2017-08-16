package org.broadinstitute.hellbender.engine.spark;

import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.ShardBoundary;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.io.Serializable;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * {@link JavaPairRDD} of Shards where the key is the shard interval and the value is
 * the record of interest at each shard.
 * @param <L>
 */
public class ShardRDD<L> implements Serializable {

    private static final long serialVersionUID = 1L;

    private JavaPairRDD<ShardBoundary, List<L>> pairRDD;

    public JavaPairRDD<ShardBoundary, List<L>> toPairRDD() {
        return pairRDD;
    }

    public interface Grouper<V, W> extends Serializable {
        List<Tuple2<V, List<W>>> group(final List<V> left, final List<W> right);
    }

    ShardRDD(final JavaPairRDD<ShardBoundary, List<L>> pairRDD) {
        this.pairRDD = Utils.nonNull(pairRDD);
    }

    public <W> ShardRDD<Tuple2<L, List<W>>> groupRight(final ShardRDD<W> right, final Grouper<L, W> grouper) {
        return new ShardRDD<>(pairRDD.join(right.pairRDD).mapValues(values -> grouper.group(values._1(), values._2())));
    }
}
