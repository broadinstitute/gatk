package org.broadinstitute.hellbender.tools.spark.utils;

import org.apache.spark.Partitioner;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.Function2;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.ShardBoundary;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipList;
import scala.Serializable;
import scala.Tuple2;

import java.util.*;
import java.util.function.BiPredicate;
import java.util.stream.Collectors;

/**
 * Created by valentin on 8/24/17.
 */
public class ShardedPairRDD<T, U> {

    private final Broadcast<IntervalsSkipList<ShardBoundary>> shards;
    private final JavaPairRDD<ShardBoundary, Tuple2<List<T>, List<U>>> rdd;
    private final Function<? super T, ? extends Iterable<SimpleInterval>> tLoc;
    private final Function<? super U, ? extends Iterable<SimpleInterval>> uLoc;

    ShardedPairRDD(final Broadcast<IntervalsSkipList<ShardBoundary>> shards, final JavaPairRDD<ShardBoundary, Tuple2<List<T>, List<U>>> rdd,
                          final Function<? super T, ? extends Iterable<SimpleInterval>> tLoc,
                          final Function<? super U, ? extends Iterable<SimpleInterval>> uLoc) {
        this.shards = Utils.nonNull(shards);
        this.rdd = Utils.nonNull(rdd);
        this.tLoc = Utils.nonNull(tLoc);
        this.uLoc = Utils.nonNull(uLoc);
    }

    /**
     * Given a mathcing predicate between the {@code <T>} elements and {@code <U>} elments, creates groups
     * of those U that match to every given T.
     * @param match
     * @return
     */
    public JavaPairRDD<T, List<U>> matchLeft(final BiPredicate<? super T, ? super U> match) {
        final BiPredicate<? super T, ? super U> overlapFilter = ( BiPredicate<? super T, ? super U> & Serializable) (t, u) -> match.test(t, u);
        final JavaPairRDD<ShardBoundary, Tuple2<T, List<U>>> x = rdd.flatMapToPair(t -> {
            final List<T> ts = t._2()._1();
            final List<U> us = t._2()._2();
            return ts.stream().map(tt -> new Tuple2<>(tt,
                    us.stream().filter(uu -> overlapFilter.test(tt,uu)).collect(Collectors.toList()))).iterator();
            }).mapToPair(tuple -> {
            final SimpleInterval first = tLoc.call(tuple._1()).iterator().next();
            return new Tuple2<SimpleInterval, Tuple2<T, List<U>>>(
                    new SimpleInterval(first.getContig(), first.getStart(), first.getEnd()), tuple);
            }).mapToPair(tuple -> new Tuple2<>(shards.getValue().getOverlapping(tuple._1()).get(0), tuple._2()));

        final JavaPairRDD<ShardBoundary, List<Tuple2<T, List<U>>>> y = x
                .aggregateByKey(new ArrayList<Tuple2<T, List<U>>>(), (l, tuple) -> { l.add(tuple); return l; }, (l, k) -> { l.addAll(k); return l; })
                .mapValues(list -> {
                    final Map<T, Set<U>> usPerT = new LinkedHashMap<>(list.size());
                    boolean multipleFound = false;
                    for (final Tuple2<T, List<U>> element : list) {
                        if (usPerT.containsKey(element._1())) {
                            multipleFound = true;
                            usPerT.get(element._1()).addAll(element._2());
                        } else {
                            usPerT.put(element._1(), new HashSet<>(element._2()));
                        }
                    }
                    return !multipleFound
                            ? list
                            : usPerT.entrySet().stream()
                                .map(entry -> new Tuple2<>(entry.getKey(), entry.getValue().stream().collect(Collectors.toList())))
                                .collect(Collectors.toList());
                });
        return y.flatMapToPair(tuple -> tuple._2().iterator());
    }
}
