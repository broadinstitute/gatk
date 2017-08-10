package org.broadinstitute.hellbender.tools.spark.sv.utils;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.function.Function;
import scala.Tuple2;

import java.io.Serializable;

public final class RDDUtils implements Serializable{
    private static final long serialVersionUID = 1L;

    /**
     * Splits the {@code input} RDD into two based on a filtering predicate.
     * @param input              input RDD to be split
     * @param predicate          filtering criteria
     * @param unpersistInputRDD  option to unpersist the input RDD or not.
     * @return a pair of RDD that where the first has all its elements evaluate to "true" given the predicate,
     *         and the second is the compliment; both are cached.
     */
    public static <T> Tuple2<JavaRDD<T>, JavaRDD<T>> split(final JavaRDD<T> input, final Function<T, Boolean> predicate,
                                                           final boolean unpersistInputRDD) {

        final JavaRDD<T> first = input.filter(predicate);
        final JavaRDD<T> second = input.filter(t -> !predicate.call(t));

        if (unpersistInputRDD) input.unpersist();
        return new Tuple2<>(first, second);
    }
}
