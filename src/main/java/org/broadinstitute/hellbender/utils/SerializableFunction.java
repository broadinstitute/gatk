package org.broadinstitute.hellbender.utils;

import java.io.Serializable;
import java.util.function.Function;

/**
 * Represents a {@link Function} that is {@link Serializable}.
 */
@FunctionalInterface
public interface SerializableFunction<T,R> extends Function<T, R>, Serializable {

    @Override
    R apply(T t);

}
