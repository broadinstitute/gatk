package org.broadinstitute.hellbender.utils;

import java.io.Serializable;
import java.util.function.BiFunction;

/**
 * Interface for serializable {@link BiFunction BiFunctions}.
 */
@FunctionalInterface
public interface SerializableBiFunction<T, U, R> extends BiFunction<T, U, R>, Serializable {

    @Override
    R apply(T t, U u);
}
