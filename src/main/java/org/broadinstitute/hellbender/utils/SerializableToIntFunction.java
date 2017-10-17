package org.broadinstitute.hellbender.utils;

import java.io.Serializable;
import java.util.function.ToIntFunction;

/**
 * Represents a {@link ToIntFunction} implementations that is {@link Serializable}.
 */
@FunctionalInterface
public interface SerializableToIntFunction<T> extends Serializable, ToIntFunction<T> {

    @Override
    int applyAsInt(T value);
}
