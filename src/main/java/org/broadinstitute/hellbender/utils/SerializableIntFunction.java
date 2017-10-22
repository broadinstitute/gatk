package org.broadinstitute.hellbender.utils;

import java.io.Serializable;
import java.util.function.IntFunction;

/**
 * Created by valentin on 10/19/17.
 */
@FunctionalInterface
public interface SerializableIntFunction<R> extends IntFunction<R>, Serializable {
}
