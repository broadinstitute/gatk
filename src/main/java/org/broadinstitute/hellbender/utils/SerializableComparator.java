package org.broadinstitute.hellbender.utils;

import java.io.Serializable;
import java.util.Comparator;
import java.util.function.Function;

/**
 * Created by valentin on 1/26/18.
 */
public interface SerializableComparator<T> extends Comparator<T>, Serializable {

    @Override
    int compare(T a, T b);
}
