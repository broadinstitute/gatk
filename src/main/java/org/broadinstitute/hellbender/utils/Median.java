package org.broadinstitute.hellbender.utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Utility class for calculating median from a data set, potentially limiting the size of data to a
 * fixed amount
 *
 */
public final class Median<T extends Comparable<T>> {
    final List<T> values;
    final int maxValuesToKeep;
    boolean sorted = false;

    public Median() {
        this(Integer.MAX_VALUE);
    }

    public Median(final int maxValuesToKeep) {
        this.maxValuesToKeep = maxValuesToKeep;
        this.values = new ArrayList<T>();
    }

    public boolean isFull() {
        return values.size() >= maxValuesToKeep;
    }

    public int size() {
        return values.size();
    }

    public boolean isEmpty() {
        return values.isEmpty();
    }

    public T getMedian() {
        if ( isEmpty() )
            throw new IllegalStateException("Cannot get median value from empty array");
        return getMedian(null);  // note that value null will never be used
    }

    /**
     * Returns the floor((n + 1) / 2) item from the list of values if the list
     * has values, or defaultValue is the list is empty.
     */
    public T getMedian(final T defaultValue) {
        if ( isEmpty() )
            return defaultValue;

        if ( ! sorted ) {
            sorted = true;
            Collections.<T>sort(values);
        }

        final int offset = (int) Math.floor((values.size() + 1) * 0.5) - 1;
        return values.get(offset);
    }

    public boolean add(final T value) {
        if ( ! isFull() ) {
            sorted = false;
            return values.add(value);
        }
        else
            return false;
    }
}
