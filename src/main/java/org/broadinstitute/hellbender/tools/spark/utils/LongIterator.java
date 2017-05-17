package org.broadinstitute.hellbender.tools.spark.utils;

/**
 * Iterator-like interface for collections of primitive long's
 */
public interface LongIterator {

    boolean hasNext();

    long next();

    default void remove() {
        throw new UnsupportedOperationException("remove");
    }

}