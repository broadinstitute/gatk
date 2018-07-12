package org.broadinstitute.hellbender.tools.spark.utils;

import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.function.Predicate;

/**
 * Applies a predicate to each element of an iterator, returning only those elements that satisfy the predicate.
 */
public final class IteratorFilter<T> implements Iterator<T> {
    private final Iterator<T> iterator;
    private final Predicate<T> predicate;
    private T next;

    public IteratorFilter( final Iterator<T> iterator, final Predicate<T> predicate ) {
        this.iterator = iterator;
        this.predicate = predicate;
    }

    @Override
    public boolean hasNext() {
        if ( next != null ) return true;
        while ( iterator.hasNext() ) {
            next = iterator.next();
            if ( predicate.test(next) ) return true;
        }
        next = null;
        return false;
    }

    @Override
    public T next() {
        if ( !hasNext() ) {
            throw new NoSuchElementException("iterator filter is exhausted.");
        }
        final T result = next;
        next = null;
        return result;
    }

    @Override
    public void remove() {
        if ( next != null ) {
            throw new IllegalStateException("remove without next in iterator filter.");
        }
        iterator.remove();
    }
}
