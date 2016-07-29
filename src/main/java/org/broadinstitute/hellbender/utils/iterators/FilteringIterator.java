package org.broadinstitute.hellbender.utils.iterators;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.function.Predicate;

/**
 * An iterator that filters objects from an existing iterator of T.
 */
public class FilteringIterator<T> implements Iterator<T>, Iterable<T> {

    private final Iterator<T> nestedIterator;
    private final Predicate<T> filter;
    private T next;

    /**
     * Create a FilteringIterator given a pre-existing iterator of T and a filter.
     * Only objects that pass the filter will be returned from this iterator.
     *
     * @param nestedIterator underlying iterator from which to pull reads (may not be null)
     * @param filter filter to apply to the iterator (may not be null)
     */
    public FilteringIterator( final Iterator<T> nestedIterator, final Predicate<T> filter ) {
        Utils.nonNull(nestedIterator);
        Utils.nonNull(filter);

        this.nestedIterator = nestedIterator;
        this.filter = filter;
        this.next = loadNextRead();
    }

    @Override
    public boolean hasNext() {
        return next != null;
    }

    @Override
    public T next() {
        if ( ! hasNext() ) {
            throw new NoSuchElementException("Iterator exhausted");
        }

        final T toReturn = next;
        next = loadNextRead();
        return toReturn;
    }

    private T loadNextRead() {
        while ( nestedIterator.hasNext() ) {
            final T candidate = nestedIterator.next();
            if ( filter.test(candidate) ) {
                return candidate;
            }
        }
        return null;
    }

    @Override
    public Iterator<T> iterator() {
        return this;
    }

}
