package org.broadinstitute.hellbender.utils.iterators;

import java.util.Iterator;
import java.util.function.Consumer;

/**
 * An {@link Iterator} that automatically closes a resource when the end of the iteration is
 * reached.
 *
 * @param <E>
 */
public class CloseAtEndIterator<E> implements Iterator<E> {

    private final Iterator<E> iterator;
    private final AutoCloseable closeable;

    public CloseAtEndIterator(Iterator<E> iterator, AutoCloseable closeable) {
        this.iterator = iterator;
        this.closeable = closeable;
    }

    @Override
    public boolean hasNext() {
        boolean hasNext = iterator.hasNext();
        if (!hasNext) {
            try {
                closeable.close();
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
        }
        return hasNext;
    }

    @Override
    public E next() {
        return iterator.next();
    }

    @Override
    public void remove() {
        iterator.remove();
    }

    @Override
    public void forEachRemaining(Consumer<? super E> action) {
        iterator.forEachRemaining(action);
    }
}