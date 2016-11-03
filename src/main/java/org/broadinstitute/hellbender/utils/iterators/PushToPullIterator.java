package org.broadinstitute.hellbender.utils.iterators;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collection;
import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Iterator wrapper around our generic {@link PushPullTransformer)} interface. Wraps an iterator of input elements,
 * and transforms the reads from that iterator using the provided transformer.
 *
 * Converts the push-style {@link PushPullTransformer)} interface to a pull model.
 *
 * @param <T> type of items to be iterated over
 */
public class PushToPullIterator<T> implements Iterator<T>, Iterable<T> {

    private final Iterator<T> inputElements;
    private final PushPullTransformer<T> transformer;
    private Iterator<T> cachedElements = null;
    private T nextElement = null;

    /**
     * @param inputElements wrapped iterator from which this iterator will pull elements
     * @param transformer transformer through which the reads from the wrapped iterator will be fed
     */
    public PushToPullIterator(Iterator<T> inputElements, PushPullTransformer<T> transformer ) {
        Utils.nonNull(inputElements, "iterator must not be null");
        Utils.nonNull(transformer, "transformer must not be null");

        this.inputElements = inputElements;
        this.transformer = transformer;

        advanceToNextElement();
    }

    @Override
    public boolean hasNext() {
        return nextElement != null;
    }

    @Override
    public T next() {
        if ( nextElement == null ) {
            throw new NoSuchElementException("next() called when there are no more items");
        }

        final T toReturn = nextElement;
        advanceToNextElement();

        return toReturn;
    }

    private void advanceToNextElement() {
        if ( readyToReleaseReads() || fillCache() ) {
            nextElement = cachedElements.next();
        }
        else {
            nextElement = null;
        }
    }

    private boolean readyToReleaseReads() {
        return cachedElements != null && cachedElements.hasNext();
    }

    private boolean fillCache() {
        while ( inputElements.hasNext() && ! transformer.hasFinalizedItems() ) {
            transformer.submit(inputElements.next());
        }

        if ( ! inputElements.hasNext() ) {
            transformer.signalEndOfInput();
        }

        final Collection<T> transformedElements = transformer.consumeFinalizedItems();
        cachedElements = transformedElements.iterator();

        return cachedElements.hasNext();
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("Cannot remove records via a Push");
    }

    @Override
    public Iterator<T> iterator() {
        return this;
    }
}