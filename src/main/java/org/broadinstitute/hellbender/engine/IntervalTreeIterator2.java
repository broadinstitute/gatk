package org.broadinstitute.hellbender.engine;

public interface IntervalTreeIterator2<E>  {

    IntervalTree2<E>.Node next();
    IntervalTree2<E>.Node peekNext();
    IntervalTree2<E>.Node peekPrevious();
    IntervalTree2<E>.Node previous();
    IntervalTreeIterator2<E> remove();
    IntervalTreeIterator2<E> seek(final int start, final int end);
    boolean hasNext();
    boolean hasPrevious();
}
