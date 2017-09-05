package org.broadinstitute.hellbender.tools.spark.utils;

import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.function.Function;

/**
 *  A little shim that let's you implement a mapPartitions operation (which takes an iterator over all items in the
 *  partition, and returns an iterator over all items to which they are mapped) in terms of a flatMap function (which
 *  takes a single input item, and returns an iterator over any number of output items).
 *  It "glues together" each of the individual output iterators from calling the flatMap function repeatedly, and
 *  presents them as if they were one grand output iterator.
 *  Sometimes it's easier and clearer to write a flatMap function that deals with just one item at a time.
 *  P.S.:  This was written before JavaRDD had flatMap and flatMapToPair methods, so it's utility is not as great
 *  as it used to be.  It survives mostly because of the sentinel functionality, and because the JavaRDD methods
 *  don't offer you the "preservesPartitioning" option.  It can also simplify serialization and makes explicit the
 *  flatMapFunction's lifecycle.
 */
public class FlatMapGluer<I,O> implements Iterator<O> {
    private final Function<I,Iterator<O>> flatMapFunc;
    private final Iterator<? extends I> inputIterator;
    private I sentinel;
    private Iterator<O> outputIterator;

    public FlatMapGluer( final Function<I,Iterator<O>> flatMapFunc,
                         final Iterator<? extends I> inputIterator ) {
        this(flatMapFunc, inputIterator, null);
    }

    public FlatMapGluer( final Function<I,Iterator<O>> flatMapFunc,
                         final Iterator<? extends I> inputIterator,
                         final I sentinel ) {
        this.flatMapFunc = flatMapFunc;
        this.inputIterator = inputIterator;
        this.sentinel = sentinel;
        this.outputIterator = Collections.emptyIterator();
    }

    @Override
    public boolean hasNext()
    {
        while ( !outputIterator.hasNext() ) {
            if ( inputIterator.hasNext() ) outputIterator = flatMapFunc.apply(inputIterator.next());
            else if ( sentinel != null ) {
                outputIterator = flatMapFunc.apply(sentinel);
                sentinel = null;
            }
            else return false;
        }
        return true;
    }

    @Override
    public O next()
    {
        if ( !hasNext() ) throw new NoSuchElementException("Iteration is exhausted.");
        return outputIterator.next();
    }

    public static <I,O> Iterator<O> applyMapFunc(final Function<I,Iterator<O>> flatMapFunc,
                                                 final Iterator<? extends I> inputIterator ) {
        return new FlatMapGluer<>(flatMapFunc,inputIterator,null);
    }

    public static <I,O> Iterator<O> applyMapFunc(final Function<I,Iterator<O>> flatMapFunc,
                                                 final Iterator<? extends I> inputIterator,
                                                 final I sentinel ) {
        return new FlatMapGluer<>(flatMapFunc,inputIterator,sentinel);
    }

    public static <T> Iterator<T> concatIterators( final Iterator<Iterator<T>> itrItr ) {
        return new FlatMapGluer<>(itr -> itr, itrItr);
    }
}
