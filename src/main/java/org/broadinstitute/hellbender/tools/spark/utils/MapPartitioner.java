package org.broadinstitute.hellbender.tools.spark.utils;

import java.util.Collections;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.function.Function;

/**
 * Class that can be used to implement mapPartition.
 * It requires a Function that transforms a single input element into an Iterator over output elements.
 * The added value is that the individual output Iterators generated from each input element get "glued together".
 *
 * NB: Spark requires this to be an Iterable, not an Iterator, which is inconsistent with functional behavior and
 * arguably a bug.  (Note that mapPartitionsWithIndex expects an Iterator return, instead!)
 * Since there's no way of replacing the input iterator with a fresh one, this is, in fact, traversable only once.
 * Fortunately, Spark traverses only once, so this works out fine.
 */
public class MapPartitioner<I,O> implements Iterable<O>
{
    private final Iterator<I> inputIterator;
    private final Function<I, Iterator<O>> outputIteratorFunction;
    private I sentinel;

    public MapPartitioner( final Iterator<I> inputIterator,
                           final Function<I, Iterator<O>> outputIteratorFunction ) {
        this(inputIterator, outputIteratorFunction, null );
    }

    public MapPartitioner( final Iterator<I> inputIterator,
                           final Function<I, Iterator<O>> outputIteratorFunction,
                           final I sentinel ) {
        this.inputIterator = inputIterator;
        this.outputIteratorFunction = outputIteratorFunction;
        this.sentinel = sentinel;
    }

    public Iterator<O> iterator() { return new FunctionApplicator(); }

    private class FunctionApplicator implements Iterator<O> {
        private Iterator<O> outputIterator = Collections.emptyIterator();

        @Override
        public boolean hasNext()
        {
            while ( !outputIterator.hasNext() ) {
                if ( inputIterator.hasNext() ) outputIterator = outputIteratorFunction.apply(inputIterator.next());
                else if ( sentinel != null ) {
                    outputIterator = outputIteratorFunction.apply(sentinel);
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
    }
}
