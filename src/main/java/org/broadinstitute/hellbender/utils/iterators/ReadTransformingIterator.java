package org.broadinstitute.hellbender.utils.iterators;

import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * An iterator that transforms read (i.e. implements a function GATKRead -> GATKRead) from an existing iterator of reads.
 * This is equivalent to a stream.map(readTransformer) operation.
 */
public class ReadTransformingIterator implements Iterator<GATKRead>, Iterable<GATKRead> {

    private final Iterator<GATKRead> nestedIterator;
    private final ReadTransformer readTransformer;

    /**
     * Create a ReadFilteringIterator given a pre-existing iterator of reads and a read filter.
     * Only reads that pass the filter will be returned from this iterator.
     *
     * @param nestedIterator underlying iterator from which to pull reads (may not be null)
     * @param readTransformer transformation to apply to the reads (may not be null)
     */
    public ReadTransformingIterator(final Iterator<GATKRead> nestedIterator, final ReadTransformer readTransformer) {
        this.nestedIterator = Utils.nonNull(nestedIterator);
        this.readTransformer = Utils.nonNull(readTransformer);
    }

    @Override
    public boolean hasNext() {
        return nestedIterator.hasNext();
    }

    @Override
    public GATKRead next() {
        if ( ! hasNext() ) {
            throw new NoSuchElementException("Iterator exhausted");
        }

        return readTransformer.apply(nestedIterator.next());
    }

    @Override
    public Iterator<GATKRead> iterator() {
        return this;
    }
}
