package org.broadinstitute.hellbender.utils.iterators;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.transformers.SamRecordTransformer;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * An iterator that transforms read (i.e. implements a function SAMRecord -> SAMRecord) from an existing iterator of reads.
 * This is equivalent to a stream.map(readTransformer) operation.
 */
public class SamRecordTransformingIterator implements Iterator<SAMRecord>, Iterable<SAMRecord> {

    private final Iterator<SAMRecord> nestedIterator;
    private final SamRecordTransformer readTransformer;

    /**
     * Create a ReadFilteringIterator given a pre-existing iterator of reads and a read filter.
     * Only reads that pass the filter will be returned from this iterator.
     *
     * @param nestedIterator underlying iterator from which to pull reads (may not be null)
     * @param readTransformer transformation to apply to the reads (may not be null)
     */
    public SamRecordTransformingIterator(final Iterator<SAMRecord> nestedIterator, final SamRecordTransformer readTransformer) {
        this.nestedIterator = Utils.nonNull(nestedIterator);
        this.readTransformer = Utils.nonNull(readTransformer);
    }

    @Override
    public boolean hasNext() {
        return nestedIterator.hasNext();
    }

    @Override
    public SAMRecord next() {
        if ( ! hasNext() ) {
            throw new NoSuchElementException("Iterator exhausted");
        }

        return readTransformer.apply(nestedIterator.next());
    }

    @Override
    public Iterator<SAMRecord> iterator() {
        return this;
    }
}
