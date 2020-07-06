package org.broadinstitute.hellbender.utils.iterators;

import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * An iterator that filters reads from an existing iterator of reads.
 */
public class ReadFilteringIterator implements Iterator<GATKRead>, Iterable<GATKRead> {

    private final Iterator<GATKRead> nestedIterator;
    private final ReadFilter readFilter;
    private GATKRead nextRead;

    /**
     * Create a ReadFilteringIterator given a pre-existing iterator of reads and a read filter.
     * Only reads that pass the filter will be returned from this iterator.
     *
     * @param nestedIterator underlying iterator from which to pull reads (may not be null)
     * @param readFilter filter to apply to the reads (may not be null)
     */
    public ReadFilteringIterator( final Iterator<GATKRead> nestedIterator, final ReadFilter readFilter ) {
        Utils.nonNull(nestedIterator);
        Utils.nonNull(readFilter);

        this.nestedIterator = nestedIterator;
        this.readFilter = readFilter;
        this.nextRead = loadNextRead();
    }

    @Override
    public boolean hasNext() {
        return nextRead != null;
    }

    @Override
    public GATKRead next() {
        if ( ! hasNext() ) {
            throw new NoSuchElementException("Iterator exhausted");
        }

        final GATKRead toReturn = nextRead;
        nextRead = loadNextRead();
        return toReturn;
    }

    private GATKRead loadNextRead() {
        while ( nestedIterator.hasNext() ) {
            final GATKRead candidate = nestedIterator.next();
            if ( readFilter.test(candidate) ) {
                return candidate;
            }
        }
        return null;
    }

    @Override
    public Iterator<GATKRead> iterator() {
        return this;
    }
}
