package org.broadinstitute.hellbender.utils.iterators;

import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.function.Predicate;

/**
 * An iterator that filters reads from an existing iterator of reads.
 */
public class ReadFilteringIterator extends FilteringIterator<GATKRead> {

    /**
     * Create a ReadFilteringIterator given a pre-existing iterator of reads and a read filter.
     * Only reads that pass the filter will be returned from this iterator.
     *
     * @param nestedIterator underlying iterator from which to pull reads (may not be null)
     * @param readFilter filter to apply to the reads (may not be null)
     */
    public ReadFilteringIterator(Iterator<GATKRead> nestedIterator, ReadFilter readFilter) {
        super(nestedIterator, readFilter);
    }

}
