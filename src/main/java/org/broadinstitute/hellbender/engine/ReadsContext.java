package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.iterators.ReadFilteringIterator;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;

/**
 * Wrapper around ReadsDataSource that presents reads overlapping a specific interval to a client,
 * without improperly exposing the entire ReadsDataSource interface.
 *
 * Reads data in the interval is lazily queried, so there's no overhead if the client chooses
 * not to examine contextual information from the reads.
 *
 * A ReadsContext may have no backing data source and/or interval. In these cases, queries on it will always
 * return empty iterators. You can determine whether there is a backing source of reads data
 * via {@link #hasBackingDataSource()}, and whether there is an interval via {@link #getInterval}.
 */
public final class ReadsContext implements Iterable<GATKRead> {

    private final GATKDataSource<GATKRead> dataSource;

    private final SimpleInterval interval;

    private final ReadFilter readFilter;

    /**
     * Create an empty ReadsContext with no backing data source or interval. Calls to
     * {@link #iterator} on this context will always return an empty iterator.
     */
    public ReadsContext() {
        this(null, null, null);
    }

    /**
     * Create a ReadsContext backed by the supplied source of reads. Calls to {@link #iterator}
     * will return reads overlapping the provided interval. The data source and/or interval
     * may be null, in which case all calls to {@link #iterator} will return an empty iterator.
     *
     * @param dataSource backing source of reads data (may be null)
     * @param interval interval over which to query (may be null)
     */
    public ReadsContext( final GATKDataSource<GATKRead> dataSource, final SimpleInterval interval ) {
        this(dataSource, interval, null);
    }

    /**
     * Create a ReadsContext backed by the supplied source of reads. Calls to {@link #iterator}
     * will return reads overlapping the provided interval, and passing the provided read filter (if non-null).
     * The data source and/or interval may be null, in which case all calls to {@link #iterator} will return an
     * empty iterator.
     *
     * @param dataSource backing source of reads data (may be null)
     * @param interval interval over which to query (may be null)
     * @param readFilter read filter to be used to filter reads during iteration (may be null)
     */
    public ReadsContext( final GATKDataSource<GATKRead> dataSource, final SimpleInterval interval, final ReadFilter readFilter ) {
        this.dataSource = dataSource;
        this.interval = interval;
        this.readFilter = readFilter;
    }

    /**
     * Does this context have a backing source of reads data?
     *
     * @return true if there is a backing ReadsDataSource, otherwise false
     */
    public boolean hasBackingDataSource() {
        return dataSource != null;
    }

    /**
     * Gets the interval spanned by this context (returned reads will overlap this interval).
     * Null if this context has no interval.
     *
     * @return the interval spanned by this context, or null if this context has no interval
     */
    public SimpleInterval getInterval() {
        return interval;
    }

    /**
     * Get an iterator over the reads in this context. Will return an empty iterator if this
     * context has no backing source of reads and/or no interval.
     *
     * @return iterator over the reads in this context
     */
    @Override
    public Iterator<GATKRead> iterator() {
        return iterator(interval);
    }

    /**
     * Get an iterator over the reads of the backing data source over a given interval. Will return an empty iterator if this
     * context has no backing source of reads.
     *
     * @return iterator over the reads in this context
     */
    public Iterator<GATKRead> iterator(final SimpleInterval interval) {
        // We can't perform a query if we lack either a dataSource or an interval to query on
        if ( dataSource == null || interval == null ) {
            return Collections.<GATKRead>emptyList().iterator();
        }

        return readFilter == null ?
                dataSource.query(interval) :
                new ReadFilteringIterator(dataSource.query(interval), readFilter);
    }

}
