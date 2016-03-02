package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.iterators.ReadFilteringIterator;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Iterator;
import java.util.List;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * A class to represent a window of reads data, optionally expanded by a configurable amount of padded data.
 *
 * The reads are lazily loaded by default (when accessing the reads via {@link #iterator}. Loading all the
 * reads in the window at once is possible via {@link #loadAllReads}.
 *
 * The reads returned will overlap the expanded padded window. It's possible to query whether they are within
 * the main part of the window via {@link #isContainedWithinWindow} and {@link #startsWithinWindow}.
 *
 * The reads in the window can be filtered via {@link #setReadFilter} (no filtering is performed by default).
 */
public final class ReadWindow implements Iterable<GATKRead>, Locatable {

    private final SimpleInterval interval;
    private final SimpleInterval paddedInterval;
    private final ReadsDataSource readsSource;
    private ReadFilter readFilter;

    /**
     * Create a new ReadWindow representing the specified interval, with the specified amount of padding.
     *
     * @param interval the genomic span covered by this window
     * @param paddedInterval the span covered by this window, plus any additional padding on each side (must contain the un-padded interval)
     * @param readsSource source of reads from which to populate this window
     */
    public ReadWindow( final SimpleInterval interval, final SimpleInterval paddedInterval, final ReadsDataSource readsSource ) {
        Utils.nonNull(interval);
        Utils.nonNull(paddedInterval);
        Utils.nonNull(readsSource);

        if ( ! paddedInterval.contains(interval) ) {
            throw new IllegalArgumentException("The padded interval must contain the un-padded interval");
        }

        this.interval = interval;
        this.paddedInterval = paddedInterval;
        this.readsSource = readsSource;
    }

    /**
     * Reads in this window will be filtered using this filter before being returned
     *
     * @param filter filter to use (may be null)
     */
    public void setReadFilter( final ReadFilter filter ) {
        this.readFilter = filter;
    }

    /**
     * @return Contig this window belongs to
     */
    @Override
    public String getContig() {
        return interval.getContig();
    }

    /**
     * @return Start position of this window
     */
    @Override
    public int getStart() {
        return interval.getStart();
    }

    /**
     * @return End position of this window
     */
    @Override
    public int getEnd() {
        return interval.getEnd();
    }

    /**
     * @return the interval this window spans
     */
    public SimpleInterval getInterval() {
        return interval;
    }

    /**
     * @return the interval this window spans, potentially with additional padding on each side
     */
    public SimpleInterval getPaddedInterval() {
        return paddedInterval;
    }

    /**
     * @return an iterator over reads in this window, as filtered using the configured read filter;
     *         reads are lazily loaded rather than pre-loaded
     */
    @Override
    public Iterator<GATKRead> iterator() {
        final Iterator<GATKRead> readsIterator = readsSource.query(paddedInterval);

        return readFilter != null ? new ReadFilteringIterator(readsIterator, readFilter) : readsIterator;
    }

    /**
     * @return a List containing all reads in this window, pre-loaded, and filtered using the configured read filter
     *
     * note: call {@link #iterator} instead to avoid pre-loading all reads at once
     */
    public List<GATKRead> loadAllReads() {
        return StreamSupport.stream(Spliterators.spliteratorUnknownSize(iterator(), Spliterator.ORDERED), false).collect(Collectors.toList());
    }

    /**
     * @param loc Locatable to test
     * @return true if loc is contained within this window's interval, otherwise false
     */
    public boolean isContainedWithinWindow( final Locatable loc ) {
        return interval.contains(loc);
    }

    /**
     * @param loc Locatable to test
     * @return true if loc starts within this window's interval, otherwise false
     */
    public boolean startsWithinWindow( final Locatable loc ) {
        return interval.contains(new SimpleInterval(loc.getContig(), loc.getStart(), loc.getStart()));
    }

    /**
     * @return number of bases of padding to the left of our interval
     */
    public int numLeftPaddingBases() {
        return interval.getStart() - paddedInterval.getStart();
    }

    /**
     * @return number of bases of padding to the right of our interval
     */
    public int numRightPaddingBases() {
        return paddedInterval.getEnd() - interval.getEnd();
    }
}
