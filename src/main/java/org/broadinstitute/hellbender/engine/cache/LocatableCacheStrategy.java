package org.broadinstitute.hellbender.engine.cache;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Deque;
import java.util.Iterator;
import java.util.List;

/**
 * Interface for implemented by cache strategy objects for {@code LocatableCache}.
 *
 * {@code LocatableCacheStrategy} implementations determine the following policies for the cache:
 *
 * <ul>
 * <li>{@link #getCacheIntervalFromQueryInterval} how to map a requested interval to a (larger) interval to be cached</li>
 * <li>{@link #refillCache}how to (re)populate the cache</li>
 * <li>{@link #queryCache}how to query the cache</li>
 * <li>{@link #trimCache}how to trim the cache</li>
 * </ul>
 */
interface LocatableCacheStrategy<CACHED_FEATURE extends Locatable> {

    /**
     * Given a query interval, return an expanded interval representing the new interval to be cached.
     * @param queryInterval the interval being queried
     * @return the new interval to be cached
     */
    SimpleInterval getCacheIntervalFromQueryInterval(final SimpleInterval queryInterval);

    /**
     * Return a {@code List} of objects from the current cache that overlap the requested interval.
     *
     * The overlapping objects should be returned but not removed from the cache (items are removed
     * by the {@link #trimCache} method).
     *
     * @param cache the cache object
     * @param queryInterval the requested query interval
     * @return {@code List} of cached objects overlapping the query interval
     */
    List<CACHED_FEATURE> queryCache(final Deque<CACHED_FEATURE> cache, final SimpleInterval queryInterval);

    /**
     * Return an iterator of items that overlap {@code queryInterval}, to be used to fill the cache.
     * @param queryInterval the query interval which returned items should overlap
     * @return {@code Iterator} of items overlapping {@code queryInterval}
     */
    Iterator<CACHED_FEATURE> refillCache(final SimpleInterval queryInterval);

    /**
     * Remove items from that are no longer needed from the cache be removing items that overlap the {@code newInterval}.
     *
     * @param cache cache object
     * @param cachedInterval the currenly cached interval
     * @param newInterval the new cached interval being requested
     * @return the new cached interval resulting from the trim operation
     */
    SimpleInterval trimCache(final Deque<CACHED_FEATURE> cache, final SimpleInterval cachedInterval, final SimpleInterval newInterval);
}