package org.broadinstitute.hellbender.engine.cache;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
 * Implementation of a locatable cache with customizable cache strategy.
 *
 * Usage:
 * -Test whether each query interval is a cache hit via {@link #cacheHit(Locatable)}
 *
 * -If it is a cache hit, trim the cache to the start position of the interval (discarding records that
 *  end before the start of the new interval) via {@link # trimCache(Locatable)}, then retrieve
 *  records up to the desired endpoint using {@link #getCachedLocatables(Locatable)}.
 *
 * -If it is a cache miss, reset the cache using {@link #refillQueryCache(Locatable)},
 *  pre-fetching a large number of records after the query interval in addition to those actually requested.
 *
 * @param <T> Type of Locatable record we are caching
 */
public class LocatableCache<T extends Locatable> {
    private static final Logger logger = LogManager.getLogger(LocatableCache.class);

    /**
     * Display name for this cache
     */
    private final String sourceDisplayName;

    /**
     * Our cache of Features, optimized for insertion/removal at both ends.
     */
    private final Deque<T> cache;

    private final LocatableCacheStrategy<T> cachingStrategy;

    /**
     * Our cache currently contains Feature records overlapping this interval
     */
    private Locatable cachedInterval;

    /**
     * Number of times we called {@link #cacheHit(Locatable)} and it returned true
     */
    private int numCacheHits = 0;

    /**
     * Number of times we called {@link #cacheHit(Locatable)} and it returned false
     */
    private int numCacheMisses = 0;

    /**
     * Initial capacity of our cache (will grow by doubling if needed)
     */
    private static final int INITIAL_CAPACITY = 1024;

    /**
     * Create an initially-empty LocatableCache with default initial capacity
     *
     * @param sourceName display name for this cache
     * @param strategy {@link LocatableCacheStrategy} for curating the cache
     */
    public LocatableCache(final String sourceName, final LocatableCacheStrategy<T> strategy) {
        Utils.nonNull(sourceName);
        Utils.nonNull(strategy);

        sourceDisplayName = sourceName;
        cachingStrategy = strategy;
        cache = new ArrayDeque<>(INITIAL_CAPACITY);
    }

    /**
     * Clear our cache and fill it with the records from the provided iterator, preserving their
     * relative ordering, and update our contig/start/stop to reflect the new interval that all
     * records in our cache overlap.
     *
     * Typically each fill operation should involve significant lookahead beyond the region
     * requested so that future queries will be cache hits.
     *
     * @param locatableIter iterator from which to pull Locatables with which to populate our cache
     *                    (replacing existing cache contents)
     * @param interval all Locatables from locatableIter overlap this interval
     */
    private void fill(final Iterator<T> locatableIter, final Locatable interval ) {
        cache.clear();
        while ( locatableIter.hasNext() ) {
            cache.add(locatableIter.next());
        }
        cachedInterval = interval;
    }

    /**
     * Returns a List of all Locatables in this data source that overlap the provided interval.
     *
     * @param interval retrieve all Locatables overlapping this interval
     * @return a {@code List} of all Locatables in this cache that overlap the provided interval
     */
    public List<T> queryAndPrefetch(final Locatable interval) {
        // If the query can be fully satisfied using existing cache contents, prepare for retrieval
        // by discarding all Locatables at the beginning of the cache that end before the start
        // of our query interval.
        if (cacheHit(interval) ) {
            cachedInterval = cachingStrategy.trimCache(cache, cachedInterval, interval);
        }
        // Otherwise, we have at least a partial cache miss, so go to disk to refill our cache.
        else {
            refillQueryCache(interval);
        }

        // Return the subset of our cache that overlaps our query interval
        return getCachedLocatables(interval);
    }

    /**
     * Refill our cache from disk after a cache miss. Will prefetch Locatables overlapping an additional
     * queryLookaheadBases bases after the end of the provided interval, in addition to those overlapping
     * the interval itself.
     *
     * Calling this has the side effect of invalidating (closing) any currently-open iteration over
     * this data source.
     *
     * @param queryInterval the query interval that produced a cache miss
     */
    private void refillQueryCache( final Locatable queryInterval)
    {
        // Expand the end of our query by the configured number of bases, in anticipation of probable future
        // queries with slightly larger start/stop positions.
        //
        // Note that it doesn't matter if we go off the end of the contig in the process, since
        // our reader's query operation is not aware of (and does not care about) contig boundaries.
        // Note: we use addExact to blow up on overflow rather than propagate negative results downstream
        final Locatable cacheInterval = cachingStrategy.getCacheIntervalFromQueryInterval(queryInterval);

        // Query iterator over our reader will be immediately closed after re-populating our cache
        Iterator<T> cacheableIterator = cachingStrategy.refillCache(cacheInterval);
        fill(cacheableIterator, cacheInterval);
    }

    /**
     * Get the name of the contig on which the Locatables in our cache are located
     *
     * @return the name of the contig on which the Locatables in our cache are located
     */
    public String getContig() {
        return cachedInterval.getContig();
    }

    /**
     * Get the start position of the interval that all Locatables in our cache overlap
     *
     * @return the start position of the interval that all Locatables in our cache overlap
     */
    public int getCacheStart() {
        return cachedInterval.getStart();
    }

    /**
     * Get the stop position of the interval that all Locatables in our cache overlap
     *
     * @return the stop position of the interval that all Locatables in our cache overlap
     */
    public int getCacheEnd() {
        return cachedInterval.getEnd();
    }

    /**
     * Does our cache currently contain no Locatables?
     *
     * @return true if our cache contains no Locatables, otherwise false
     */
    public boolean isEmpty() {
        return cache.isEmpty();
    }

    /**
     * @return Number of times we called {@link #cacheHit(Locatable)} and it returned true
     */
    private int getNumCacheHits() {
        return numCacheHits;
    }

    /**
     * @return Number of times we called {@link #cacheHit(Locatable)} and it returned false
     */
    private int getNumCacheMisses() {
        return numCacheMisses;
    }

    /**
     * Determines whether all records overlapping the provided interval are already contained in our cache.
     *
     * @param interval the interval to check against the contents of our cache
     * @return true if all records overlapping the provided interval are already contained in our cache, otherwise false
     */
    @VisibleForTesting
    boolean cacheHit( final Locatable interval ) {
        final boolean cacheHit = cachedInterval != null && cachedInterval.contains(interval);

        if ( cacheHit ) {
            ++numCacheHits;
        }
        else {
            ++numCacheMisses;
        }

        return cacheHit;
    }

    /**
     * Returns (but does not remove) all cached Locatables that overlap the region from the start
     * of our cache (cacheStart) to the specified stop position.
     *
     * @param interval Endpoint of the interval that returned Locatables must overlap
     * @return all cached Locatables that overlap the region from the start of our cache to the specified stop position
     */
    @VisibleForTesting
    List<T> getCachedLocatables(final Locatable interval ) {
        return cachingStrategy.queryCache(cache, interval);
    }

    /**
     * Print statistics about the cache hit rate for debugging.
     */
    public String getCacheStatistics() {

        final int totalQueries = getNumCacheHits() + getNumCacheMisses();
        return String.format("Cache hit rate %s was %.2f%% (%d out of %d total queries)",
                sourceDisplayName,
                totalQueries > 0 ? ((double)getNumCacheHits() / totalQueries) * 100.0 : 0.0,
                getNumCacheHits(),
                totalQueries);
    }
}

