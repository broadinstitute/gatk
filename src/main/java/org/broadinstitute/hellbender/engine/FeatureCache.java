package org.broadinstitute.hellbender.engine;

import htsjdk.tribble.Feature;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.*;

/**
 * FeatureCache: helper class for {@link FeatureDataSource} to manage the cache of Feature records used
 * during query operations initiated via {@link FeatureDataSource#query(org.broadinstitute.hellbender.utils.SimpleInterval)}
 * and/or {@link FeatureDataSource#queryAndPrefetch(org.broadinstitute.hellbender.utils.SimpleInterval)}.
 *
 * Strategy is to pre-fetch a large number of records AFTER each query interval that produces
 * a cache miss. This optimizes for the use case of intervals with gradually increasing start
 * positions, as many subsequent queries will find their records wholly contained in the cache
 * before we have another cache miss. Performance will be poor for random/non-localized access
 * patterns, or intervals with decreasing start positions.
 *
 * Usage:
 * -Test whether each query interval is a cache hit via {@link #cacheHit(org.broadinstitute.hellbender.utils.SimpleInterval)}
 *
 * -If it is a cache hit, trim the cache to the start position of the interval (discarding records that
 *  end before the start of the new interval) via {@link #trimToNewStartPosition(int)}, then retrieve
 *  records up to the desired endpoint using {@link #getCachedFeaturesUpToStopPosition(int)}.
 *
 * -If it is a cache miss, reset the cache using {@link #fill(java.util.Iterator, org.broadinstitute.hellbender.utils.SimpleInterval)}, pre-fetching
 *  a large number of records after the query interval in addition to those actually requested.
 *
 * @param <CACHED_FEATURE> Type of Feature record we are caching
 */
class FeatureCache<CACHED_FEATURE extends Feature> {
    private static final Logger logger = LogManager.getLogger(FeatureCache.class);

    /**
     * Our cache of Features, optimized for insertion/removal at both ends.
     */
    private final Deque<CACHED_FEATURE> cache;

    /**
     * Our cache currently contains Feature records overlapping this interval
     */
    private SimpleInterval cachedInterval;

    /**
     * Number of times we called {@link #cacheHit(SimpleInterval)} and it returned true
     */
    private int numCacheHits = 0;

    /**
     * Number of times we called {@link #cacheHit(SimpleInterval)} and it returned false
     */
    private int numCacheMisses = 0;

    /**
     * Initial capacity of our cache (will grow by doubling if needed)
     */
    private static final int INITIAL_CAPACITY = 1024;

    /**
     * When we trim our cache to a new start position, this is the maximum number of
     * Features we expect to need to place into temporary storage for the duration of
     * the trim operation. Performance only suffers slightly if our estimate is wrong.
     */
    private static final int EXPECTED_MAX_OVERLAPPING_FEATURES_DURING_CACHE_TRIM = 128;

    /**
     * Create an initially-empty FeatureCache with default initial capacity
     */
    public FeatureCache() {
        cache = new ArrayDeque<>(INITIAL_CAPACITY);
    }

    /**
     * Get the name of the contig on which the Features in our cache are located
     *
     * @return the name of the contig on which the Features in our cache are located
     */
    public String getContig() {
        return cachedInterval.getContig();
    }

    /**
     * Get the start position of the interval that all Features in our cache overlap
     *
     * @return the start position of the interval that all Features in our cache overlap
     */
    public int getCacheStart() {
        return cachedInterval.getStart();
    }

    /**
     * Get the stop position of the interval that all Features in our cache overlap
     *
     * @return the stop position of the interval that all Features in our cache overlap
     */
    public int getCacheEnd() {
        return cachedInterval.getEnd();
    }

    /**
     * Does our cache currently contain no Features?
     *
     * @return true if our cache contains no Features, otherwise false
     */
    public boolean isEmpty() {
        return cache.isEmpty();
    }

    /**
     * @return Number of times we called {@link #cacheHit(SimpleInterval)} and it returned true
     */
    public int getNumCacheHits() {
        return numCacheHits;
    }

    /**
     * @return Number of times we called {@link #cacheHit(SimpleInterval)} and it returned false
     */
    public int getNumCacheMisses() {
        return numCacheMisses;
    }

    /**
     * Clear our cache and fill it with the records from the provided iterator, preserving their
     * relative ordering, and update our contig/start/stop to reflect the new interval that all
     * records in our cache overlap.
     *
     * Typically each fill operation should involve significant lookahead beyond the region
     * requested so that future queries will be cache hits.
     *
     * @param featureIter iterator from which to pull Features with which to populate our cache
     *                    (replacing existing cache contents)
     * @param interval all Features from featureIter overlap this interval
     */
    public void fill( final Iterator<CACHED_FEATURE> featureIter, final SimpleInterval interval ) {
        cache.clear();
        while ( featureIter.hasNext() ) {
            cache.add(featureIter.next());
        }

        cachedInterval = interval;
    }

    /**
     * Determines whether all records overlapping the provided interval are already contained in our cache.
     *
     * @param interval the interval to check against the contents of our cache
     * @return true if all records overlapping the provided interval are already contained in our cache, otherwise false
     */
    public boolean cacheHit( final SimpleInterval interval ) {
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
     * Trims the cache to the specified new start position by discarding all records that end before it
     * while preserving relative ordering of records.
     *
     * @param newStart new start position on the current contig to which to trim the cache
     */
    public void trimToNewStartPosition( final int newStart ) {
        if ( newStart > cachedInterval.getEnd() ) {
            throw new GATKException(String.format("BUG: attempted to trim Feature cache to an improper new start position (%d). Cache stop = %d",
                    newStart, cachedInterval.getEnd()));
        }

        List<CACHED_FEATURE> overlappingFeaturesBeforeNewStart = new ArrayList<>(EXPECTED_MAX_OVERLAPPING_FEATURES_DURING_CACHE_TRIM);

        // In order to trim the cache to the new start position, we need to find
        // all Features in the cache that start before the new start position,
        // and discard those that don't overlap the new start while keeping those
        // that do overlap. We can stop once we find a Feature that starts on or
        // after the new start position, since the Features are assumed to be sorted
        // by start position.
        while ( ! cache.isEmpty() && cache.getFirst().getStart() < newStart ) {
            CACHED_FEATURE featureBeforeNewStart = cache.removeFirst();

            if ( featureBeforeNewStart.getEnd() >= newStart ) {
                overlappingFeaturesBeforeNewStart.add(featureBeforeNewStart);
            }
        }

        // Add back the Features that started before the new start but overlapped it
        // in the reverse of the order in which we encountered them so that their original
        // relative ordering in the cache is restored.
        for ( int i = overlappingFeaturesBeforeNewStart.size() - 1; i >= 0; --i ) {
            cache.addFirst(overlappingFeaturesBeforeNewStart.get(i));
        }

        // Record our new start boundary
        cachedInterval = new SimpleInterval(cachedInterval.getContig(), newStart, cachedInterval.getEnd());
    }

    /**
     * Returns (but does not remove) all cached Features that overlap the region from the start
     * of our cache (cacheStart) to the specified stop position.
     *
     * @param stopPosition Endpoint of the interval that returned Features must overlap
     * @return all cached Features that overlap the region from the start of our cache to the specified stop position
     */
    public List<CACHED_FEATURE> getCachedFeaturesUpToStopPosition( final int stopPosition ) {
        List<CACHED_FEATURE> matchingFeatures = new ArrayList<>(cache.size());

        // Find (but do not remove from our cache) all Features that start before or on the provided stop position
        for ( CACHED_FEATURE candidateFeature : cache ) {
            if ( candidateFeature.getStart() > stopPosition ) {
                break; // No more possible matches among the remaining cached Features, so stop looking
            }
            matchingFeatures.add(candidateFeature);
        }
        return matchingFeatures;
    }

    /**
     * Print statistics about the cache hit rate for debugging.
     */
    public void printCacheStatistics() {
        printCacheStatistics("");
    }

    /**
     * Print statistics about the cache hit rate for debugging.
     * @param sourceName The source for the features in this cache.
     */
    public void printCacheStatistics(final String sourceName) {

        final String sourceNameString = sourceName.isEmpty() ? "" : "for data source " + sourceName;

        final int totalQueries = getNumCacheHits() + getNumCacheMisses();
        logger.debug(String.format("Cache hit rate %s was %.2f%% (%d out of %d total queries)",
                sourceNameString,
                totalQueries > 0 ? ((double)getNumCacheHits() / totalQueries) * 100.0 : 0.0,
                getNumCacheHits(),
                totalQueries));
    }
}

