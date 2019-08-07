package org.broadinstitute.hellbender.engine.cache;

import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.ArrayList;
import java.util.Deque;
import java.util.Iterator;
import java.util.List;
import java.util.function.Function;

/**
 * A {@code LocatableCacheStrategy} used to cache primary Feature inputs.
 *
 * Strategy is to pre-fetch a large number of records AFTER each query interval that produces
 * a cache miss. This optimizes for the use case of intervals with gradually increasing start
 * positions, as many subsequent queries will find their records wholly contained in the cache
 * before we have another cache miss. Performance will be poor for random/non-localized access
 * patterns, or intervals with decreasing start positions.
 *
 * @param <CACHED_FEATURE> type Feature being cached
 */
public class DrivingFeatureInputCacheStrategy<CACHED_FEATURE extends Feature> implements LocatableCacheStrategy<CACHED_FEATURE> {

    /**
     * When we trimCache our cache to a new start position, this is the maximum number of
     * {@ocde CACHED_FEATURE} objects we expect to need to place into temporary storage for the duration of
     * the trimCache operation. Performance only suffers slightly if our estimate is wrong.
     */
    private static final int EXPECTED_MAX_OVERLAPPING_FEATURES_DURING_CACHE_TRIM = 128;

    /**
     * When we experience a cache miss (ie., a query interval not fully contained within our cache) and need
     * to re-populate the Feature cache from disk to satisfy a query, this controls the number of extra bases
     * AFTER the end of our interval to fetch. Should be sufficiently large so that typically a significant number
     * of subsequent queries will be cache hits (ie., query intervals fully contained within our cache) before
     * we have another cache miss and need to go to disk again.
     */
    private final int queryLookaheadBases;

    /**
     * Function called on a cache miss to provide results for a cached interval to populate the cache.
     */
    private final Function<SimpleInterval, Iterator<CACHED_FEATURE>> queryResultsProvider;

    /**
     * @param queryLookaheadBases number of a additional base positions beyond the requested query interval to cache
     * @param queryResultsProvider a {@code Function} that can be called get an iterator for a given interval suitable
     *                            for re-populating the cache
     */
    public DrivingFeatureInputCacheStrategy(
            final int queryLookaheadBases,
            final Function<SimpleInterval, Iterator<CACHED_FEATURE>> queryResultsProvider) {
        this.queryLookaheadBases = queryLookaheadBases;
        this.queryResultsProvider = queryResultsProvider;
    }

    /**
     * Given a requested queryInterval, return a new (expanded) interval representing the interval for which items
     * should be cached.
     *
     * @param queryInterval the interval being queried
     * @return the interval to be cached
     */
    @Override
    public SimpleInterval getCacheIntervalFromQueryInterval(final SimpleInterval queryInterval) {
        return new SimpleInterval(queryInterval.getContig(), queryInterval.getStart(), Math.addExact(queryInterval.getEnd(), queryLookaheadBases));
    }

    /**
     * Return a set of features from the cache that satisfy overlap the given interval
     * @param cache the cache from which to pull items
     * @param queryInterval the interval being queried
     * @return {@ocde List} of {@ocde CACHED_FEATURE} objects that overlap {@ocde queryInterval}
     */
    public List<CACHED_FEATURE> queryCache(final Deque<CACHED_FEATURE> cache, final SimpleInterval queryInterval) {
        List<CACHED_FEATURE> matchingFeatures = new ArrayList<>(cache.size());

        // Find (but do not remove from our cache) all Features that start before or on the provided stop position
        for ( CACHED_FEATURE candidateFeature : cache ) {
            if ( candidateFeature.getStart() > queryInterval.getEnd() ) {
                break; // No more possible matches among the remaining cached Features, so stop looking
            }
            matchingFeatures.add(candidateFeature);
        }
        return matchingFeatures;
    }

    /**
     * Refill the cache with items overlapping {@ocde queryInterval}.
     * @param queryInterval the query interval which returned items should overlap
     * @return Iterator of items matching the query interval
     */
    @Override
    public Iterator<CACHED_FEATURE> refillCache(SimpleInterval queryInterval) {
        return queryResultsProvider.apply(queryInterval);
    }

    /**
     * Trims the cache to the specified new start position by discarding all records that end before it
     * while preserving relative ordering of records.
     *
     * @param cache the cache from which to pull items
     * @param cachedInterval the currently cached interval
     * @param interval new interval to which to trim the cache
     * @return the newly cached interval
     */
    @Override
    public SimpleInterval trimCache(final Deque<CACHED_FEATURE> cache, final SimpleInterval cachedInterval, final SimpleInterval interval ) {
        if ( interval.getStart() > cachedInterval.getEnd() ) {
            throw new GATKException(String.format("BUG: attempted to trimCache Feature cache to an improper new start position (%d). Cache stop = %d",
                    interval.getStart(), cachedInterval.getEnd()));
        }

        final List<CACHED_FEATURE> overlappingFeaturesBeforeNewStart = new ArrayList<>(EXPECTED_MAX_OVERLAPPING_FEATURES_DURING_CACHE_TRIM);

        // In order to trimCache the cache to the new start position, we need to find
        // all Features in the cache that start before the new start position,
        // and discard those that don't overlap the new start while keeping those
        // that do overlap. We can stop once we find a Feature that starts on or
        // after the new start position, since the Features are assumed to be sorted
        // by start position.
        while ( ! cache.isEmpty() && cache.getFirst().getStart() < interval.getStart() ) {
            CACHED_FEATURE featureBeforeNewStart = cache.removeFirst();

            if ( featureBeforeNewStart.getEnd() >= interval.getStart() ) {
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
        return new SimpleInterval(cachedInterval.getContig(), interval.getStart(), cachedInterval.getEnd());
    }

}
