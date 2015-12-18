package org.broadinstitute.hellbender.utils.recalibration.covariates;

import org.broadinstitute.hellbender.utils.LRUCache;
import org.broadinstitute.hellbender.utils.Utils;

/*
 * Use an LRU cache to keep cache of keys (int[][][]) arrays for each read length we've seen.
 * The cache allows us to avoid the expense of recreating these arrays for every read.  The LRU
 * keeps the total number of cached arrays to less than LRU_CACHE_SIZE.
 */
public final class CovariateKeyCache {

    /**
     * How big should we let the LRU cache grow
     */
    private static final int LRU_CACHE_SIZE = 500;

    private final LRUCache<Integer, int[][][]> keysCache = new LRUCache<>(LRU_CACHE_SIZE);

    /**
     * Get the cached value for the given readlength or null is no value is cached.
     */
    public int[][][] get(final int readLength) {
        return keysCache.get(readLength);
    }

    /**
     * Store the given array in the cache.
     */
    public void put(final int readLength, final int[][][] keys) {
        Utils.nonNull(keys);
        keysCache.put(readLength, keys);
    }

    /**
     * Returns the size of this cache.
     */
    public int size() {
        return keysCache.size();
    }
}
