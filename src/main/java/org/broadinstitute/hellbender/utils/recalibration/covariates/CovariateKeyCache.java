package org.broadinstitute.hellbender.utils.recalibration.covariates;

import it.unimi.dsi.fastutil.ints.Int2ObjectLinkedOpenHashMap;
import org.broadinstitute.hellbender.utils.Utils;

/*
 * Use an LRU cache to keep cache of keys (int[][][]) arrays for each read length we've seen.
 * The cache allows us to avoid the expense of recreating these arrays for every read.  The LRU
 * keeps the total number of cached arrays to less than LRU_CACHE_SIZE.
 * 
 * Uses FastUtil's Int2ObjectLinkedOpenHashMap for primitive int keys to avoid boxing overhead.
 */
public final class CovariateKeyCache {

    /**
     * How big should we let the LRU cache grow
     */
    private static final int LRU_CACHE_SIZE = 500;

    private final Int2ObjectLinkedOpenHashMap<int[][][]> keysCache = new Int2ObjectLinkedOpenHashMap<>();

    /**
     * Get the cached value for the given readlength or null is no value is cached.
     * Uses primitive int key to avoid boxing overhead.
     */
    public int[][][] get(final int readLength) {
        return keysCache.getAndMoveToLast(readLength);
    }

    /**
     * Store the given array in the cache.
     */
    public void put(final int readLength, final int[][][] keys) {
        Utils.nonNull(keys);
        
        // If the key already exists, update it and move to end
        if (keysCache.containsKey(readLength)) {
            keysCache.putAndMoveToLast(readLength, keys);
        } else {
            // If at capacity, remove the first (least recently used) entry
            if (keysCache.size() >= LRU_CACHE_SIZE) {
                keysCache.removeFirst();
            }
            // Add to end (most recently used)
            keysCache.putAndMoveToLast(readLength, keys);
        }
    }

    /**
     * Returns the size of this cache.
     */
    public int size() {
        return keysCache.size();
    }
}
