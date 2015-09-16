package org.broadinstitute.hellbender.utils.recalibration;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.LRUCache;

/**
 * The object temporarily held by a read that describes all of its covariates.
 *
 * In essence, this is an array of CovariateValues, but it also has some functionality to deal with the optimizations of the NestedHashMap
 */
public final class ReadCovariates {
    private static final Logger logger = LogManager.getLogger(ReadCovariates.class);

    /**
     * How big should we let the LRU cache grow
     */
    private static final int LRU_CACHE_SIZE = 500;

    /**
     * Use an LRU cache to keep cache of keys (int[][][]) arrays for each read length we've seen.
     * The cache allows us to avoid the expense of recreating these arrays for every read.  The LRU
     * keeps the total number of cached arrays to less than LRU_CACHE_SIZE.
     *
     */
    private static final LRUCache<Integer, int[][][]> keysCache = new LRUCache<>(LRU_CACHE_SIZE);

    /**
     * The keys cache is only valid for a single covariate count.  Normally this will remain constant for the analysis.
     * If running multiple analyses (or the unit test suite), it's necessary to clear the cache.
     */
    public static void clearKeysCache() {
        keysCache.clear();
    }

    /**
     * Our keys, indexed by event type x read length x covariate
     */
    private final int[][][] keys;

    /**
     * The index of the current covariate, used by addCovariate
     */
    private int currentCovariateIndex = 0;

    public ReadCovariates(final int readLength, final int numberOfCovariates) {
        final LRUCache<Integer, int[][][]> cache = keysCache;
        final int[][][] cachedKeys = cache.get(readLength);
        if ( cachedKeys == null ) {
            // There's no cached value for read length so we need to create a new int[][][] array
            if ( logger.isDebugEnabled() ) logger.debug("Keys cache miss for length " + readLength + " cache size " + cache.size());
            keys = new int[EventType.values().length][readLength][numberOfCovariates];
            cache.put(readLength, keys);
        } else {
            keys = cachedKeys;
        }
    }

    public void setCovariateIndex(final int index) {
        currentCovariateIndex = index;
    }

    /**
     * Update the keys for mismatch, insertion, and deletion for the current covariate at read offset
     *
     * NOTE: no checks are performed on the number of covariates, for performance reasons.  If the count increases
     * after the keysCache has been accessed, this method will throw an ArrayIndexOutOfBoundsException.  This currently
     * only occurs in the testing harness, and we don't anticipate that it will become a part of normal runs.
     *
     * @param mismatch the mismatch key value
     * @param insertion the insertion key value
     * @param deletion the deletion key value
     * @param readOffset the read offset, must be >= 0 and <= the read length used to create this ReadCovariates
     */
    public void addCovariate(final int mismatch, final int insertion, final int deletion, final int readOffset) {
        keys[EventType.BASE_SUBSTITUTION.ordinal()][readOffset][currentCovariateIndex] = mismatch;
        keys[EventType.BASE_INSERTION.ordinal()][readOffset][currentCovariateIndex] = insertion;
        keys[EventType.BASE_DELETION.ordinal()][readOffset][currentCovariateIndex] = deletion;
    }

    /**
     * Get the keys for all covariates at read position for error model
     *
     * @param readPosition
     * @param errorModel
     * @return
     */
    public int[] getKeySet(final int readPosition, final EventType errorModel) {
        return keys[errorModel.ordinal()][readPosition];
    }

    public int[][] getKeySet(final EventType errorModel) {
        return keys[errorModel.ordinal()];
    }

    // ----------------------------------------------------------------------
    //
    // routines for testing
    //
    // ----------------------------------------------------------------------

    protected int[][] getMismatchesKeySet() { return getKeySet(EventType.BASE_SUBSTITUTION); }
    protected int[][] getInsertionsKeySet() { return getKeySet(EventType.BASE_INSERTION); }
    protected int[][] getDeletionsKeySet() { return getKeySet(EventType.BASE_DELETION); }

    protected int[] getMismatchesKeySet(final int readPosition) {
        return getKeySet(readPosition, EventType.BASE_SUBSTITUTION);
    }

    protected int[] getInsertionsKeySet(final int readPosition) {
        return getKeySet(readPosition, EventType.BASE_INSERTION);
    }

    protected int[] getDeletionsKeySet(final int readPosition) {
        return getKeySet(readPosition, EventType.BASE_DELETION);
    }
}
