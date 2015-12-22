package org.broadinstitute.hellbender.utils.recalibration.covariates;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.recalibration.EventType;

/**
 * The object temporarily held by a read that describes all of its covariates.
 */
public final class ReadCovariates {
    private static final Logger logger = LogManager.getLogger(ReadCovariates.class);

    /**
     * Our keys, indexed by event type x read length x covariate
     */
    private final int[][][] keys;

    /**
     * The index of the current covariate, used by addCovariate
     */
    private int currentCovariateIndex = 0;

    /**
     * Use an LRU cache to keep cache of keys (int[][][]) arrays for each read length we've seen.
     * The cache allows us to avoid the expense of recreating these arrays for every read.  The LRU
     * keeps the total number of cached arrays to less than LRU_CACHE_SIZE.
     */
    public ReadCovariates(final int readLength, final int numberOfCovariates, final CovariateKeyCache keysCache) {
        Utils.nonNull(keysCache);
        final int[][][] cachedKeys = keysCache.get(readLength);
        if ( cachedKeys == null ) {
            if ( logger.isDebugEnabled() ) logger.debug("Keys cache miss for length " + readLength + " cache size " + keysCache.size());
            keys = new int[EventType.values().length][readLength][numberOfCovariates];
            keysCache.put(readLength, keys);
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
