package org.broadinstitute.hellbender.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * A helper class to maintain a cache of an int to double function defined on n = 0, 1, 2. . .
 * The cache expands when a number is not available.
 * NOTE: this cache is thread safe and it may be accessed from multiple threads.
 */
public abstract class IntToDoubleFunctionCache {
    private static final Logger logger = LogManager.getLogger(IntToDoubleFunctionCache.class);

    private double[] cache = new double[] { };

    protected abstract int maxSize();

    protected abstract double compute(final int n);

    public IntToDoubleFunctionCache() { }

    /**
     * Get the value of the function, expanding the cache as necessary.  The cache only applies to non-negative values.
     * @param i operand
     * @return the value of the cached function at {@code i}
     */
    public double get(final int i) {
        Utils.validateArg(i >= 0, () -> String.format("Cache doesn't apply to negative number %d", i));
        if (i >= cache.length) {
            if (i >= maxSize()) {
                return compute(i);
            }
            final int newCapacity = Math.max(i + 10, 2 * cache.length);
            logger.debug("cache miss " + i + " > " + (cache.length-1) + " expanding to " + newCapacity);
            expandCache(newCapacity);
        }
        /*
           Array lookups are not atomic.  It's possible that the reference to cache could be
           changed between the time the reference is loaded and the data is fetched from the correct
           offset.  However, the value retrieved can't change, and it's guaranteed to be present in the
           old reference by the conditional above.
         */
        return cache[i];
    }

    /**
     * Ensures that the cache contains a value for n.  After completion of expandCache(n),
     * get(n) is guaranteed to return without causing a cache expansion
     * @param newCapacity desired value to be precomputed
     */
    public synchronized void expandCache(final int newCapacity) {
        if (newCapacity < cache.length) {
            //prevents a race condition when multiple threads want to expand the cache at the same time.
            //in that case, one of them will be first to enter the synchronized method expandCache and
            //so the others may end up in this method even if n < cache.length
            return;
        }
        final double[] newCache = new double[newCapacity + 1];
        System.arraycopy(cache, 0, newCache, 0, cache.length);
        for (int i = cache.length; i < newCache.length; i++) {
            newCache[i] = compute(i);
        }
        cache = newCache;
    }

    public int size() {
        return cache.length;
    }

}
