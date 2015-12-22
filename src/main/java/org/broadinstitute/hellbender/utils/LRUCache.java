package org.broadinstitute.hellbender.utils;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * An LRU cache implemented as an extension to LinkedHashMap
 */
public final class LRUCache<K,V> extends LinkedHashMap<K,V> {

    private static final long serialVersionUID = 1L;
    private final int maxCapacity; // Maximum number of items in the cache.

    public LRUCache(int maxCapacity) {
        super(maxCapacity+1, 1.0f, true); // Pass 'true' for accessOrder.
        this.maxCapacity = maxCapacity;
    }

    @Override
    protected boolean removeEldestEntry(final Map.Entry<K,V> entry) {
        return size() > this.maxCapacity;
    }
}
