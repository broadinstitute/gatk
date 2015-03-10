package org.broadinstitute.hellbender.utils;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * An LRU cache implemented as an extension to LinkedHashMap
 */
public class LRUCache<K,V> extends LinkedHashMap<K,V> {

    private static final long serialVersionUID = 1L;
    private int capacity; // Maximum number of items in the cache.

    public LRUCache(int capacity) {
        super(capacity+1, 1.0f, true); // Pass 'true' for accessOrder.
        this.capacity = capacity;
    }

    protected boolean removeEldestEntry(final Map.Entry<K,V> entry) {
        return (size() > this.capacity);
    }
}
