package org.broadinstitute.hellbender.utils.collections;

import java.util.Collection;
import java.util.Map;
import java.util.function.BiPredicate;

/**
 * A map that can contain multiple values for a given key, but distinct entries.
 */
public final class HopscotchUniqueMultiMap<K, V, T extends Map.Entry<K, V>>  extends HopscotchMultiMap<K, V, T> {
    public HopscotchUniqueMultiMap() {}
    public HopscotchUniqueMultiMap( final int capacity ) { super(capacity); }
    public HopscotchUniqueMultiMap( final Collection<? extends T> collection ) { super(collection); }

    /** in a unique multimap, uniqueness is on the entry, so we just compare entries to see if we already have one */
    @Override
    protected BiPredicate<T, T> entryCollides() { return Object::equals; }
}
