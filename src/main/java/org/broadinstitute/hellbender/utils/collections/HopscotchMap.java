package org.broadinstitute.hellbender.utils.collections;

import java.util.Collection;
import java.util.Map;
import java.util.Objects;
import java.util.function.BiPredicate;
import java.util.function.Function;

/**
 * A uniquely keyed map with O(1) operations.
 * Sadly, it's not a java.util.Map, but it does behave like a java.util.Map's entrySet.
 */
public final class HopscotchMap<K, V, T extends Map.Entry<K, V>>  extends HopscotchSet<T> {
    public HopscotchMap() {}
    public HopscotchMap( final int capacity ) { super(capacity); }
    public HopscotchMap( final Collection<? extends T> collection ) { super(collection); }

    /** in a map, uniqueness is on the key value, so we compare keys to see if we already have an entry for that key */
    @Override
    protected BiPredicate<T, T> entryCollides() { return (t1, t2) -> Objects.equals(t1.getKey(), t2.getKey()); }

    /** getKey returns the key part of a Map.Entry */
    @Override
    protected Function<T, Object> toKey() { return T::getKey; }
}
