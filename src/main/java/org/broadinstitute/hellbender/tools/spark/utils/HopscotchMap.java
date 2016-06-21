package org.broadinstitute.hellbender.tools.spark.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;

import java.util.Collection;
import java.util.Map;
import java.util.Objects;
import java.util.function.BiPredicate;
import java.util.function.Function;

/**
 * A uniquely keyed map with O(1) operations.
 * Sadly, it's not a java.util.Map, but it does behave like a java.util.Map's entrySet.
 */
@DefaultSerializer(HopscotchMap.Serializer.class)
public final class HopscotchMap<K, V, T extends Map.Entry<K, V>>  extends HopscotchSet<T> {
    public HopscotchMap() {}
    public HopscotchMap( final int capacity ) { super(capacity); }
    public HopscotchMap( final Collection<? extends T> collection ) { super(collection); }
    protected HopscotchMap( final Kryo kryo, final Input input ) { super(kryo, input); }

    /** in a map, uniqueness is on the key value, so we compare keys to see if we already have an entry for that key */
    @Override
    protected BiPredicate<T, T> entryCollides() { return (t1, t2) -> Objects.equals(t1.getKey(), t2.getKey()); }

    /** getKey returns the key part of a Map.Entry */
    @Override
    protected Function<T, Object> toKey() { return T::getKey; }

    @SuppressWarnings("rawtypes")
    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<HopscotchMap> {
        @Override
        public void write( final Kryo kryo, final Output output, final HopscotchMap hopscotchMap ) {
            hopscotchMap.serialize(kryo, output);
        }

        @Override
        public HopscotchMap read( final Kryo kryo, final Input input, final Class<HopscotchMap> klass ) {
            return new HopscotchMap(kryo, input);
        }
    }
}
