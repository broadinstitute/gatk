package org.broadinstitute.hellbender.tools.spark.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;

import java.util.Collection;
import java.util.Map;
import java.util.function.BiPredicate;
import java.util.function.Function;

/**
 * A map that can contain multiple values for a given key, but distinct entries.
 */
@DefaultSerializer(HopscotchUniqueMultiMap.Serializer.class)
public final class HopscotchUniqueMultiMap<K, V, T extends Map.Entry<K, V>>  extends HopscotchMultiMap<K, V, T> {
    public HopscotchUniqueMultiMap() {}
    public HopscotchUniqueMultiMap( final int capacity ) { super(capacity); }
    public HopscotchUniqueMultiMap( final Collection<? extends T> collection ) { super(collection); }
    protected HopscotchUniqueMultiMap( final Kryo kryo, final Input input ) { super(kryo, input); }

    /** in a unique multimap, uniqueness is on the entry, so we just compare entries to see if we already have one */
    @Override
    protected BiPredicate<T, T> entryCollides() { return Object::equals; }

    @SuppressWarnings("rawtypes")
    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<HopscotchUniqueMultiMap> {
        @Override
        public void write( final Kryo kryo, final Output output, final HopscotchUniqueMultiMap hopscotchMultiMap ) {
            hopscotchMultiMap.serialize(kryo, output);
        }

        @Override
        public HopscotchUniqueMultiMap read( final Kryo kryo, final Input input, final Class<HopscotchUniqueMultiMap> klass ) {
            return new HopscotchUniqueMultiMap(kryo, input);
        }
    }
}
