package org.broadinstitute.hellbender.tools.spark.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;

import java.util.Collection;
import java.util.Map;
import java.util.function.BiPredicate;

/**
 * A map that can contain multiple values for a given key, but distinct entries.
 */
@DefaultSerializer(HopscotchUniqueMultiMapSpark.Serializer.class)
public final class HopscotchUniqueMultiMapSpark<K, V, T extends Map.Entry<K, V>>  extends HopscotchMultiMapSpark<K, V, T> {
    public HopscotchUniqueMultiMapSpark() {}
    public HopscotchUniqueMultiMapSpark( final int capacity ) { super(capacity); }
    public HopscotchUniqueMultiMapSpark( final Collection<? extends T> collection ) { super(collection); }
    private HopscotchUniqueMultiMapSpark( final Kryo kryo, final Input input ) { super(kryo, input); }

    /** in a unique multimap, uniqueness is on the entry, so we just compare entries to see if we already have one */
    @Override
    protected BiPredicate<T, T> entryCollides() { return Object::equals; }

    @SuppressWarnings("rawtypes")
    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<HopscotchUniqueMultiMapSpark> {
        @Override
        public void write( final Kryo kryo, final Output output, final HopscotchUniqueMultiMapSpark hopscotchMultiMap ) {
            hopscotchMultiMap.serialize(kryo, output);
        }

        @Override
        public HopscotchUniqueMultiMapSpark read( final Kryo kryo, final Input input, final Class<HopscotchUniqueMultiMapSpark> klass ) {
            return new HopscotchUniqueMultiMapSpark(kryo, input);
        }
    }
}
