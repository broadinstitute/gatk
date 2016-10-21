package org.broadinstitute.hellbender.tools.spark.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;

import java.util.Collection;
import java.util.Map;
import java.util.function.Function;

/**
 * A map that can contain multiple values for a given key.
 */
@DefaultSerializer(HopscotchMultiMap.Serializer.class)
public class HopscotchMultiMap<K, V, T extends Map.Entry<K, V>>  extends HopscotchCollection<T> {
    public HopscotchMultiMap() {}
    public HopscotchMultiMap( final int capacity ) { super(capacity); }
    public HopscotchMultiMap( final Collection<? extends T> collection ) { super(collection); }
    protected HopscotchMultiMap( final Kryo kryo, final Input input ) { super(kryo, input); }

    /** getKey returns the key part of a Map.Entry */
    @Override
    protected Function<T, Object> toKey() { return T::getKey; }

    @SuppressWarnings("rawtypes")
    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<HopscotchMultiMap> {
        @Override
        public void write(final Kryo kryo, final Output output, final HopscotchMultiMap hopscotchMultiMap ) {
            hopscotchMultiMap.serialize(kryo, output);
        }

        @Override
        public HopscotchMultiMap read(final Kryo kryo, final Input input, final Class<HopscotchMultiMap> klass ) {
            return new HopscotchMultiMap(kryo, input);
        }
    }
}
