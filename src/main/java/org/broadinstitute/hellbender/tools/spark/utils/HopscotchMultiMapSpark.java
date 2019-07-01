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
@DefaultSerializer(HopscotchMultiMapSpark.Serializer.class)
public class HopscotchMultiMapSpark<K, V, T extends Map.Entry<K, V>>  extends HopscotchCollectionSpark<T> {
    public HopscotchMultiMapSpark() {}
    public HopscotchMultiMapSpark( final int capacity ) { super(capacity); }
    public HopscotchMultiMapSpark( final Collection<? extends T> collection ) { super(collection); }
    protected HopscotchMultiMapSpark( final Kryo kryo, final Input input ) { super(kryo, input); }

    /** getKey returns the key part of a Map.Entry */
    @Override
    protected Function<T, Object> toKey() { return T::getKey; }

    @SuppressWarnings("rawtypes")
    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<HopscotchMultiMapSpark> {
        @Override
        public void write(final Kryo kryo, final Output output, final HopscotchMultiMapSpark hopscotchMultiMapSpark ) {
            hopscotchMultiMapSpark.serialize(kryo, output);
        }

        @Override
        public HopscotchMultiMapSpark read( final Kryo kryo, final Input input, final Class<HopscotchMultiMapSpark> klass ) {
            return new HopscotchMultiMapSpark(kryo, input);
        }
    }
}
