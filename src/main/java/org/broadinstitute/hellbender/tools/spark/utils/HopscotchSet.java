package org.broadinstitute.hellbender.tools.spark.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;

import java.util.Collection;
import java.util.Objects;
import java.util.Set;
import java.util.function.BiPredicate;

/**
 * Implements Set by imposing a unique-element constraint on HopscotchCollection.
 * Also implements equals and hashCode to be consistent with the documented requirements of the java Set interface.
 */
@DefaultSerializer(HopscotchSet.Serializer.class)
public class HopscotchSet<T> extends HopscotchCollection<T> implements Set<T> {
    public HopscotchSet() {}
    public HopscotchSet( final int capacity ) { super(capacity); }
    public HopscotchSet( final Collection<? extends T> collection ) { super(collection); }
    protected HopscotchSet( final Kryo kryo, final Input input ) { super(kryo, input); }

    @Override
    public final boolean equals( final Object obj ) {
        if ( this == obj ) return true;
        if ( !(obj instanceof Set) ) return false;
        @SuppressWarnings("rawtypes")
        final Set that = (Set)obj;
        return this.size() == that.size() && this.containsAll(that);
    }

    @Override
    public final int hashCode() { return stream().mapToInt(Objects::hashCode).sum(); }

    /** in a set, uniqueness is on the entry, so we just compare entries to see if we already have an identical one */
    @Override
    protected BiPredicate<T, T> entryCollides() { return Object::equals; }

    @SuppressWarnings("rawtypes")
    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<HopscotchSet> {
        @Override
        public void write( final Kryo kryo, final Output output, final HopscotchSet hopscotchSet ) {
            hopscotchSet.serialize(kryo, output);
        }

        @Override
        public HopscotchSet read( final Kryo kryo, final Input input, final Class<HopscotchSet> klass ) {
            return new HopscotchSet(kryo, input);
        }
    }
}
