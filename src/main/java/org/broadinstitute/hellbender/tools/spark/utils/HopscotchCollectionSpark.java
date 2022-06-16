package org.broadinstitute.hellbender.tools.spark.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.utils.collections.HopscotchCollection;

import java.util.Collection;

@DefaultSerializer(HopscotchCollectionSpark.Serializer.class)
public class HopscotchCollectionSpark<T> extends HopscotchCollection<T> {
    /** make a small HopscotchCollectionSpark */
    public HopscotchCollectionSpark() {}

    /** make a HopscotchCollectionSpark for a specified capacity (or good guess) */
    public HopscotchCollectionSpark( final int capacity ) { super(capacity); }

    /** make a HopscotchCollectionSpark from a collection */
    public HopscotchCollectionSpark( final Collection<? extends T> collection ) { super(collection); }

    @SuppressWarnings("unchecked")
    protected HopscotchCollectionSpark( final Kryo kryo, final Input input ) {
        super((int)(input.readInt() * HopscotchCollection.LOAD_FACTOR));

        final boolean oldReferences = kryo.getReferences();
        kryo.setReferences(false);

        int nElements = input.readInt();
        while ( nElements-- > 0 ) {
            add((T)kryo.readClassAndObject(input));
        }

        kryo.setReferences(oldReferences);
    }

    @SuppressWarnings("unchecked")
    protected void serialize( final Kryo kryo, final Output output ) {
        final boolean oldReferences = kryo.getReferences();
        kryo.setReferences(false);

        final int capacity = capacity();
        output.writeInt(capacity);
        output.writeInt(size());

        // write the chain heads, and then the squatters
        chainHeads().forEach(t -> kryo.writeClassAndObject(output, t));
        squatters().forEach(t -> kryo.writeClassAndObject(output, t));

        kryo.setReferences(oldReferences);
    }

    public static final class Serializer<T> extends com.esotericsoftware.kryo.Serializer<HopscotchCollectionSpark<T>> {
        @Override
        public void write( final Kryo kryo, final Output output, final HopscotchCollectionSpark<T> hopscotchCollection ) {
            hopscotchCollection.serialize(kryo, output);
        }

        @Override
        public HopscotchCollectionSpark<T> read( final Kryo kryo, final Input input,
                                            final Class<HopscotchCollectionSpark<T>> klass ) {
            return new HopscotchCollectionSpark<>(kryo, input);
        }
    }
}
