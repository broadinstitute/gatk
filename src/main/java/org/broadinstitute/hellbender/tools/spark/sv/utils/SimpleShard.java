package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.engine.Shard;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Iterator;

/**
 * Created by valentin on 11/5/17.
 */
@DefaultSerializer(SimpleShard.Serializer.class)
public final class SimpleShard<T> implements Shard<T> {

    private final SimpleInterval interval;
    private Iterable<T> elements;

    private SimpleShard(final SimpleInterval interval, final Iterable<T> elements) {
        this.interval = interval;
        this.elements = elements;
    }

    public static <T>  SimpleShard<T> of(final SimpleInterval interval, final Iterable<T> elements) {
        return new SimpleShard<>(interval, elements);
    }

    @Override
    public SimpleInterval getInterval() {
        return interval;
    }

    @Override
    public SimpleInterval getPaddedInterval() {
        return interval;
    }

    @Override
    public Iterator<T> iterator() {
        return elements.iterator();
    }

    public final static class Serializer<T> extends com.esotericsoftware.kryo.Serializer<SimpleShard<T>> {

        @Override
        public void write(Kryo kryo, Output output, SimpleShard<T> object) {
            output.writeString(object.interval.getContig());
            output.writeInt(object.interval.getStart());
            output.writeInt(object.interval.getEnd());
            kryo.writeClassAndObject(output, object.elements);
        }

        @Override
        public SimpleShard<T> read(Kryo kryo, Input input, Class<SimpleShard<T>> type) {
            final SimpleInterval interval = new SimpleInterval(
                    input.readString(),
                    input.readInt(),
                    input.readInt()
            );
            @SuppressWarnings("unchecked")
            final Iterable<T> elements = (Iterable<T>) kryo.readClassAndObject(input);
            return new SimpleShard<>(interval, elements);
        }
    }
}
