package org.broadinstitute.hellbender.tools.spark.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;

import java.util.*;
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

    /**
     * Returns an iterator that is guaranteed to iterate over all entries that share the same key in one consecutive
     * stretch.
     */
    public Iterator<T> getGroupingIterator() {
        return new GroupingIterator();
    }

    protected class GroupingIterator extends CompleteIterator {

        private Set<K> keysInBucket = new HashSet<>();
        private ElementIterator elementIterator = null;

        public GroupingIterator() {
            super();
        }

        @Override
        protected void nextBucketHead() {
            super.nextBucketHead();
            if (bucketHeadIndex != NO_ELEMENT_INDEX) {
                keysInBucket = findKeysInBucket();
                elementIterator = getNextElementIterator();
            }
        }

        @Override
        public boolean hasNext() {
            return ( elementIterator != null && elementIterator.hasNext() );
        }

        @Override
        public T next() {
            final T next = elementIterator.next();
            if (! elementIterator.hasNext()) {
                if (! keysInBucket.isEmpty()) {
                    elementIterator = getNextElementIterator();
                } else {
                    nextBucketHead();
                }
            }
            return next;
        }

        private Set<K> findKeysInBucket() {
            final HashSet<K> keys = new HashSet<>();
            final BucketIterator bucketIterator = new BucketIterator(bucketHeadIndex);
            while (bucketIterator.hasNext()) {
                T next = bucketIterator.next();
                keys.add(next.getKey());
            }
            return keys;
        }

        private ElementIterator getNextElementIterator() {
            if (keysInBucket.isEmpty()) {
                throw new NoSuchElementException("No more keys in bucket.");
            }
            final Iterator<K> keyIterator = keysInBucket.iterator();
            final K nextKey = keyIterator.next();
            keyIterator.remove();
            return new ElementIterator(nextKey);
        }
    }

    private class BucketIterator extends CompleteIterator {

        private BucketIterator(int bucketHeadIndex) {
            this.bucketHeadIndex = bucketHeadIndex;
            currentIndex = bucketHeadIndex;
        }

        @Override
        protected void nextBucketHead() {
            bucketHeadIndex = NO_ELEMENT_INDEX;
        }
    }

}
