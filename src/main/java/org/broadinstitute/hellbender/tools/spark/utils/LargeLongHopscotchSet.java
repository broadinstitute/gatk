package org.broadinstitute.hellbender.tools.spark.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;
import java.util.*;

/**
 * Set of longs that is larger than the max Java array size ( ~ 2^31 ~ 2 billion) and therefore cannot fit into a
 * single LongHopscotchSet. Maintains partitions of LongHopscotchSets given the number of elements to be added.
 * Bins each entry into a LongHopscotchSet using the entry hash value. Note the number of partitions is based on
 * the size estimate passed to the constructor and does not resize dynamically.
 */
@DefaultSerializer(LargeLongHopscotchSet.Serializer.class)
public final class LargeLongHopscotchSet implements Serializable, QueryableLongSet {

    private static final long serialVersionUID = 1L;
    private final List<LongHopscotchSet> sets;
    private final int numSets;

    public LargeLongHopscotchSet(final long numElements) {
        Utils.validateArg(numElements > 0, "Number of elements must be greater than 0");

        int elementsPerPartition = (int) Math.sqrt(numElements);
        int partitions;
        try {
            partitions = SetSizeUtils.getLegalSizeBelow(elementsPerPartition);
        } catch (final IllegalArgumentException e) {
            //If there were no legal sizes small enough, just use 1 set
            partitions = 1;
        }
        elementsPerPartition = (int) ((numElements / partitions) + 1);

        sets = new ArrayList<>();
        for (int i = 0; i < partitions; i++) {
            sets.add(new LongHopscotchSet(elementsPerPartition));
        }
        numSets = sets.size();
    }

    @SuppressWarnings("unchecked")
    protected LargeLongHopscotchSet(final Kryo kryo, final Input stream) {
        final boolean oldReferences = kryo.getReferences();
        kryo.setReferences(false);

        numSets = stream.readInt();
        sets = new ArrayList<>(numSets);
        for (int i = 0; i < numSets; i++) {
            sets.add(kryo.readObject(stream, LongHopscotchSet.class));
        }

        kryo.setReferences(oldReferences);
    }

    private static int longHash(final long entryVal) {
        return (int) SVUtils.fnvLong64(entryVal);
    }

    protected void serialize(final Kryo kryo, final Output stream) {
        final boolean oldReferences = kryo.getReferences();
        kryo.setReferences(false);

        stream.writeInt(sets.size());
        for (final LongHopscotchSet set : sets) {
            kryo.writeObject(stream, set);
        }

        kryo.setReferences(oldReferences);
    }

    public boolean add(final long entryValue) {
        final int hashValue = longHash(entryValue);
        final int setIndex = setIndexOf(hashValue);
        return sets.get(setIndex).add(entryValue, hashValue);
    }

    public void addAll(final long[] entryValues) {
        for (final long val : entryValues) {
            add(val);
        }
    }

    public long size() {
        long sum = 0;
        for (final LongHopscotchSet s : sets) {
            sum += s.size();
        }
        return sum;
    }

    @VisibleForTesting
    public Collection<LongHopscotchSet> getSets() {
        return sets;
    }

    public long capacity() {
        long sum = 0;
        for (final LongHopscotchSet s : sets) {
            sum += s.capacity();
        }
        return sum;
    }

    public boolean contains(final long key) {
        final int hash = longHash(key);
        return sets.get(setIndexOf(hash)).contains(key, hash);
    }

    public boolean containsAll(final long[] vals) {
        for (final long val : vals) {
            if (!contains(val))
                return false;
        }
        return true;
    }

    public LongIterator iterator() {
        return new LargeLongHopscotchSetIterator();
    }

    public boolean remove(final long key) {
        final int hash = longHash(key);
        return sets.get(setIndexOf(hash)).remove(key, hash);
    }

    public boolean removeAll(final long[] set) {
        boolean result = false;
        for (final long val : set) {
            if (remove(val)) result = true;
        }
        return result;
    }

    public boolean isEmpty() {
        return this.size() == 0;
    }

    private int setIndexOf(final int hash) {
        return Integer.remainderUnsigned(hash, numSets);
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (!(o instanceof LargeLongHopscotchSet)) return false;

        final LargeLongHopscotchSet that = (LargeLongHopscotchSet) o;

        return numSets == that.numSets && sets.equals(that.sets);

    }

    @Override
    public int hashCode() {
        return sets.stream().mapToInt(LongHopscotchSet::hashCode).sum();
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<LargeLongHopscotchSet> {
        @Override
        public void write(final Kryo kryo, final Output output, final LargeLongHopscotchSet hopscotchSet) {
            hopscotchSet.serialize(kryo, output);
        }

        @Override
        public LargeLongHopscotchSet read(final Kryo kryo, final Input input, final Class<LargeLongHopscotchSet> klass) {
            return new LargeLongHopscotchSet(kryo, input);
        }
    }

    private final class LargeLongHopscotchSetIterator implements LongIterator {

        //Iterator over partitions
        private final Iterator<LongHopscotchSet> outerIterator;
        //Iterator over current partition. If null, then no elements are left. The converse is true after hasNext().
        private LongIterator innerIterator;

        public LargeLongHopscotchSetIterator() {
            outerIterator = sets.iterator();
            if (outerIterator.hasNext()) {
                //We have at least 1 partition
                innerIterator = outerIterator.next().iterator();
            } else {
                //Empty set of partitions
                innerIterator = null;
            }
        }

        public boolean hasNext() {
            if (innerIterator == null) return false;
            while (!innerIterator.hasNext()) {
                //While we are at the end of a partition, try to move on
                if (outerIterator.hasNext()) {
                    //A next partition exists, so move to it
                    innerIterator = outerIterator.next().iterator();
                } else {
                    //There are no partitions left, so kill the iterator
                    innerIterator = null;
                    return false;
                }
            }
            return true;
        }

        public long next() {
            if (!hasNext())
                throw new NoSuchElementException("LargeLongHopscotchSetIterator is exhausted.");
            return innerIterator.next();
        }

    }

}
