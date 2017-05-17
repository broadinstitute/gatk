package org.broadinstitute.hellbender.tools.spark.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.tools.spark.sv.SVUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * Set of longs that is larger than the max Java array size ( ~ 2^31 ~ 2 billion) and therefore cannot fit into a
 * single LongHopscotchSet. Maintains partitions of LongHopscotchSets given a maximum partition size and the number
 * of elements to be added. Bins each entry into a corresponding LongHopscotchSet using the entry hash value. Does
 * not support dynamic resizing.
 */
@DefaultSerializer(LargeLongHopscotchSet.Serializer.class)
public class LargeLongHopscotchSet implements Serializable, QueryableLongSet {

    private static final long serialVersionUID = 1L;
    private final List<LongHopscotchSet> sets;
    private int numSets;

    public LargeLongHopscotchSet(final long maxPartitionBytes, final long numElements) {
        sets = new ArrayList<>();
        initialize(maxPartitionBytes, numElements);
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
        return (int) SVUtils.fnvLong64(0x4c79ac6a, entryVal);
    }

    private void initialize(final long maxPartitionBytes, final long numElements) {
        Utils.validateArg(maxPartitionBytes > 0, "Max partition size must be greater than 0");
        Utils.validateArg(numElements > 0, "Number of elements must be greater than 0");

        //Largest partition according to the max bytes per partition constraint
        final int maxPartitionSize = LongHopscotchSet.getLegalSizeBelow(maxPartitionBytes / LongHopscotchSet.bytesPerBucket);
        //Truncate numElements to largest int possible
        final int idealPartitionSize = LongHopscotchSet.getPartitionSize(numElements);
        //Number of elements per partition
        final int elementsPerPartitionGuess = Math.min(idealPartitionSize, maxPartitionSize);
        final int partitions = (int) Math.ceil(((double) numElements) / elementsPerPartitionGuess);
        final int elementsPerPartition = (int)Math.ceil(numElements / (double)partitions);

        for (int i = 0; i < partitions; i++) {
            sets.add(new LongHopscotchSet(elementsPerPartition));
        }
        numSets = sets.size();
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
        final int setIndex = setIndexOf(entryValue);
        final LongHopscotchSet setRef = sets.get(setIndex);
        return setRef.add(entryValue);
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
        return sets.get(setIndexOf(key)).contains(key);
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
        return sets.get(setIndexOf(key)).remove(key);
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

    private int setIndexOf(final long val) {
        return Integer.remainderUnsigned(longHash(val), numSets);
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final LargeLongHopscotchSet that = (LargeLongHopscotchSet) o;

        return numSets == that.numSets && sets.equals(that.sets);

    }

    @Override
    public int hashCode() {
        final LongIterator itr = iterator();
        int sum = 0;
        while (itr.hasNext()) {
            sum += LongHopscotchSet.longHash(itr.next());
        }
        return sum;
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

        private static final int NO_ELEMENT_INDEX = -1;
        private int currentSetIndex;
        private int nextSetIndex;
        private LongIterator itr;

        public LargeLongHopscotchSetIterator() {
            nextSetIndex = numSets > 1 ? 1 : NO_ELEMENT_INDEX;
            if (isEmpty()) {
                currentSetIndex = NO_ELEMENT_INDEX;
            } else {
                currentSetIndex = 0;
                itr = sets.get(currentSetIndex).iterator();
            }
        }

        public boolean hasNext() {
            return itr != null && itr.hasNext();
        }

        public long next() {
            if (!hasNext())
                throw new NoSuchElementException("LargeLongHopscotchSet.LargeLongHopscotchSetIterator is exhausted.");
            final long result = itr.next();
            while (itr != null && !itr.hasNext()) {
                currentSetIndex = nextSetIndex;
                if (currentSetIndex != NO_ELEMENT_INDEX) {
                    itr = sets.get(currentSetIndex).iterator();
                    nextSetIndex = currentSetIndex < numSets - 1 ? nextSetIndex + 1 : NO_ELEMENT_INDEX;
                } else {
                    itr = null;
                }
            }
            return result;
        }

    }

}
