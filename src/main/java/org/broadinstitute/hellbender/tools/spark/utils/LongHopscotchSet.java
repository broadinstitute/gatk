package org.broadinstitute.hellbender.tools.spark.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.tools.spark.sv.SVUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;
import java.util.NoSuchElementException;
import java.util.Objects;
import java.util.stream.LongStream;

/**
 * This class is based on the HopscotchCollection and HopscotchSet classes for storing Objects. Saves memory used by
 * object headers by using primitive longs.
 * <p>
 * Non-negative longs only! We set the MSB to say that a bin is not empty.
 */
@DefaultSerializer(LongHopscotchSet.Serializer.class)
public final class LongHopscotchSet implements Serializable {

    static final int bytesPerEntry = 9;
    @VisibleForTesting
    static final double LOAD_FACTOR = .85;
    private static final long serialVersionUID = 1L;
    private static final int NO_ELEMENT_INDEX = -1;
    private int capacity;
    private int size;

    // buckets have the most significant bit set to 0 if empty and 1 otherwise (entries must be non-negative)
    private long[] buckets;

    // format of the status bytes:
    // high bit set indicates that the bucket contains a "chain head" (i.e., an entry that naturally belongs in the
    // corresponding bucket).  high bit not set indicates a "squatter" (i.e., an entry that got placed here through the
    // collision resolution methodology).  we use Byte.MIN_VALUE (i.e., 0x80) to pick off this bit.
    // low 7 bits give the (unsigned) offset from the current entry to the next entry in the collision resolution chain.
    // if the low 7 bits are 0, then we'd be pointing at ourselves, which is nonsense, so that particular value marks
    // "end of chain" instead.  we use Byte.MAX_VALUE (i.e., 0x7f) to pick off these bits.
    // If the bucket is unused, the corresponding status byte is irrelevant, but is always set to 0.
    private byte[] status;

    /**
     * make a small LongHopscotchSet
     */
    public LongHopscotchSet() {
        this(12000);
    }

    /**
     * make a LongHopscotchSet for a specified minimum capacity
     */
    public LongHopscotchSet(final int capacity) {
        this.capacity = SetSizeUtils.getLegalSizeAbove(capacity);
        this.size = 0;
        this.buckets = new long[this.capacity];
        this.status = new byte[this.capacity];
    }

    /**
     * make a LongHopscotchSet containing an initial array of values
     */
    public LongHopscotchSet(final long[] values) {
        this.capacity = SetSizeUtils.getLegalSizeAbove(values.length);
        this.size = 0;
        this.buckets = new long[this.capacity];
        this.status = new byte[this.capacity];
        for (final long val : values) {
            add(val);
        }
    }

    protected LongHopscotchSet(final Kryo kryo, final Input stream) {
        capacity = stream.readInt();
        size = 0;
        buckets = new long[capacity];
        status = new byte[capacity];
        int nElements = stream.readInt();
        while (nElements-- > 0) {
            add(stream.readLong());
        }
    }

    /**
     * Returns the value for a given bucket entry (zeroing the bit used for declaring the bucket null)
     */
    private static long getValue(final long entry) {
        return entry & Long.MAX_VALUE;
    }

    public static int longHash(final long entryVal) {
        return (int) SVUtils.fnvLong64(entryVal);
    }

    protected void serialize(final Kryo kryo, final Output stream) {

        stream.writeInt(capacity);
        stream.writeInt(size);

        // write the chain heads, and then the squatters
        int count = 0;
        for (int idx = 0; idx != capacity; ++idx) {
            if (isChainHead(idx)) {
                stream.writeLong(getValue(buckets[idx]));
                count += 1;
            }
        }
        for (int idx = 0; idx != capacity; ++idx) {
            final long val = buckets[idx];
            if (!isUnusedValue(val) && !isChainHead(idx)) {
                stream.writeLong(getValue(val));
                count += 1;
            }
        }
        if (count != size) {
            throw new IllegalStateException("Failed to serialize the expected number of objects: expected=" + size + " actual=" + count + ".");
        }
    }

    public final boolean add(final long entryValue) {
        final int hashValue = longHash(entryValue);
        return add(entryValue, hashValue);
    }

    public final boolean add(final long entryValue, final int hashValue) {
        Utils.validateArg(isValidKey(entryValue), "Tried to add negative entry to LongHopScotchSet");
        if (size == capacity) resize();
        try {
            return insert(entryValue, hashValue);
        } catch (final IllegalStateException ise) {
            resize();
            return insert(entryValue, hashValue);
        }
    }

    public final void clear() {
        for (int idx = 0; idx != capacity; ++idx) {
            buckets[idx] = 0;
            status[idx] = 0;
        }
        size = 0;
    }

    /**
     * maximum number of elements that can be held without resizing. (but we may have to resize earlier.)
     */
    public final long capacity() {
        return capacity;
    }

    public final boolean contains(final long key) {
        return contains(key, longHash(key));
    }

    public final boolean contains(final long key, final int hash) {
        int bucketIndex = hashToIndex(hash);
        if (!isChainHead(bucketIndex)) return false;
        long entryVal = getValue(buckets[bucketIndex]);
        if (entryVal == key) return true;
        int offset;
        while ((offset = getOffset(bucketIndex)) != 0) {
            bucketIndex = getIndex(bucketIndex, offset);
            entryVal = getValue(buckets[bucketIndex]);
            if (entryVal == key) return true;
        }
        return false;
    }

    public final boolean isEmpty() {
        return size == 0;
    }

    // -------- internal methods ----------

    public final LongIterator iterator() {
        return new LongCompleteIterator();
    }

    public final boolean remove(final long key) {
        return remove(key, longHash(key));
    }

    public final boolean remove(final long key, final int hash) {
        Utils.validateArg(isValidKey(key), "Tried to remove by negative key in LongHopScotchSet");
        int bucketIndex = hashToIndex(hash);
        if (isUnusedValue(buckets[bucketIndex]) || !isChainHead(bucketIndex)) return false;
        int predecessorIndex = NO_ELEMENT_INDEX;
        while (getValue(buckets[bucketIndex]) != key) {
            final int offset = getOffset(bucketIndex);
            if (offset == 0) return false;
            predecessorIndex = bucketIndex;
            bucketIndex = getIndex(bucketIndex, offset);
        }
        removeAtIndex(bucketIndex, predecessorIndex);
        return true;
    }

    public final int size() {
        return size;
    }

    private int valueToIndex(final long entryVal) {
        return hashToIndex(longHash(entryVal));
    }

    private int hashToIndex(final int hashVal) {
        int result = hashVal % capacity;
        if (result < 0) result += capacity;
        return result;
    }

    private boolean isValidKey(final long key) {
        return key >= 0;
    }

    private boolean insert(final long entryValue) {
        return insert(entryValue, longHash(entryValue));
    }

    private boolean insert(final long entryValue, final int hashValue) {
        final int bucketIndex = hashToIndex(hashValue);

        // if there's a squatter where the new entry should go, move it elsewhere and put the entry there
        if (!isUnusedValue(buckets[bucketIndex]) && !isChainHead(bucketIndex)) evict(bucketIndex);

        // if the place where it should go is empty, just put the new entry there
        if (isUnusedValue(buckets[bucketIndex])) {
            buckets[bucketIndex] = entryValue;
            setOccupied(bucketIndex);
            status[bucketIndex] = Byte.MIN_VALUE;
            size += 1;
            return true;
        }

        // walk to end of chain
        // along the way, make sure the entry isn't already present if necessary
        int endOfChainIndex = bucketIndex;
        while (true) {
            // if entry is already in the set
            final long tableEntry = buckets[endOfChainIndex];
            if (getValue(tableEntry) == entryValue) return false;
            final int offset = getOffset(endOfChainIndex);
            if (offset == 0) break;
            endOfChainIndex = getIndex(endOfChainIndex, offset);
        }

        // find a place for the new entry
        final int emptyBucketIndex = insertIntoChain(bucketIndex, endOfChainIndex);

        // put the new entry into the empty bucket
        buckets[emptyBucketIndex] = entryValue;
        setOccupied(emptyBucketIndex);
        size += 1;
        return true;
    }

    private void removeAtIndex(final int bucketIndex, final int predecessorIndex) {
        final int offset = getOffset(bucketIndex);
        if (offset == 0) { // if end of chain
            buckets[bucketIndex] = 0;
            status[bucketIndex] = 0;
            if (predecessorIndex != NO_ELEMENT_INDEX) { // fix up offset of previous element in chain if there is one
                status[predecessorIndex] -= getOffset(predecessorIndex);
            }
        } else {
            // move the item at the end of the chain into the hole we're creating by deleting this entry
            int prevIndex = bucketIndex;
            int nextIndex = getIndex(prevIndex, offset);
            int offsetToNext;
            while ((offsetToNext = getOffset(nextIndex)) != 0) {
                prevIndex = nextIndex;
                nextIndex = getIndex(nextIndex, offsetToNext);
            }
            buckets[bucketIndex] = buckets[nextIndex];
            buckets[nextIndex] = 0;
            status[prevIndex] -= getOffset(prevIndex);
        }
        size -= 1;
    }

    private int insertIntoChain(final int bucketIndex, final int endOfChainIndex) {
        final int offsetToEndOfChain = getIndexDiff(bucketIndex, endOfChainIndex);

        // find an empty bucket for the new entry
        int emptyBucketIndex = findEmptyBucket(bucketIndex);

        // if the distance to the empty bucket is larger than this, we'll have to hopscotch
        final int maxOffset = offsetToEndOfChain + Byte.MAX_VALUE;

        // hopscotch the empty bucket into range if it's too far away
        int offsetToEmpty;
        while ((offsetToEmpty = getIndexDiff(bucketIndex, emptyBucketIndex)) > maxOffset) {
            emptyBucketIndex = hopscotch(bucketIndex, emptyBucketIndex);
        }

        // if the new entry lies downstream of the current chain end, just link it in
        if (offsetToEmpty > offsetToEndOfChain) {
            status[endOfChainIndex] += offsetToEmpty - offsetToEndOfChain;
        } else {
            linkIntoChain(bucketIndex, emptyBucketIndex);
        }

        return emptyBucketIndex;
    }

    // walk the chain until we find where the new slot gets linked in
    private void linkIntoChain(final int bucketIndex, final int emptyBucketIndex) {
        int offsetToEmpty = getIndexDiff(bucketIndex, emptyBucketIndex);
        int tmpIndex = bucketIndex;
        int offset;
        while ((offset = getOffset(tmpIndex)) < offsetToEmpty) {
            tmpIndex = getIndex(tmpIndex, offset);
            offsetToEmpty -= offset;
        }
        offset -= offsetToEmpty;
        status[tmpIndex] -= offset;
        status[emptyBucketIndex] = (byte) offset;
    }

    private void evict(final int bucketToEvictIndex) {
        final int bucketIndex = valueToIndex(getValue(buckets[bucketToEvictIndex]));
        final int offsetToEvictee = getIndexDiff(bucketIndex, bucketToEvictIndex);
        int emptyBucketIndex = findEmptyBucket(bucketIndex);
        int fromIndex = bucketIndex;
        while (true) {
            while (getIndexDiff(bucketIndex, emptyBucketIndex) > offsetToEvictee) {
                emptyBucketIndex = hopscotch(fromIndex, emptyBucketIndex);
            }
            if (emptyBucketIndex == bucketToEvictIndex) return;
            fromIndex = emptyBucketIndex;
            linkIntoChain(bucketIndex, emptyBucketIndex);
            int prevIndex = bucketIndex;
            int offsetToNext = getOffset(prevIndex);
            int nextIndex = getIndex(prevIndex, offsetToNext);
            while ((offsetToNext = getOffset(nextIndex)) != 0) {
                prevIndex = nextIndex;
                nextIndex = getIndex(nextIndex, offsetToNext);
            }
            buckets[emptyBucketIndex] = buckets[nextIndex];
            buckets[nextIndex] = 0;
            status[nextIndex] = 0;
            status[prevIndex] -= getOffset(prevIndex);
            emptyBucketIndex = nextIndex;
        }
    }

    private int findEmptyBucket(int bucketIndex) {
        do {
            bucketIndex = getIndex(bucketIndex, 1);
        }
        while (!isUnusedValue(buckets[bucketIndex]));
        return bucketIndex;
    }

    private boolean isChainHead(final int bucketIndex) {
        return (status[bucketIndex] & Byte.MIN_VALUE) != 0;
    }

    private int getOffset(final int bucketIndex) {
        return status[bucketIndex] & Byte.MAX_VALUE;
    }

    private boolean isUnusedValue(final long val) {
        return val == 0L;
    }

    private void setOccupied(final int bucketIndex) {
        buckets[bucketIndex] |= Long.MIN_VALUE;
    }

    private int getIndex(final int bucketIndex, final int offset) {
        int result = bucketIndex + offset;
        if (result >= capacity) result -= capacity;
        else if (result < 0) result += capacity;
        return result;
    }

    // bucket1 is assumed to be upstream of bucket2 (even if bucket2's index has wrapped)
    // i.e., the result is always positive
    private int getIndexDiff(final int bucketIndex1, final int bucketIndex2) {
        int result = bucketIndex2 - bucketIndex1;
        if (result < 0) result += capacity;
        return result;
    }

    private int hopscotch(final int fromIndex, final int emptyBucketIndex) {
        final int fromToEmptyDistance = getIndexDiff(fromIndex, emptyBucketIndex);
        int offsetToEmpty = Byte.MAX_VALUE;
        while (offsetToEmpty > 1) {
            final int bucketIndex = getIndex(emptyBucketIndex, -offsetToEmpty);
            final int offsetInBucket = getOffset(bucketIndex);
            if (offsetInBucket != 0 &&
                    offsetInBucket < offsetToEmpty &&
                    offsetToEmpty - offsetInBucket < fromToEmptyDistance) {
                final int bucketToMoveIndex = getIndex(bucketIndex, offsetInBucket);
                move(bucketIndex, bucketToMoveIndex, emptyBucketIndex);
                return bucketToMoveIndex;
            }
            offsetToEmpty -= 1;
        }
        // this happens now and then, but is usually caught and remedied by a resize
        throw new IllegalStateException("Hopscotching failed at load factor " + (1. * size / capacity));
    }

    private void move(int predecessorBucketIndex, final int bucketToMoveIndex, final int emptyBucketIndex) {
        int toEmptyDistance = getIndexDiff(bucketToMoveIndex, emptyBucketIndex);
        int nextOffset = getOffset(bucketToMoveIndex);
        if (nextOffset == 0 || nextOffset > toEmptyDistance) {
            status[predecessorBucketIndex] += toEmptyDistance;
        } else {
            status[predecessorBucketIndex] += nextOffset;
            toEmptyDistance -= nextOffset;
            predecessorBucketIndex = getIndex(bucketToMoveIndex, nextOffset);
            while ((nextOffset = getOffset(predecessorBucketIndex)) != 0 && nextOffset < toEmptyDistance) {
                toEmptyDistance -= nextOffset;
                predecessorBucketIndex = getIndex(predecessorBucketIndex, nextOffset);
            }
            status[predecessorBucketIndex] = (byte) toEmptyDistance;
        }
        if (nextOffset != 0) {
            status[emptyBucketIndex] = (byte) (nextOffset - toEmptyDistance);
        }
        buckets[emptyBucketIndex] = buckets[bucketToMoveIndex];
        buckets[bucketToMoveIndex] = 0;
        status[bucketToMoveIndex] = 0;
    }

    private void resize() {
        if (buckets == null) {
            throw new IllegalStateException("Someone must be doing something ugly with reflection -- I have no buckets.");
        }
        final int oldCapacity = capacity;
        final int oldSize = size;
        final long[] oldBuckets = buckets;
        final byte[] oldStatus = status;

        capacity = SetSizeUtils.getLegalSizeAbove(capacity);
        size = 0;
        buckets = new long[capacity];
        status = new byte[capacity];

        try {
            int idx = 0;
            do {
                final long entry = oldBuckets[idx];
                if (!isUnusedValue(entry)) insert(getValue(entry));
            }
            while ((idx = (idx + 127) % oldCapacity) != 0);
        } catch (final IllegalStateException ise) {
            capacity = oldCapacity;
            size = oldSize;
            buckets = oldBuckets;
            status = oldStatus;
            // this shouldn't happen except in the case of really bad hashCode implementations
            throw new IllegalStateException("Hopscotching failed at load factor " + 1. * size / capacity + ", and resizing didn't help.");
        }

        if (size != oldSize) {
            // this should never happen, period.
            throw new IllegalStateException("Lost some elements during resizing.");
        }
    }

    public boolean containsAll(final long[] vals) {
        for (final long val : vals) {
            if (!contains(val))
                return false;
        }
        return true;
    }

    public void addAll(final long[] entryValues) {
        for (final long val : entryValues) {
            add(val);
        }
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final LongHopscotchSet that = (LongHopscotchSet) o;

        if (size != that.size) return false;
        final LongIterator itr = this.iterator();
        while (itr.hasNext()) {
            if (!that.contains(itr.next())) return false;
        }
        return true;
    }

    @Override
    public int hashCode() {
        return LongStream.of(buckets).filter(val -> !isUnusedValue(val)).map(LongHopscotchSet::getValue).mapToInt(Objects::hashCode).sum();
    }

    public boolean removeAll(final LongHopscotchSet collection) {
        boolean result = false;
        final LongIterator itr = collection.iterator();
        while (itr.hasNext()) {
            if (remove(itr.next())) result = true;

        }
        return result;
    }

    public boolean removeAll(final long[] set) {
        boolean result = false;
        for (final long val : set) {
            if (remove(val)) result = true;
        }
        return result;
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<LongHopscotchSet> {
        @Override
        public void write(final Kryo kryo, final Output output, final LongHopscotchSet hopscotchSet) {
            hopscotchSet.serialize(kryo, output);
        }

        @Override
        public LongHopscotchSet read(final Kryo kryo, final Input input, final Class<LongHopscotchSet> klass) {
            return new LongHopscotchSet(kryo, input);
        }
    }

    private abstract class LongBaseIterator implements LongIterator {
        protected int currentIndex;
        protected int prevIndex;
        protected int removeIndex;
        protected int removePrevIndex;

        LongBaseIterator() {
            currentIndex = prevIndex = removeIndex = removePrevIndex = NO_ELEMENT_INDEX;
        }

        @Override
        public void remove() {
            if (removeIndex == NO_ELEMENT_INDEX) throw new IllegalStateException("Remove without next.");

            removeAtIndex(removeIndex, removePrevIndex);

            // If we haven't deleted the end of a chain, we'll now have an unseen element under previousElementIndex.
            // So we need to back up and let the user know about it at the next call to the next method.  If we
            // have deleted an end of chain, then the bucket will be empty and no adjustment needs to be made.
            if (!isUnusedValue(buckets[removeIndex])) {
                currentIndex = removeIndex;
                prevIndex = removePrevIndex;
            }

            // Set state to "invalid to call remove again".
            removeIndex = NO_ELEMENT_INDEX;
        }
    }

    private final class LongCompleteIterator extends LongBaseIterator {
        // Class Invariants:
        //  bucketHeadIndex is a valid bucket head until the iteration is complete.
        //    When iteration is complete it has the value buckets.length.
        //  currentIndex points to the element that the next method will return.
        //    It is always primed and ready to go until iteration is complete.
        //    When iteration is complete its value should be ignored (it's not set to anything in particular).
        //  removeIndex is set by calling next.  It is set to the invalid value following a call to remove,
        //    as well as at the beginning of iteration (before the first call to next).  It remains valid when iteration
        //    is complete, unless and until remove is called.
        private int bucketHeadIndex;

        LongCompleteIterator() {
            bucketHeadIndex = NO_ELEMENT_INDEX;
            nextBucketHead();
        }

        @Override
        public boolean hasNext() {
            return bucketHeadIndex != NO_ELEMENT_INDEX;
        }

        @Override
        public long next() {
            if (!hasNext()) throw new NoSuchElementException("Iterator exhausted.");

            removeIndex = currentIndex;
            removePrevIndex = prevIndex;

            final int offset = getOffset(currentIndex);
            // if we're at the end of a chain, advance to the next bucket
            if (offset == 0) nextBucketHead();
            else { // otherwise step to the next item in the chain
                prevIndex = currentIndex;
                currentIndex = getIndex(currentIndex, offset);
            }

            return getValue(buckets[removeIndex]);
        }

        private void nextBucketHead() {
            while (++bucketHeadIndex < buckets.length) {
                if (isChainHead(bucketHeadIndex)) {
                    currentIndex = bucketHeadIndex;
                    prevIndex = NO_ELEMENT_INDEX;
                    return;
                }
            }
            bucketHeadIndex = NO_ELEMENT_INDEX;
        }
    }
}
