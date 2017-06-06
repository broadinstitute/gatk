package org.broadinstitute.hellbender.tools.spark.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;

import java.util.*;
import java.util.function.BiPredicate;
import java.util.function.Function;

/**
 * Multiset implementation that provides low memory overhead with a high load factor by using the hopscotch algorithm.
 * Retrieval times are usually a little slower than for the JDK's HashSet (neglecting GC time), but this class uses much
 * less memory.  Insertion and deletion times are typically a little better than HashSet.
 * It's probably a nice choice for very large collections.
 *
 * You can make it faster than HashSet by replacing the prime-number table sizing with power-of-2 table
 * sizing, but then it'll behave just as badly as HashSet given poor hashCode implementations.  (The extra time in
 * retrieval is mostly due to having to take the hash value mod the capacity, rather than masking off a few bits.)
 *
 * Very rarely, and usually when your element type's hashCode implementation is poor, the hopscotching will fail even
 * after the table is (automatically) resized to be 1.4 times larger.  In this case an IllegalStateException will be
 * thrown.  To be honest, I haven't ever seen it happen with this version of the code, but it might.  If it happens to
 * you, you might want to take a look at your T's hashCode implementation to make certain it has good avalanche
 * characteristics.  The ones you see recommended everywhere that add in the final piece of state, don't have good
 * avalanche.  We try to take care of this common defect with the SPREADER, but that's not a perfect solution.
 */
@DefaultSerializer(HopscotchCollection.Serializer.class)
public class HopscotchCollection<T> extends AbstractCollection<T> {
    // For power-of-2 table sizes add this line
    //private static final int maxCapacity = Integer.MAX_VALUE/2 + 1;

    private int capacity;
    // unused buckets contain null.  (this data structure does not support null entries.)
    // if the bucket is unused, the corresponding status byte is irrelevant, but is always set to 0.
    private int size;
    private T[] buckets;
    // format of the status bytes:
    // high bit set indicates that the bucket contains a "chain head" (i.e., an entry that naturally belongs in the
    // corresponding bucket).  high bit not set indicates a "squatter" (i.e., an entry that got placed here through the
    // collision resolution methodology).  we use Byte.MIN_VALUE (i.e., 0x80) to pick off this bit.
    // low 7 bits give the (unsigned) offset from the current entry to the next entry in the collision resolution chain.
    // if the low 7 bits are 0, then we'd be pointing at ourselves, which is nonsense, so that particular value marks
    // "end of chain" instead.  we use Byte.MAX_VALUE (i.e., 0x7f) to pick off these bits.
    private byte[] status;

    private static final double LOAD_FACTOR = .85;
    private static final int NO_ELEMENT_INDEX = -1;
    private static final int SPREADER = 241;

    /** make a small HopscotchCollection */
    public HopscotchCollection() { this(12000); }

    /** make a HopscotchCollection for a specified capacity (or good guess) */
    @SuppressWarnings("unchecked")
    public HopscotchCollection( final int capacity ) {
        this.capacity = computeCapacity(capacity);
        this.size = 0;
        // Not unsafe, because the remainder of the API allows only elements known to be T's to be assigned to buckets.
        this.buckets = (T[])new Object[this.capacity];
        this.status = new byte[this.capacity];
    }

    /** make a HopscotchCollection from a collection */
    public HopscotchCollection( final Collection<? extends T> collection ) {
        this(collection.size());
        addAll(collection);
    }

    @SuppressWarnings("unchecked")
    protected HopscotchCollection( final Kryo kryo, final Input input ) {
        final boolean oldReferences = kryo.getReferences();
        kryo.setReferences(false);

        capacity = input.readInt();
        size = 0;
        buckets = (T[])new Object[capacity];
        status = new byte[capacity];
        int nElements = input.readInt();
        while ( nElements-- > 0 ) {
            add((T)kryo.readClassAndObject(input));
        }

        kryo.setReferences(oldReferences);
    }

    protected void serialize( final Kryo kryo, final Output output ) {
        final boolean oldReferences = kryo.getReferences();
        kryo.setReferences(false);

        output.writeInt(capacity);
        output.writeInt(size);

        // write the chain heads, and then the squatters
        int count = 0;
        for ( int idx = 0; idx != capacity; ++idx ) {
            if ( isChainHead(idx) ) {
                kryo.writeClassAndObject(output, buckets[idx]);
                count += 1;
            }
        }
        for ( int idx = 0; idx != capacity; ++idx ) {
            final T val = buckets[idx];
            if ( val != null && !isChainHead(idx) ) {
                kryo.writeClassAndObject(output, val);
                count += 1;
            }
        }

        kryo.setReferences(oldReferences);

        if ( count != size ) {
            throw new IllegalStateException("Failed to serialize the expected number of objects: expected="+size+" actual="+count+".");
        }
    }

    @Override
    public final boolean add( final T entry ) {
        if ( entry == null ) throw new UnsupportedOperationException("This collection cannot contain null.");
        if ( size == capacity ) resize();
        final BiPredicate<T, T> collision = entryCollides();
        try {
            return insert(entry, collision);
        } catch ( final IllegalStateException ise ) {
            resize();
            return insert(entry, collision);
        }
    }

    @Override
    public final void clear() {
        for ( int idx = 0; idx != capacity; ++idx ) {
            buckets[idx] = null;
            status[idx] = 0;
        }
        size = 0;
    }

    /** maximum number of elements that can be held without resizing. (but we may have to resize earlier.) */
    public final int capacity() { return capacity; }

    @Override
    public final boolean contains( final Object key ) { return find(key) != null; }

    /** find an entry equivalent to the key, or return null */
    public final T find( final Object key ) {
        int bucketIndex = hashToIndex(key);
        if ( !isChainHead(bucketIndex) ) return null;
        T entry = buckets[bucketIndex];
        if ( equivalent(entry, key) ) return entry;
        int offset;
        while ( (offset = getOffset(bucketIndex)) != 0 ) {
            bucketIndex = getIndex(bucketIndex, offset);
            entry = buckets[bucketIndex];
            if ( equivalent(entry, key) ) return entry;
        }
        return null;
    }

    /** get an iterator over each of the elements equivalent to the key */
    public final Iterator<T> findEach( final Object key ) {
        return new ElementIterator(key);
    }

    @Override
    public final boolean isEmpty() { return size == 0; }

    @Override
    public final Iterator<T> iterator() { return new CompleteIterator(); }

    @Override
    public final boolean remove( final Object key ) {
        int bucketIndex = hashToIndex(key);
        if ( buckets[bucketIndex] == null || !isChainHead(bucketIndex) ) return false;
        int predecessorIndex = NO_ELEMENT_INDEX;
        while ( !equivalent(buckets[bucketIndex], key) ) {
            final int offset = getOffset(bucketIndex);
            if ( offset == 0 ) return false;
            predecessorIndex = bucketIndex;
            bucketIndex = getIndex(bucketIndex, offset);
        }
        removeAtIndex(bucketIndex, predecessorIndex);
        return true;
    }

    /** remove each of the elements equivalent to the key */
    public final boolean removeEach( final Object key ) {
        boolean result = false;
        final Iterator<T> elementItr = new ElementIterator(key);
        while ( elementItr.hasNext() ) {
            elementItr.next();
            elementItr.remove();
            result = true;
        }
        return result;
    }

    /** unlike the AbstractCollection implementation, this one iterates over the supplied collection */
    @Override
    public final boolean removeAll( final Collection<?> collection ) {
        boolean result = false;
        for ( final Object entry : collection ) {
            if ( removeEach(entry) ) result = true;
        }
        return result;
    }

    @Override
    public final int size() { return size; }


    /** in a general multiset there is no uniqueness criterion, so there is never a collision */
    protected BiPredicate<T, T> entryCollides() { return (t1, t2) -> false; }

    /** in a general multiset, the entry is the key */
    protected Function<T, Object> toKey() { return i -> i; }

    // -------- internal methods ----------

    private boolean equivalent( final T entry, final Object key ) {
        return Objects.equals(toKey().apply(entry), key);
    }

    private int hashToIndex( final Object entry ) {
        // For power-of-2 table sizes substitute this line
        // return (SPREADER*hashVal)&(capacity-1);
        int result = entry == null ? 0 : ((SPREADER * entry.hashCode()) % capacity);
        if ( result < 0 ) result += capacity;
        return result;
    }

    private boolean insert( final T entry, final BiPredicate<T, T> collision ) {
        final int bucketIndex = hashToIndex(toKey().apply(entry));

        // if there's a squatter where the new entry should go, move it elsewhere and put the entry there
        if ( buckets[bucketIndex] != null && !isChainHead(bucketIndex) ) evict(bucketIndex);

        // if the place where it should go is empty, just put the new entry there
        if ( buckets[bucketIndex] == null ) {
            buckets[bucketIndex] = entry;
            status[bucketIndex] = Byte.MIN_VALUE;
            size += 1;
            return true;
        }

        // walk to end of chain
        // along the way, make sure the entry isn't already present if necessary
        int endOfChainIndex = bucketIndex;
        while ( true ) {
            // if entry is already in the set
            final T tableEntry = buckets[endOfChainIndex];
            if ( collision.test(tableEntry, entry) ) return false;
            final int offset = getOffset(endOfChainIndex);
            if ( offset == 0 ) break;
            endOfChainIndex = getIndex(endOfChainIndex, offset);
        }

        // find a place for the new entry
        final int emptyBucketIndex = insertIntoChain(bucketIndex, endOfChainIndex);

        // put the new entry into the empty bucket
        buckets[emptyBucketIndex] = entry;
        size += 1;
        return true;
    }

    private void removeAtIndex( final int bucketIndex, final int predecessorIndex ) {
        final int offset = getOffset(bucketIndex);
        if ( offset == 0 ) { // if end of chain
            buckets[bucketIndex] = null;
            status[bucketIndex] = 0;
            if ( predecessorIndex != NO_ELEMENT_INDEX ) { // fix up offset of previous element in chain if there is one
                status[predecessorIndex] -= getOffset(predecessorIndex);
            }
        } else {
            // move the item at the end of the chain into the hole we're creating by deleting this entry
            int prevIndex = bucketIndex;
            int nextIndex = getIndex(prevIndex, offset);
            int offsetToNext;
            while ( (offsetToNext = getOffset(nextIndex)) != 0 ) {
                prevIndex = nextIndex;
                nextIndex = getIndex(nextIndex, offsetToNext);
            }
            buckets[bucketIndex] = buckets[nextIndex];
            buckets[nextIndex] = null;
            status[prevIndex] -= getOffset(prevIndex);
        }
        size -= 1;
    }

    private int insertIntoChain( final int bucketIndex, final int endOfChainIndex ) {
        final int offsetToEndOfChain = getIndexDiff(bucketIndex, endOfChainIndex);

        // find an empty bucket for the new entry
        int emptyBucketIndex = findEmptyBucket(bucketIndex);

        // if the distance to the empty bucket is larger than this, we'll have to hopscotch
        final int maxOffset = offsetToEndOfChain + Byte.MAX_VALUE;

        // hopscotch the empty bucket into range if it's too far away
        int offsetToEmpty;
        while ( (offsetToEmpty = getIndexDiff(bucketIndex, emptyBucketIndex)) > maxOffset ) {
            emptyBucketIndex = hopscotch(bucketIndex, emptyBucketIndex);
        }

        // if the new entry lies downstream of the current chain end, just link it in
        if ( offsetToEmpty > offsetToEndOfChain ) {
            status[endOfChainIndex] += offsetToEmpty - offsetToEndOfChain;
        } else {
            linkIntoChain(bucketIndex, emptyBucketIndex);
        }

        return emptyBucketIndex;
    }

    // walk the chain until we find where the new slot gets linked in
    private void linkIntoChain( final int bucketIndex, final int emptyBucketIndex ) {
        int offsetToEmpty = getIndexDiff(bucketIndex, emptyBucketIndex);
        int tmpIndex = bucketIndex;
        int offset;
        while ( (offset = getOffset(tmpIndex)) < offsetToEmpty ) {
            tmpIndex = getIndex(tmpIndex, offset);
            offsetToEmpty -= offset;
        }
        offset -= offsetToEmpty;
        status[tmpIndex] -= offset;
        status[emptyBucketIndex] = (byte) offset;
    }

    private void evict( final int bucketToEvictIndex ) {
        final int bucketIndex = hashToIndex(toKey().apply(buckets[bucketToEvictIndex]));
        final int offsetToEvictee = getIndexDiff(bucketIndex, bucketToEvictIndex);
        int emptyBucketIndex = findEmptyBucket(bucketIndex);
        int fromIndex = bucketIndex;
        while ( true ) {
            while ( getIndexDiff(bucketIndex, emptyBucketIndex) > offsetToEvictee ) {
                emptyBucketIndex = hopscotch(fromIndex, emptyBucketIndex);
            }
            if ( emptyBucketIndex == bucketToEvictIndex ) return;
            fromIndex = emptyBucketIndex;
            linkIntoChain(bucketIndex, emptyBucketIndex);
            int prevIndex = bucketIndex;
            int offsetToNext = getOffset(prevIndex);
            int nextIndex = getIndex(prevIndex, offsetToNext);
            while ( (offsetToNext = getOffset(nextIndex)) != 0 ) {
                prevIndex = nextIndex;
                nextIndex = getIndex(nextIndex, offsetToNext);
            }
            buckets[emptyBucketIndex] = buckets[nextIndex];
            buckets[nextIndex] = null;
            status[nextIndex] = 0;
            status[prevIndex] -= getOffset(prevIndex);
            emptyBucketIndex = nextIndex;
        }
    }

    private int findEmptyBucket( int bucketIndex ) {
        do {
            bucketIndex = getIndex(bucketIndex, 1);
        }
        while ( buckets[bucketIndex] != null );
        return bucketIndex;
    }

    private boolean isChainHead( final int bucketIndex ) {
        return (status[bucketIndex] & Byte.MIN_VALUE) != 0;
    }

    private int getOffset( final int bucketIndex ) {
        return status[bucketIndex] & Byte.MAX_VALUE;
    }

    private int getIndex( final int bucketIndex, final int offset ) {
        int result = bucketIndex + offset;
        if ( result >= capacity ) result -= capacity;
        else if ( result < 0 ) result += capacity;
        return result;
    }

    // bucket1 is assumed to be upstream of bucket2 (even if bucket2's index has wrapped)
    // i.e., the result is always positive
    private int getIndexDiff( final int bucketIndex1, final int bucketIndex2 ) {
        int result = bucketIndex2 - bucketIndex1;
        if ( result < 0 ) result += capacity;
        return result;
    }

    private int hopscotch( final int fromIndex, final int emptyBucketIndex ) {
        final int fromToEmptyDistance = getIndexDiff(fromIndex, emptyBucketIndex);
        int offsetToEmpty = Byte.MAX_VALUE;
        while ( offsetToEmpty > 1 ) {
            final int bucketIndex = getIndex(emptyBucketIndex, -offsetToEmpty);
            final int offsetInBucket = getOffset(bucketIndex);
            if ( offsetInBucket != 0 &&
                    offsetInBucket < offsetToEmpty &&
                    offsetToEmpty-offsetInBucket < fromToEmptyDistance ) {
                final int bucketToMoveIndex = getIndex(bucketIndex, offsetInBucket);
                move(bucketIndex, bucketToMoveIndex, emptyBucketIndex);
                return bucketToMoveIndex;
            }
            offsetToEmpty -= 1;
        }
        // this happens now and then, but is usually caught and remedied by a resize
        throw new IllegalStateException("Hopscotching failed at load factor "+(1.*size/capacity));
    }

    private void move( int predecessorBucketIndex, final int bucketToMoveIndex, final int emptyBucketIndex ) {
        int toEmptyDistance = getIndexDiff(bucketToMoveIndex, emptyBucketIndex);
        int nextOffset = getOffset(bucketToMoveIndex);
        if ( nextOffset == 0 || nextOffset > toEmptyDistance ) {
            status[predecessorBucketIndex] += toEmptyDistance;
        } else {
            status[predecessorBucketIndex] += nextOffset;
            toEmptyDistance -= nextOffset;
            predecessorBucketIndex = getIndex(bucketToMoveIndex, nextOffset);
            while ( (nextOffset = getOffset(predecessorBucketIndex)) != 0 && nextOffset < toEmptyDistance ) {
                toEmptyDistance -= nextOffset;
                predecessorBucketIndex = getIndex(predecessorBucketIndex, nextOffset);
            }
            status[predecessorBucketIndex] = (byte) toEmptyDistance;
        }
        if ( nextOffset != 0 ) {
            status[emptyBucketIndex] = (byte) (nextOffset - toEmptyDistance);
        }
        buckets[emptyBucketIndex] = buckets[bucketToMoveIndex];
        buckets[bucketToMoveIndex] = null;
        status[bucketToMoveIndex] = 0;
    }

    @SuppressWarnings("unchecked")
    private void resize() {
        if ( buckets == null ) {
            throw new IllegalStateException("Someone must be doing something ugly with reflection -- I have no buckets.");
        }
        final int oldCapacity = capacity;
        final int oldSize = size;
        final T[] oldBuckets = buckets;
        final byte[] oldStatus = status;

        capacity = SetSizeUtils.getLegalSizeAbove(capacity);
        size = 0;
        buckets = (T[])new Object[capacity];
        status = new byte[capacity];

        try {
            int idx = 0;
            do {
                final T entry = oldBuckets[idx];
                if ( entry != null ) insert(entry, (t1, t2) -> false );
            }
            while ( (idx = (idx+127)%oldCapacity) != 0 );
        } catch ( final IllegalStateException ise ) {
            capacity = oldCapacity;
            size = oldSize;
            buckets = oldBuckets;
            status = oldStatus;
            // this shouldn't happen except in the case of really bad hashCode implementations
            throw new IllegalStateException("Hopscotching failed at load factor "+1.*size/capacity+", and resizing didn't help.");
        }

        if ( size != oldSize ) {
            // this should never happen, period.
            throw new IllegalStateException("Lost some elements during resizing.");
        }
    }

    private static int computeCapacity( final int size ) {
        // For power-of-2 table sizes substitute these lines
        /*
        if ( size > maxCapacity ) throw new IllegalArgumentException("Table can't be that big.");
        size = (int)(size/LOAD_FACTOR) - 1;
        size |= size >>> 1;
        size |= size >>> 2;
        size |= size >>> 4;
        size |= size >>> 8;
        size |= size >>> 16;
        return (size < 256) ? 256 : (size >= maxCapacity) ? maxCapacity : size + 1;
        */
        if ( size < LOAD_FACTOR*Integer.MAX_VALUE ) {
            final int augmentedSize = (int) (size / LOAD_FACTOR);
            for ( final int legalSize : SetSizeUtils.legalSizes ) {
                if ( legalSize >= augmentedSize ) return legalSize;
            }
        }
        return SetSizeUtils.legalSizes[SetSizeUtils.legalSizes.length-1];
    }

    private abstract class BaseIterator implements Iterator<T> {
        protected int currentIndex;
        protected int prevIndex;
        protected int removeIndex;
        protected int removePrevIndex;

        BaseIterator() { currentIndex = prevIndex = removeIndex = removePrevIndex = NO_ELEMENT_INDEX; }

        @Override
        public void remove() {
            if ( removeIndex == NO_ELEMENT_INDEX ) throw new IllegalStateException("Remove without next.");

            removeAtIndex(removeIndex, removePrevIndex);

            // If we haven't deleted the end of a chain, we'll now have an unseen element under previousElementIndex.
            // So we need to back up and let the user know about it at the next call to the next method.  If we
            // have deleted an end of chain, then the bucket will be empty and no adjustment needs to be made.
            if ( buckets[removeIndex] != null ) {
                currentIndex = removeIndex;
                prevIndex = removePrevIndex;
            }

            // Set state to "invalid to call remove again".
            removeIndex = NO_ELEMENT_INDEX;
        }
    }

    private final class ElementIterator extends BaseIterator {
        private final Object key;

        ElementIterator( final Object key ) {
            this.key = key;
            int bucketIndex = hashToIndex(key);
            if ( !isChainHead(bucketIndex) ) return;
            currentIndex = bucketIndex;
            ensureEquivalence();
        }

        @Override
        public boolean hasNext() { return currentIndex != NO_ELEMENT_INDEX; }

        @Override
        public T next() {
            if ( !hasNext() ) throw new NoSuchElementException("HopscotchCollection.ElementIterator is exhausted.");
            final T result = buckets[currentIndex];
            removeIndex = currentIndex;
            removePrevIndex = prevIndex;
            advanceUntilEquivalent();
            return result;
        }

        @Override
        public void remove() {
            super.remove();
            if ( currentIndex != NO_ELEMENT_INDEX ) {
                ensureEquivalence();
            }
        }

        private void ensureEquivalence() {
            if ( !equivalent(buckets[currentIndex], key) ) {
                advanceUntilEquivalent();
            }
        }

        private void advanceUntilEquivalent() {
            do {
                final int offset = getOffset(currentIndex);
                if ( offset == 0 ) {
                    currentIndex = NO_ELEMENT_INDEX;
                    break;
                }
                prevIndex = currentIndex;
                currentIndex = getIndex(currentIndex, offset);
            }
            while ( !equivalent(buckets[currentIndex], key) );
        }
    }

    private final class CompleteIterator extends BaseIterator {
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

        CompleteIterator() {
            bucketHeadIndex = NO_ELEMENT_INDEX;
            nextBucketHead();
        }

        @Override
        public boolean hasNext() { return bucketHeadIndex != NO_ELEMENT_INDEX; }

        @Override
        public T next() {
            if ( !hasNext() ) throw new NoSuchElementException("Iterator exhausted.");

            removeIndex = currentIndex;
            removePrevIndex = prevIndex;

            final int offset = getOffset(currentIndex);
            // if we're at the end of a chain, advance to the next bucket
            if ( offset == 0 ) nextBucketHead();
            else { // otherwise step to the next item in the chain
                prevIndex = currentIndex;
                currentIndex = getIndex(currentIndex, offset);
            }

            return buckets[removeIndex];
        }

        private void nextBucketHead() {
            while ( ++bucketHeadIndex < buckets.length ) {
                if ( isChainHead(bucketHeadIndex) ) {
                    currentIndex = bucketHeadIndex;
                    prevIndex = NO_ELEMENT_INDEX;
                    return;
                }
            }
            bucketHeadIndex = NO_ELEMENT_INDEX;
        }
    }

    public static final class Serializer<T> extends com.esotericsoftware.kryo.Serializer<HopscotchCollection<T>> {
        @Override
        public void write( final Kryo kryo, final Output output, final HopscotchCollection<T> hopscotchCollection ) {
            hopscotchCollection.serialize(kryo, output);
        }

        @Override
        public HopscotchCollection<T> read( final Kryo kryo, final Input input,
                                            final Class<HopscotchCollection<T>> klass ) {
            return new HopscotchCollection<>(kryo, input);
        }
    }
}
