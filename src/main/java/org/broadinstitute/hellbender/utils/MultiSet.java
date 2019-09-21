package org.broadinstitute.hellbender.utils;

import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;

import java.lang.reflect.Array;
import java.util.*;
import java.util.stream.Collectors;

/**
 * A class to represent a set where each element might be present an arbitrary number of times 1 being the minimum.
 *
 * Created by gauthier on 9/25/15.
 */
public final class MultiSet<T>  implements Collection<T> {

    /**
     * Default separator between distinct element value-count pairs in {@link #toString} operations.
     */
    public static final String DEFAULT_ELEMENT_SEPARATOR = ",";


    /**
     * Default separator between element value and count in {@link #toString} operations.
     */
    public static final String DEFAULT_COUNT_SEPARATOR = ":";

    /**
     * Holds the number of repeats (values) for each element (key).
     * <p>
     *     The implementation code guarateens that every element has a repeat count of at least 1, so there is no "phantom"
     *     elements with repeats of 0 or negative values. As soon as these would occur they are removed from the map.
     * </p>
     *
     */
    private final Object2IntMap<T> valueCounts;

    /**
     * Holds the total number of elements in the multiset considering repeat counts.
     *
     * <p>
     *     This number could be much larger than the max possible integer value since
     *     is relatively cheap memory wise to have multi-sets with just a few element but very large number of repeats.
     *     For this reason we use a {@code long} type field.
     * </p>
     */
    private long size;


    public MultiSet() {
        this.valueCounts = new Object2IntOpenHashMap<>();
    }

    public MultiSet(final int distinctCapacity) {
        this.valueCounts = new Object2IntOpenHashMap<>(distinctCapacity);
    }

    public MultiSet(final Collection<? extends T> start) {
        this.valueCounts = new Object2IntOpenHashMap<>(start.size());
        addAll(start);
    }

    public MultiSet(final MultiSet<? extends T> start) {
        this.valueCounts = new Object2IntOpenHashMap<>(start.valueCounts.size());
        for (final Object2IntMap.Entry<? extends T> entry :start.valueCounts.object2IntEntrySet()) {
            this.valueCounts.put(entry.getKey(), entry.getIntValue());
        }
        size = start.size;
    }

    public MultiSet(final Iterable<? extends T> start) {
        this.valueCounts = new Object2IntOpenHashMap<>();
        addAll(start);
    }

    /**
     * Returns the number of occurrence of an element in the multi-set.
     * <p>
     *     This method will return 0 for an input value that is not an element of the multi-set.
     * </p>
     * <p>
     *     For members it will always return 1 or greater.
     * </p>
     * @param value the target element.
     * @return 0 or greater.
     */
    public int multificity(final Object value) {
        return valueCounts.getOrDefault(value,0);
    }

    /**
     * Returns the total size of this multi-set considering all element multiplicities.
     * <p>
     *     The result size is down-casted version of the more precise {@link #longSize()} into an {@code int} for compatibility with Java standard library collection
     *     library.
     * </p>
     * @return 0 or greater, but could be negative in case of integer overflow.
     */
    @Override
    public int size() {
        return (int) size;
    }

    /**
     * Returns the total (long) size of this multi-set considering all element multiplicities.
     *
     * @return 0 or greater.
     */
    public long longSize() {
        return size;
    }

    /**
     * Returns the number of distinct elements in the multi-set.
     * @return 0 or greater.
     */
    public int distinctSize() {
        return valueCounts.size();
    }

    public boolean isEmpty(){
        return valueCounts.isEmpty();
    }

    @Override
    public boolean contains(final Object o) {
        return valueCounts.containsKey(o);
    }

    /**
     * Returns a "live" set with the distinct values present in this multi-set.
     * <p>
     *     Changes in this multi-set would be visible in the returned set.
     * </p>
     * @return never {@code null}
     */
    public Set<T> distinct() {
        return Collections.unmodifiableSet(valueCounts.keySet());
    }

    @Override
    public Iterator<T> iterator(){

        return new Iterator<T>() {

            private final Iterator<Object2IntMap.Entry<T>> countsIterator = valueCounts.object2IntEntrySet().iterator();
            private int currentRemaining = 0;
            private T current = null;

            @Override
            public boolean hasNext() {
                return currentRemaining > 0 || countsIterator.hasNext();
            }

            @Override
            public T next() {
                if (currentRemaining > 0) {
                    currentRemaining--;
                    return current;
                } else if (countsIterator.hasNext()) {
                    final Object2IntMap.Entry<T> nextEntry = countsIterator.next();
                    currentRemaining = nextEntry.getIntValue() - 1;
                    return current = nextEntry.getKey();
                } else {
                    throw new NoSuchElementException("trying to iterate beyond the end of this collection");
                }
            }
        };
    }

    @Override
    public Object[] toArray() {
        return new Object[0];
    }

    @Override
    @SuppressWarnings("unchecked")
    public <T> T[] toArray(T[] a) {
        final T[] result;
        if (a.length >= size) {
            result = a;
        } else if (size > Integer.MAX_VALUE) {
            throw new UnsupportedOperationException("cannot create arrays larger than the max integer value");
        } else {
            result = (T[]) Array.newInstance( a.getClass().getComponentType(), (int) size);
        }
        int i = 0;
        for (final Object2IntMap.Entry<?> entry : valueCounts.object2IntEntrySet()) {
            final T value = (T) entry.getKey();
            for (int j = entry.getIntValue(); j > 0; j--) {
                result[i++] = value;
            }
        }
        return result;
    }

    @Override
    public String toString()  {
        return toString(DEFAULT_ELEMENT_SEPARATOR, DEFAULT_COUNT_SEPARATOR);
    }

    public String toString(final String sep, final String countSep)  {
        return valueCounts.object2IntEntrySet().stream()
                .map(entry -> "" + entry.getKey() + countSep + entry.getIntValue())
                .collect(Collectors.joining(sep));
    }

    public String toString(final Comparator<? super T> comparator)  {
        return toString(comparator, DEFAULT_ELEMENT_SEPARATOR, DEFAULT_COUNT_SEPARATOR);
    }

    public String toString(final Comparator<? super T> comparator, final String sep, final String countSep) {
        return valueCounts.object2IntEntrySet().stream()
                .sorted(Comparator.comparing(Object2IntMap.Entry::getKey,comparator))
                .map(entry -> "" + entry.getKey() + countSep + entry.getIntValue())
                .collect(Collectors.joining(sep));
    }

    public boolean add(final T val) {
        return add(val, 1);
    }

    @Override
    @SuppressWarnings("unchecked")
    public boolean remove(final Object o) {
        if (valueCounts.containsKey(o)) {
            final int newCount = valueCounts.get(o) - 1;
            if (newCount <= 0) {
                valueCounts.remove(o);
            } else {
                valueCounts.put((T) o, newCount);
            }
        }
        return false;
    }

    public boolean containsAll(final MultiSet<?> other) {
        if (other.size > size) {
            return false;
        } else {
            for (final Object2IntMap.Entry<?> otherEntry : other.valueCounts.object2IntEntrySet()) {
                if (!valueCounts.containsKey(otherEntry.getKey())) {
                    return false;
                } else {
                    final int otherMultiplicity = otherEntry.getIntValue();
                    final int thisMultiplicity = valueCounts.getInt(otherEntry.getKey());
                    if (otherMultiplicity > thisMultiplicity) {
                        return false;
                    }
                }
            }
            return true;
        }
    }

    @Override
    public boolean containsAll(final Collection<?> other) {
        return containsAll(new MultiSet<>(other));
    }

    @Override
    public boolean removeAll(final Collection<?> c) {
        boolean result = false;
        for (final Object other : c) {
            result |= remove(other);
        }
        return result;
    }

    public boolean removeAll(final MultiSet<?> other) {
        boolean result = false;
        for (final Object2IntMap.Entry<?> entry : other.valueCounts.object2IntEntrySet()) {
            final Object otherKey = entry.getKey();
            if (valueCounts.containsKey(otherKey)) {
                final int oldCount = valueCounts.getInt(otherKey);
                final int newCount = oldCount - entry.getIntValue();
                if (newCount <= 0) {
                    size -= oldCount;
                    valueCounts.remove(otherKey);
                } else {
                    @SuppressWarnings("unchecked")
                    final T otherKeyCasted = (T) otherKey;
                    valueCounts.put(otherKeyCasted, newCount);
                    size -= entry.getIntValue();
                }
                result = true;
            }
        }
        return result;
    }

    public boolean retainAll(final MultiSet<?> other) {
        boolean result = false;
        final Iterator<Object2IntMap.Entry<T>> it = valueCounts.object2IntEntrySet().iterator();
        while (it.hasNext()) {
            final Object2IntMap.Entry<T> entry = it.next();
            final T key = entry.getKey();
            final int otherMultiplicity = other.multificity(key);
            if (otherMultiplicity == 0) {
                it.remove();
                size -= entry.getIntValue();
                result = true;
            } else if (otherMultiplicity < entry.getIntValue()) {
                size -= entry.getIntValue() - otherMultiplicity;
                entry.setValue(otherMultiplicity);
                result = true;
            }
        }
        return result;
    }

    @Override
    public boolean retainAll(final Collection<?> c) {
        return retainAll(new MultiSet<>(c));
    }

    @Override
    public void clear() {
       size = 0;
       valueCounts.clear();
    }

    public boolean add(final T val, final int count){
        if (count == 0) { // adding 0 is same as doing nothing.
            return false;
        } else if (count < 0) {
            throw new IllegalArgumentException("the count provided cannot be negative");
        } else if (valueCounts.containsKey(val)) {
            valueCounts.put(val, valueCounts.getInt(val) + count);
        } else {
            valueCounts.put(val, count);
        }
        size += count;
        return true;
    }

    /**
     * Add all the elements in the input multi-set into this multi-set.
     * <p>
     *     Repeated elements will be handled appropriately by increasing the count
     *     of those element in the multi-set.
     * </p>
     * @param other the (repeated) elements to add.
     */
    public boolean addAll(final MultiSet<? extends T> other) {
        if (other.size == 0) {
            return false;
        } else {
            for (final Object2IntMap.Entry<? extends T> otherEntry : other.valueCounts.object2IntEntrySet()) {
                this.add(otherEntry.getKey(), otherEntry.getIntValue());
            }
            return true;
        }
    }

    /**
     * Add all the elements in the input collection into this multi-set.
     * <p>
     *     Repeated elements will be handled appropriately by increasing the count
     *     of those element in the multi-set.
     * </p>
     * @param other the elements to add.
     */
    public boolean addAll(final Iterable<? extends T> other) {
        boolean result = false;
        for (final T otherElement : other) {
            add(otherElement);
            result = true;
        }
        return result;
    }

   @Override
   public boolean addAll(final Collection<? extends T> other) {
      if (other.isEmpty()) {
          return false;
      } else {
          return addAll((Iterable<? extends T>) other);
      }
   }

}
