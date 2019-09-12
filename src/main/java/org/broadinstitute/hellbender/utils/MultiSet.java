package org.broadinstitute.hellbender.utils;

import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;

import java.util.*;
import java.util.stream.Collectors;

/**
 * A class to represent a set where each element might be present a arbitrary number of times (at least 1).
 *
 * Created by gauthier on 9/25/15.
 */
public final class MultiSet<T>  implements Iterable<T> {

    private Object2IntMap<T> valueCounts = new Object2IntOpenHashMap<>();

    private long size;

    /**
     * Returns the number of ocurrence of an element in the multi-set.
     * <p>
     *     This method will return 0 for an input value that is not an element of the multi-set.
     * </p>
     * <p>
     *     For members it will always return 1 or greater.
     * </p>
     * @param value the target element.
     * @return 0 or greater.
     */
    public int multificity(final T value) {
        return valueCounts.getOrDefault(value,0);
    }

    /**
     * Returns the total size of this multi-set considering all element multiplicities.
     * <p>
     *     The result size is down-casted into an {@code int} for compatibility with Java standard library collection
     *     library.
     * </p>
     * @return 0 or greater, but could be negative in case of integer overflow.
     */
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
    public String toString()  {
        return valueCounts.object2IntEntrySet().stream()
                .map(entry -> "" + entry.getKey() + "," + entry.getIntValue())
                .collect(Collectors.joining(","));
    }


    public void add(final T val) {
        add(val, 1);
    }

    public void add(final T val, final int count){
        if (count == 0) { // adding 0 is same as doing nothing.
            return;
        } else if (count < 0) {
            throw new IllegalArgumentException("the count provided cannot be negative");
        } else if (valueCounts.containsKey(val)) {
            valueCounts.put(val, valueCounts.getInt(val) + count);
        } else {
            valueCounts.put(val, count);
        }
        size += count;
    }

    /**
     * Add all the elements in the input multi-set into this multi-set.
     * <p>
     *     Repeated elements will be handled appropriately by increasing the count
     *     of those element in the multi-set.
     * </p>
     * @param other the (repeated) elements to add.
     * @param <U> the type-param for the input element type/subtype.
     */
    public <U extends T> void addAll(final MultiSet<U> other) {
        for (final Object2IntMap.Entry<U> otherEntry : other.valueCounts.object2IntEntrySet()) {
            this.add(otherEntry.getKey(), otherEntry.getIntValue());
        }
    }

    /**
     * Add all the elements in the input collection into this multi-set.
     * <p>
     *     Repeated elements will be handled appropriately by increasing the count
     *     of those element in the multi-set.
     * </p>
     * @param other the elements to add.
     * @param <U> the type-param for the input element type/subtype.
     */
    public <U extends T> void addAll(final Iterable<U> other) {
        for (final U otherElement : other) {
            add(otherElement);
        }
    }
}
