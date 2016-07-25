package org.broadinstitute.hellbender.utils.collections;

import java.util.Arrays;
import java.util.NoSuchElementException;

/**
 * Implementation of a small set of integer primitive values.
 * <p>
 * It does not provide a O(1) look-up of its elements. Elements are kept in a sorted array so look up is implemented
 * using a binary search O(log n). Therefore it might not be optimal for problems that require large integer sets.
 * </p>
 * <p>
 * Also note that addition can be costly for large sets unless done in order: O(n).
 * </p>
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class CountSet {

    /**
     * The size of the set.
     */
    private int size;

    /**
     * Holds the element of the set within the subrange <code>[0 .. size - 1]</code> in ascending order.
     */
    private int[] elements;

    /**
     * Creates a new set indicating the expected maximum number of elements it will contain.
     * @param initialCapacity the desired initial capacity of the set.
     * @throws IllegalArgumentException if <code>initialCapacity</code> is negative.
     */
    public CountSet(final int initialCapacity) {
        if (initialCapacity < 0) {
            throw new IllegalArgumentException("initialCapacity must be non-negative but was " + initialCapacity);
        }
        elements = new int[initialCapacity];
        size = 0;
    }

    /**
     * Returns the smallest integer value in the set.
     *
     * @throws NoSuchElementException if the set is empty (thus there is no minimum).
     * @return the smallest integer value in the set.
     */
    public int min() {
        if (size == 0) {
            throw new NoSuchElementException("cannot have a min from an empty set");
        }
        return elements[0];
    }

    /**
     * Returns the largest integer value in the set.
     *
     * @throws NoSuchElementException if the set is empty (thus there is no maximum).
     * @return the largest integer value in the set.
     */
    public int max() {
        if (size == 0) {
            throw new NoSuchElementException("cannot have a max from an empty set");
        }
        return elements[size - 1];
    }

    /**
     * Add an integer value to the set.
     * @param value to add to the set.
     * @return <code>true</code> if the set changed as a result of this invocation, <code>false</code> otherwise.
     */
    public boolean add(final int value) {
        final int pos = Arrays.binarySearch(elements, 0, size, value);
        if (pos >= 0) return false;
        final int insertPos = - pos - 1;
        ensureCapacity(size + 1);
        System.arraycopy(elements, insertPos, elements, insertPos + 1, size - insertPos);
        elements[insertPos] = value;
        size++;
        return true;
    }

    /**
     * Make sure that this int-set has capacity to handle a number of elements.
     * <p/>
     * If the set has already that or greater capacity nothing would be changed.
     *
     * @param capacity the requested capacity.
     */
    private void ensureCapacity(final int capacity) {
        if (elements.length >= capacity) return;
        final int newLength = Math.max(elements.length << 1, capacity);
        elements = Arrays.copyOf(elements, newLength);
    }

    public int size() {
        return size;
    }

    public boolean isEmpty() {
        return size() == 0;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder(2 + size() * 10);
        sb.append('{');
        for (int i = 0; i < size; i++) {
            sb.append(elements[i]).append(',');
        }
        sb.replace(sb.length()-1,sb.length(),"}");
        return sb.toString();
    }
}
