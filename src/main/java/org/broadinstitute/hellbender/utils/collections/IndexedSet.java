package org.broadinstitute.hellbender.utils.collections;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
* Set where each element can be reference by a unique integer index that runs from
*     0 to the size of the set - 1.
*
* @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
*/
public final class IndexedSet<E> extends AbstractSet<E> {

    /**
     * Elements stored in an array-list by their index.
     */
    private final List<E> elements;

    /**
     * A unmodifiable view to the element list. Initially {@code null} it is thread-unsafe lazy instantiated
     * when requested first time through {@link #asList}. Therefore typically it is shared by invoking code but
     * there could be some extra copies (rare though) in multi-thread runs.
     */
    private List<E> unmodifiableElementsListView;

    /**
     * Quick element to index lookup map.
     * <p>
     *  Uses a primitive int value map for efficiency sake.
     * </p>
     */
    private final Map<E, Integer> indexByElement;

    /**
     * Creates an empty indexed set indicating the expected number of elements.
     *
     * @param initialCapacity the initial number of elements.
     */
    public IndexedSet(final int initialCapacity) {
        elements = new ArrayList<>(initialCapacity);
        indexByElement = new LinkedHashMap<>(initialCapacity);
    }

    /**
     * Creates a new sample list from a existing collection of elements.
     *
     * <p>
     *     Elements will be indexed as they appear in the input array. Repeats will be ignored.
     * </p>
     *
     * @param values the original sample list.
     *
     * @throws IllegalArgumentException
     * if {@code values} array is {@code null} itself, or it contains {@code null}.
     */
    public IndexedSet(final Collection<E> values) {
        Utils.nonNull(values, "input values cannot be null");
        Utils.containsNoNull(values, "values can't contain nulls");

        final int initialCapacity = values.size();
        elements = new ArrayList<>(initialCapacity);
        indexByElement = new LinkedHashMap<>(initialCapacity);
        int nextIndex = 0;
        for (final E value : values) {
            if (indexByElement.containsKey(value)) {
                continue;
            }
            indexByElement.put(value, nextIndex++);
            elements.add(value);
        }
    }

    /**
     * Creates a new sample list from a existing array of elements.
     *
     * <p>
     *     Elements will be indexed as they appear in the collection. Repeats will be ignored.
     * </p>
     *
     * @param values the original sample list.
     *
     * @throws IllegalArgumentException
     * if {@code values} collection is {@code null} itself, or it contains {@code null}.
     */
    @SafeVarargs
    @SuppressWarnings("varargs")
    public IndexedSet(final E ... values) {
        Utils.nonNull(values, "input values cannot be null");

        final int initialCapacity = values.length;
        elements = new ArrayList<>(initialCapacity);
        indexByElement = new LinkedHashMap<>(initialCapacity);
        int nextIndex = 0;
        for (final E value : values) {
            Utils.nonNull(value, "null element not allowed: index == " + nextIndex);

            if (indexByElement.containsKey(value)) {
                continue;
            }
            indexByElement.put(value, nextIndex++);
            elements.add(value);
        }
    }

    /**
     * Returns a list view of the elements in the set.
     *
     * <p>
     *     Elements are sorted by their index within the set.
     * </p>
     *
     * <p>
     *     This view changes as the indexed set changes but it cannot be used to update its contents.
     *     In such case a {@link UnsupportedOperationException} exception will be thrown if the calling
     *     code tries to tho just that.
     * </p>
     *
     * @return never {@code null}.
     */
    public List<E> asList() {
        if (unmodifiableElementsListView == null) {
            unmodifiableElementsListView = Collections.unmodifiableList(elements);
        }
        return unmodifiableElementsListView;
    }

    @Override
    public Iterator<E> iterator() {
        return asList().iterator();
    }

    /**
     * Returns number of elements in the set.
     * @return 0 or greater
     */
    @Override
    public int size() {
        return elements.size();
    }

    /**
     *
     * @param o
     * @return {@code true} iff {@code o} is in
     */
    @Override
    public boolean contains(final Object o) {
        return o != null && indexByElement.containsKey(o);
    }

    /**
     * Adds a new element to the set.
     *
     * <p>
     *     If the element was already in th set nothing will happen and the method will return {@code false}. However,
     *     if the element is new to this set, it will assigned the next index available (equal to the size before addition).
     *     The method will return {@code true} in this case.
     * </p>
     *
     * @param o the object to add.
     *
     * @throw IllegalArgumentException if {@code o} is {@code null}.
     *
     * @return {@code true} iff the set was modified by this operation.
     */
    @Override
    public boolean add(final E o) {
        Utils.nonNull(o, "the input argument cannot be null");
        if (contains(o)) {
            return false;
        }
        final int nextIndex = size();
        elements.add(o);
        indexByElement.put(o, nextIndex);
        return true;
    }

    /**
     * Removes an element from the set.
     *
     * <p>
     *     If the element was not present in the set, nothing happens and the method return false. However,
     *     if the element is new to this set, it will be assigned the next index available (equal to the size
     *     before addition).
     *     The method will return {@code true} in this case.
     * </p>
     *
     * @param o the object to remove.
     *
     * @throw IllegalArgumentException if {@code o} is {@code null}.
     *
     * @return {@code true} iff the set was modified by this operation.
     */
    @Override
    public boolean remove(final Object o) {
        final int index = indexOf(o);
        if (index == -1) {
            return false;
        }
        elements.remove(index);
        indexByElement.remove(o);
        final ListIterator<E> it = elements.listIterator(index);
        int nextIndex = index;
        while (it.hasNext()) {
            indexByElement.put(it.next(), nextIndex++);
        }
        return true;
    }

    /**
     * Removes all elements in the set.
     */
    @Override
    public void clear() {
        elements.clear();
        indexByElement.clear();
    }

    /**
     * Compares this with another indexed set.
     * @param o the other object to compare to.
     * @return {@code false} unless {@code o} is a indexed-set that contains the same elements in the same order.
     */
    @Override
    public boolean equals(final Object o) {
        if (o == this) {
            return true;
        }
        if (!(o instanceof IndexedSet<?>)) {
            return false;
        }

        final IndexedSet<?> other = (IndexedSet<?>)o;

        return equals(other);
    }

    /**
     * Compare to another indexed set.
     *
     * @param other the target indexed set.
     *
     * @throws IllegalArgumentException if {@code other} is {@code null}.
     *
     * @return {@code true} iff {@code other} is not {@code null}, and contains exactly the same elements
     * (as compared using {@link Object#equals} a this set with matching indices.
     */
    public boolean equals(final IndexedSet<?> other) {
        Utils.nonNull(other, "other cannot be null");
        final List<?> otherElements = other.elements;

        final int elementCount = elements.size();
        if (otherElements.size() != elementCount) {
            return false;
        }
        for (int i = 0; i < elementCount; i++) {
            if (!elements.get(i).equals(otherElements.get(i))) {
                return false;
            }
        }
        return true;
    }

    @Override
    public int hashCode() {
        int result = 1;

        for (final E element : elements) {
            result = 31 * result + (element == null ? 0 : element.hashCode());
        }
        return result;
    }

    /**
     * Returns the element given its index within the set.
     * @param index the target element's index.
     *
     * @throws IllegalArgumentException if {@code index} is not valid; in [0,{@link #size()}).
     *
     * @return never {@code null}; as null is not a valid element.
     */
    public E get(final int index) {
        Utils.validIndex(index, size());
        return elements.get(index);
    }

    /**
     * Returns the index of an object.
     * @param o the object of interest.
     *
     * @throws IllegalArgumentException if {@code o} is {@code null}.
     *
     * @return {@code -1} if such an object is not an element of this set, otherwise is index in the set thus a
     * values within [0,{@link #size()}).
     */
    public int indexOf(final Object o) {
        Utils.nonNull(o, "the query object cannot be null");
        return indexByElement.getOrDefault(o, -1);
    }

}
