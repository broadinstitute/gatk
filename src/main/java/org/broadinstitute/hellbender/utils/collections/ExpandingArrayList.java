package org.broadinstitute.hellbender.utils.collections;

import java.util.ArrayList;

/**
 * An expanding list is a list that grows when you try to insert an element at a position that does not exist yet.
 * @param <E> type of the elements.
 */
public final class ExpandingArrayList<E> extends ArrayList<E> {

    private static final long serialVersionUID = 1;

    public ExpandingArrayList() { super(); }
    public ExpandingArrayList(int initialCapacity) { super(initialCapacity); }

    /**
     * Returns the element at the specified position in this list.  If index > size,
     * returns null.  Otherwise tries to access the array
     * @throws IndexOutOfBoundsException in index < 0
     */
    public E get(int index) throws IndexOutOfBoundsException {
        if ( index < size() )
            return super.get(index);
        else
            return null;
    }

    private void maybeExpand(int index, E value) {
        if ( index >= size() ) {
            ensureCapacity(index+1); // make sure we have space to hold at least index + 1 elements
            // We need to add null items until we can safely set index to element
            for ( int i = size(); i < index; i++ ){
                add(null);
            }
            //here, size() == index
            add(value);
        }
    }

    /**
     * Sets the element at a given index. Grows the array if needed.
     */
    public E set(int index, E element) {
        maybeExpand(index, null);
        return super.set(index, element);
    }
}