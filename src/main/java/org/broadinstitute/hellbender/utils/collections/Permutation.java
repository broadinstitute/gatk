package org.broadinstitute.hellbender.utils.collections;

import java.util.List;

/**
 * Represent a permutation of a ordered set or list of elements.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public interface Permutation<E> {

    /**
     * Checks whether this permutation is a partial one of the original list.
     *
     * <p>
     *     A partial permutation is one in that not all original elements take part in.
     * </p>
     *
     * @return {@code true} iff this is a partial permutation.
     */
    public boolean isPartial();

    /**
     * Checks whether this is a trivial permutation where the resulting element list is the same as original.
     *
     * @return {@code true} iff the resulting element list is the same as the original.
     */
    public boolean isNonPermuted();

    /**
     * Given an index on the original list, returns the position of tha element in the resulting list.
     *
     * @param fromIndex the query original element index.
     *
     * @throws IllegalArgumentException if {@code fromIndex} is not a valid index within the original list.
     *
     * @return -1 if that element is not part of the result (partial) permutation, otherwise some number between
     *   0 and {@link #toSize()} - 1.
     */
    public int toIndex(final int fromIndex);

    /**
     * Given an index on the resulting list, it gives you the index of that element on the original list.
     * @param toIndex the query resulting list index.
     *
     * @throws IllegalArgumentException if {@code toIndex} is not a valid index, i.e. in [0,{@link #toSize()}-1).
     *
     * @return a value between 0 and {@link #fromSize()} - 1.
     */
    public int fromIndex(final int toIndex);

    /**
     * Given an index of the original list, return whether this index is found at any position of the permuted list.
     * This is trivial if the permutation is not partial.
     *
     * @param fromIndex
     * @return
     */
    public boolean isKept(final int fromIndex);

    /**
     * Length of the original element list.
     *
     * @return 0 or greater.
     */
    public int fromSize();

    /**
     * Length of the resulting element list.
     *
     * @return 0 or greater.
     */
    public int toSize();

    /**
     * Returns an unmodifiable view to the original element list.
     * @return never {@code null}.
     */
    public List<E> fromList();

    /**
     * Returns an unmodifiable view to the original element list.
     *
     * @return never {@code null}.
     */
    public List<E> toList();
}
