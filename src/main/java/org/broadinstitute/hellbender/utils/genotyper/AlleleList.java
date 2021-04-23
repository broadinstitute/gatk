package org.broadinstitute.hellbender.utils.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.AbstractList;
import java.util.List;

/**
 * Minimal interface for random access to a collection of Alleles.
 */
//Note: Names in this interface are unusual because of name clash in a subclass.
// For example the name of AlleleList.alleleCount() cannot be simply size(), as would be usual,
// because {@link ReadLikelihoods} implements AlleleList and SampleList and then size() would be ambiguous.
public interface AlleleList<A extends Allele> {

    static <A extends Allele> AlleleList<A> newList(final List<A> alleles) {
        return new IndexedAlleleList<A>(alleles);
    }

    /**
     * Returns the number of alleles in this AlleleList.
     */
    int numberOfAlleles();

    /**
     * Returns the index of the given Allele in this AlleleList.
     * Returns a negative number if the given allele is not present in this AlleleList.
     * @throws IllegalArgumentException if allele is null.
     */
    int indexOfAllele(final Allele allele);

    /**
     * Returns the allele at the given index in this AlleleList.
     * @throws IllegalArgumentException if index is negative or equal
     * to or higher than the number of elements in this AlleleList {@link AlleleList#numberOfAlleles()}).
     */
    A getAllele(final int index);

    /**
     * Returns <code>true</code> if this AlleleList contains the specified allele
     * and <code>false</code> otherwise.
     */
    default boolean containsAllele(final Allele allele) {
        return indexOfAllele(allele) >= 0;
    }

    AlleleList<Allele> EMPTY_LIST = new AlleleList<Allele>() {
        @Override
        public int numberOfAlleles() {
            return 0;
        }

        @Override
        public int indexOfAllele(final Allele allele) {
            Utils.nonNull(allele);
            return -1;
        }

        @Override
        public Allele getAllele(final int index) {
            throw new IllegalArgumentException("allele index is out of range");  //we know this without checking because it's an empty list
        }
    };

    /**
     * Returns an unmodifiable empty allele-list.
     * @param <A> the allele class.
     * @return never {@code null}.
     */
    @SuppressWarnings("unchecked")
    static <A extends Allele> AlleleList<A> emptyAlleleList() {
        return (AlleleList<A>) EMPTY_LIST;
    }

    /**
     * Checks whether two allele lists are in fact the same.
     * @param first one list to compare.
     * @param second another list to compare.
     *
     * @throws IllegalArgumentException if if either list is {@code null}.
     *
     * @return {@code true} iff both list are equal.
     */
    static <A extends Allele> boolean equals(final AlleleList<A> first, final AlleleList<A> second) {
        if (first == null || second == null) {
            throw new IllegalArgumentException("no null list allowed");
        }
        final int alleleCount = first.numberOfAlleles();
        if (alleleCount != second.numberOfAlleles()) {
            return false;
        }

        for (int i = 0; i < alleleCount; i++) {
            final A firstSample = first.getAllele(i);
            Utils.nonNull(firstSample, "no null samples allowed in sample-lists: first list at " + i);
            final A secondSample = second.getAllele(i);
            Utils.nonNull(secondSample,"no null samples allowed in sample-list: second list at " + i);
            if (!firstSample.equals(secondSample)) {
                return false;
            }
        }

        return true;
    }

    /**
     * Resolves the index of the reference allele in an allele-list.
     *
     * <p>
     *     If there is no reference allele, it returns -1. If there is more than one reference allele,
     *     it returns the first occurrence (lowest index).
     * </p>
     *
     *
     * @throws IllegalArgumentException if {@code list} is {@code null}.
     *
     * @return -1 if there is no reference allele, or a values in [0,{@code list.alleleCount()}).
     */
    default int indexOfReference() {
        final int alleleCount = this.numberOfAlleles();
        for (int i = 0; i < alleleCount; i++) {
            if (this.getAllele(i).isReference()) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Returns a {@link List} unmodifiable view of this allele-list
     *
     * @return never {@code null}.
     */
    default List<A> asListOfAlleles() {
        return new AbstractList<A>() {

            @Override
            public A get(final int index) {
                return AlleleList.this.getAllele(index);
            }

            @Override
            public int size() {
                return AlleleList.this.numberOfAlleles();
            }
        };
    }

    /**
     * Returns a permutation between two allele lists.
     * @param target the target allele list.
     *
     * @throws IllegalArgumentException if {@code target} is {@code null}, or
     * elements in {@code target} is not contained in {@code this}
     *
     * @return never {@code null}
     */
    default AlleleListPermutation<A> permutation(final AlleleList<A> target) {
        if (equals(this, target)) {
            return new NonPermutation<>(this);
        } else {
            return new ActualPermutation<>(this, target);
        }
    }

    /**
     * This is the identity permutation.
     */
    final class NonPermutation<A extends Allele> implements AlleleListPermutation<A> {

        private final AlleleList<A> list;

        NonPermutation(final AlleleList<A> original) {
            list = original;
        }

        @Override
        public boolean isPartial() {
            return false;
        }

        @Override
        public boolean isNonPermuted() {
            return true;
        }

        @Override
        public int toIndex(final int fromIndex) {
            return fromIndex;
        }

        @Override
        public int fromIndex(final int toIndex) {
            return toIndex;
        }

        @Override
        public boolean isKept(final int fromIndex) { return true; }

        @Override
        public int fromSize() {
            return list.numberOfAlleles();
        }

        @Override
        public int toSize() {
            return list.numberOfAlleles();
        }

        @Override
        public List<A> fromList() {
            return list.asListOfAlleles();
        }

        @Override
        public List<A> toList() {
            return list.asListOfAlleles();
        }

        @Override
        public int numberOfAlleles() {
            return list.numberOfAlleles();
        }

        @Override
        public int indexOfAllele(final Allele allele) {
            return list.indexOfAllele(allele);
        }

        @Override
        public A getAllele(final int index) {
            return list.getAllele(index);
        }
    }

    final class ActualPermutation<A extends Allele> implements AlleleListPermutation<A> {

        private final AlleleList<A> from;

        private final AlleleList<A> to;

        private final int[] fromIndex;

        private final boolean[] keptFromIndices;
        
        private final boolean nonPermuted;

        private final boolean isPartial;

        private ActualPermutation(final AlleleList<A> original, final AlleleList<A> target) {
            this.from = original;
            this.to = target;
            keptFromIndices = new boolean[original.numberOfAlleles()];
            final int toSize = target.numberOfAlleles();
            final int fromSize = original.numberOfAlleles();
            if (fromSize < toSize) {
                throw new IllegalArgumentException("target allele list is not a permutation of the original allele list");
            }

            fromIndex = new int[toSize];
            boolean nonPermuted = fromSize == toSize;
            this.isPartial = !nonPermuted;
            for (int i = 0; i < toSize; i++) {
                final int originalIndex = original.indexOfAllele(target.getAllele(i));
                if (originalIndex < 0) {
                    throw new IllegalArgumentException("target allele list is not a permutation of the original allele list");
                }
                keptFromIndices[originalIndex] = true;
                fromIndex[i] = originalIndex;
                nonPermuted &= originalIndex == i;
            }

            this.nonPermuted = nonPermuted;
        }

        @Override
        public boolean isPartial() {
            return isPartial;
        }

        @Override
        public boolean isNonPermuted() {
            return nonPermuted;
        }

        @Override
        public int toIndex(final int fromIndex) {
            return to.indexOfAllele(from.getAllele(fromIndex));
        }

        @Override
        public int fromIndex(final int toIndex) {
            return fromIndex[toIndex];
        }

        @Override
        public boolean isKept(final int fromIndex) {
            return keptFromIndices[fromIndex];
        }

        @Override
        public int fromSize() {
            return from.numberOfAlleles();
        }

        @Override
        public int toSize() {
            return to.numberOfAlleles();
        }

        @Override
        public List<A> fromList() {
            return from.asListOfAlleles();
        }

        @Override
        public List<A> toList() {
            return to.asListOfAlleles();
        }

        @Override
        public int numberOfAlleles() {
            return to.numberOfAlleles();
        }

        @Override
        public int indexOfAllele(final Allele allele) {
            return to.indexOfAllele(allele);
        }

        @Override
        public A getAllele(final int index) {
            return to.getAllele(index);
        }
    }
}
