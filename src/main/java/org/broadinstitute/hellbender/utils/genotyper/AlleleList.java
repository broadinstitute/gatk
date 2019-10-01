package org.broadinstitute.hellbender.utils.genotyper;

import htsjdk.variant.variantcontext.Allele;
import it.unimi.dsi.fastutil.ints.IntList;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
 * Minimal interface for random access to a collection of Alleles.
 */
//Note: Names in this interface are unusual because of name clash in a subclass.
// For example the name of AlleleList.alleleCount() cannot be simply size(), as would be usual,
// because {@link ReadLikelohoods} implements AlleleList and SampleList and then size() would be ambiguous.
public interface AlleleList<A extends Allele> {

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
     * Returns the index of the given Allele in the AlleleList ignoring its reference status or returned
     * the index of NON-REF when it applies.
     *
     * <p>
     *     When {@code failOverToNonRef} is {@code true}, this method will conside to return its index iff
     *     the input allele is not a reference allele or {@code ignoreRefStatus} is {@code true}. So just to make
     *     this clear, if {@code ignoreRefStatus} is {@code false} this method won't ever fail over to NON-REF
     *     with reference input alleles.
     * </p>
     *
     * @param allele
     * @param ignoreRefStatus
     * @param failOverToNonRef
     * @return
     */
    default int indexOfAllele(final Allele allele, final boolean ignoreRefStatus, final boolean failOverToNonRef) {
        int result = indexOfAllele(allele);
        if (result == -1 && ignoreRefStatus && !allele.isSymbolic() && !allele.isNoCall()) {
            final Allele secondAttemptAllele = Allele.create(allele.getBases(), !allele.isReference());
            result = indexOfAllele(secondAttemptAllele);
        }
        if (result == -1 && failOverToNonRef && (ignoreRefStatus || !allele.isReference())) {
            result = nonRefAlleleIndex();
        }
        return result;
    }

    default int indexOfAllele(final Allele allele, final Allele ref, final boolean ignoreRefStatus, final boolean failOverToNonRef) {
        if (ref == null || allele.isSymbolic() || allele.isNoCall()) {
            return indexOfAllele(allele, ignoreRefStatus, failOverToNonRef);
        } else {
            int thisRefAlleleIndex = indexOfUniqueReference();
            final Allele searchAllele;
            if (thisRefAlleleIndex == -1) {
                searchAllele = allele;
            } else {
                final Allele thisRefAllele = getAllele(thisRefAlleleIndex);
                final byte[] thisRefBases = thisRefAllele.getBases();
                final byte[] thatRefBases = ref.getBases();
                if (Nucleotide.same(thisRefBases, thatRefBases)) {
                    searchAllele = allele;
                } else if (Nucleotide.startsWith(thisRefBases, thatRefBases)) {
                    final byte[] extensionBases = Arrays.copyOfRange(thisRefBases, thatRefBases.length, thisRefBases.length);
                    searchAllele = Allele.extend(allele, extensionBases);
                } else if (Nucleotide.startsWith(thatRefBases, thisRefBases)) {
                    final int extensionSize = thisRefBases.length - thatRefBases.length;
                    final int to = allele.getBases().length + extensionSize;
                    searchAllele = to > 0 ? Allele.create(Arrays.copyOfRange(allele.getBases(), 0, to), allele.isReference()) : null;
                } else {
                    throw new IllegalArgumentException("cannot match reference");
                }
            }
            return searchAllele == null ? -1 : indexOfAllele(searchAllele, ignoreRefStatus, failOverToNonRef);
        }
    }
    /**
     * Returns an list of ints indicate the location of an allele in the ith position of this list in
     * the input allele list. It will contain -1, in case such allele is missing.
     * <p>
     *     If {@code matchReference} is set to true, it try to match reference allele in this and the input list
     *     and apply allele base extension or contraction before matching alleles.
     * </p>
     * <p>
     *     If {@code useNonRefFailOver} when an allele in this list cannot be matched with an allele in the input list
     *     the NON_REF allele (if present) will be used instead.
     * </p>
     * @param originalList where the returned indexes point to.
     * @return never {@code null}, a list with exactly one position per allele in this list.
     */
    default IntList alleleIndexMap(final AlleleList<?> originalList, final boolean matchReferences,
                                   final boolean useNonRefFailOver) {
        Utils.nonNull(originalList);
        final int basesLengthDifference;
        final byte[] extensionBases;
        if (matchReferences) {
            final int thisRefIdx = indexOfUniqueReference();
            final int thatRefIdx = originalList.indexOfUniqueReference();
            if (thisRefIdx == -1 || thatRefIdx == -1) {
                basesLengthDifference = 0;
                extensionBases = null;
            } else {
                final Allele thisRefAllele = getAllele(thisRefIdx);
                final Allele thatRefAllele = originalList.getAllele(thatRefIdx);
                final byte[] thisRefBases = thisRefAllele.getBases();
                final byte[] thatRefBases = thatRefAllele.getBases();
                if (Nucleotide.same(thisRefBases, thatRefBases)) {
                    basesLengthDifference = 0;
                    extensionBases = null;
                } else if (Nucleotide.startsWith(thatRefBases, thisRefBases)) {
                    basesLengthDifference = thatRefBases.length - thisRefBases.length;
                    extensionBases = Arrays.copyOfRange(thatRefBases, thisRefBases.length, thatRefBases.length);
                } else if (Nucleotide.startsWith(thisRefBases, thatRefBases)) {
                    basesLengthDifference = thatRefBases.length - thatRefBases.length;
                    extensionBases = null;
                } else {
                    throw new IllegalArgumentException("unmatchable reference alleles: " + thisRefAllele + " vs " + thatRefAllele);
                }
            }
        } else {
            basesLengthDifference = 0;
            extensionBases = null;
        }
        final int numberOfAlleles = numberOfAlleles();
        final int[] result = new int[numberOfAlleles];
        if (basesLengthDifference == 0) {
            for (int i = 0; i < numberOfAlleles; i++) {
                result[i] = originalList.indexOfAllele(getAllele(i), true, true);
            }
        } else if (basesLengthDifference > 0) {
            for (int i = 0; i < numberOfAlleles; i++) {
                result[i] = originalList.indexOfAllele(Allele.extend(getAllele(i), extensionBases), true, true);
            }
        } else {
            for (int i = 0; i < numberOfAlleles; i++) {
                final Allele thisAllele = getAllele(i);
                final byte[] thisBases = thisAllele.getBases();
                final Allele searchAllele = thisAllele.isSymbolic() || thisAllele.isNoCall() ? thisAllele :
                        Allele.create(Arrays.copyOfRange(thisBases, 0, thisBases.length + basesLengthDifference));
                result[i] = originalList.indexOfAllele(Allele.create(searchAllele, thisAllele.isReference()), true, true);
            }
        }
        return new IntArrayList(result);
    }

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
    default boolean containsAllele(final A allele){
        return indexOfAllele(allele) >= 0;
    }

    @SuppressWarnings({"rawtypes"})
    static final AlleleList EMPTY_LIST = new AlleleList() {
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
        return (AlleleList<A>)EMPTY_LIST;
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

    static class MergeResult<A extends Allele> {
        public final AlleleList<A> left;
        public final AlleleList<A> right;
        public final AlleleList<A> result;
        public final IntList leftIndexesOnResult;
        public final IntList rightIndexesOnResult;

        private MergeResult(final AlleleList<A> left, final AlleleList<A> right, final AlleleList<A> result,
                            final IntList leftIndexesOnResult, final IntList rightIndexesOnResult) {
            this.left = left;
            this.right = right;
            this.result = result;
            this.leftIndexesOnResult = leftIndexesOnResult;
            this.rightIndexesOnResult = rightIndexesOnResult;
        }
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
    default public List<A> asListOfAlleles() {
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

        public NonPermutation(final AlleleList<A> original) {
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
        public int indexOfAllele(final A allele) {
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
        public int indexOfAllele(final A allele) {
            return to.indexOfAllele(allele);
        }

        @Override
        public A getAllele(final int index) {
            return to.getAllele(index);
        }
    }
}
