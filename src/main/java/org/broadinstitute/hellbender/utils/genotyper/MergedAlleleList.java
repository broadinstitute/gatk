package org.broadinstitute.hellbender.utils.genotyper;

import htsjdk.variant.variantcontext.Allele;
import it.unimi.dsi.fastutil.ints.AbstractIntList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class MergedAlleleList<A extends Allele> implements AlleleList<A> {
    private final AlleleList<A> left;
    private final AlleleList<A> right;
    private final AlleleList<A> result;
    private final IntList indexOfLeftAllelesOnResult;
    private final IntList indexOfRightAllelesOnResult;


    private MergedAlleleList(final AlleleList<A> left, final AlleleList<A> right, final AlleleList<A> result, final IntList indexOfLeftAllelesOnResult,
                             final IntList indexOfRightAllelesOnResult) {
        this.left = left;
        this.right = right;
        this.result = result;
        this.indexOfLeftAllelesOnResult = indexOfLeftAllelesOnResult;
        this.indexOfRightAllelesOnResult = indexOfRightAllelesOnResult;
    }

    @SuppressWarnings("unchecked")
    private static <A extends Allele> A extendWithBasesIfPossible(final A original, final byte[] bases) {
        if (original.getClass() != Allele.class) {
            return original;
        }
        for (final byte base : original.getBases()) {
            switch (base) {
                case 'a':
                case 'A':
                case 'c':
                case 'C':
                case 't':
                case 'T':
                case 'n':
                case 'N':
                    continue;
                default:
                    return original;
            }
        }
        return (A) Allele.create(Utils.concat(original.getBases(), bases), original.isReference());
    }

    @Override
    public int numberOfAlleles() {
        return result.numberOfAlleles();
    }

    @Override
    public int indexOfAllele(A allele) {
        return result.indexOfAllele(allele);
    }

    @Override
    public A getAllele(int index) {
        return result.getAllele(index);
    }

    public static <A extends Allele> MergedAlleleList<A> merge(final AlleleList<A> left, final AlleleList<A> right) {
        if (left == right) { // may happen often in some cases.
            final IntList indexOnList = new FirstNIntegers(left.numberOfAlleles());
            return new MergedAlleleList<>(left, left, left, indexOnList, indexOnList);
        } else {
            final int leftSize = left.numberOfAlleles();
            int leftRefIdx = -1;
            int rightRefIdx = -1;
            for (int i = 0; i < leftSize; i++) {
                final A allele = left.getAllele(i);
                if (allele.isReference()) {
                    if (leftRefIdx != -1) {
                        throw new IllegalArgumentException("more than one reference in input left allele list");
                    } else {
                        leftRefIdx = i;
                    }
                }
            }
            final int rightSize = right.numberOfAlleles();
            boolean sameList = leftSize == rightSize;
            for (int j = 0; j < rightSize; j++) {
                final A allele = right.getAllele(j);
                if (allele.isReference()) {
                    if (rightRefIdx != -1) {
                        throw new IllegalArgumentException("more than one reference in input left allele list");
                    } else {
                        rightRefIdx = j;
                    }
                }
                sameList &= allele.equals(left.getAllele(j));
            }
            if (sameList) {
                final IntList indexOnList = new FirstNIntegers(leftSize);
                return new MergedAlleleList<>(left, left, left, indexOnList, indexOnList);
            } else {
                return mergeDifferentAlleleList(left, right, leftRefIdx, rightRefIdx);
            }
        }
    }

    private static <A extends Allele> MergedAlleleList<A> mergeDifferentAlleleList(final AlleleList<A> originalLeft, final AlleleList<A> originalRight,
                                                                final int leftRefIdx, final int rightRefIdx) {
        final AlleleList<A> left;
        final AlleleList<A> right;
        if (leftRefIdx >= 0 && rightRefIdx >= 0) {
            final A leftRef = originalLeft.getAllele(leftRefIdx);
            final A rightRef = originalRight.getAllele(rightRefIdx);
            if (leftRef.equals(rightRef)) {
                left = originalLeft;
                right = originalRight;
            } else if (leftRef.getDisplayString().startsWith(rightRef.getDisplayString())) {
                left = originalLeft;
                final byte[] extension = Arrays.copyOfRange(leftRef.getBases(), rightRef.getBases().length, leftRef.getBases().length);
                right = new IndexedAlleleList<>(originalRight.asListOfAlleles().stream()
                        .map(a -> extendWithBasesIfPossible(a, extension))
                        .collect(Collectors.toList()));
            } else if (rightRef.getDisplayString().startsWith(leftRef.getDisplayString())) {
                right = originalRight;
                final byte[] extension = Arrays.copyOfRange(rightRef.getBases(), leftRef.getBases().length, rightRef.getBases().length);
                left = new IndexedAlleleList<>(originalLeft.asListOfAlleles().stream()
                        .map(a -> extendWithBasesIfPossible(a, extension))
                        .collect(Collectors.toList()));
            } else {
                throw new IllegalArgumentException("reference alleles are not compatible: " + leftRef + " " + rightRef);
            }
        } else {
            left = originalLeft;
            right = originalRight;
        }
        final AlleleList<A> result = new IndexedAlleleList<A>(Stream.concat(left.asListOfAlleles().stream(),
                right.asListOfAlleles().stream())
                .sorted().distinct().collect(Collectors.toList()));
        final IntList leftIndexes = new IntArrayList(left.asListOfAlleles().stream().mapToInt(result::indexOfAllele).toArray());
        final IntList rightIndexes = new IntArrayList(right.asListOfAlleles().stream().mapToInt(result::indexOfAllele).toArray());
        return new MergedAlleleList<>(originalLeft, originalRight, result, leftIndexes, rightIndexes);
    }


    private static class FirstNIntegers extends AbstractIntList {

        private final int size;

        private FirstNIntegers(final int n) {
            size = n;
        }

        @Override
        public int getInt(int i) {
            if (i < 0 || size <= i) {
                throw new IndexOutOfBoundsException();
            } else {
                return i;
            }
        }

        @Override
        public int size() {
            return size;
        }
    }


}
