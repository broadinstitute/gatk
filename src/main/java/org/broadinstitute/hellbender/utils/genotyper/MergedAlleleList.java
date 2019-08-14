package org.broadinstitute.hellbender.utils.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import it.unimi.dsi.fastutil.ints.AbstractIntList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntLists;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAlleleCounts;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculator;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculators;
import org.broadinstitute.hellbender.utils.Log10Cache;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.ToDoubleBiFunction;
import java.util.function.ToDoubleFunction;
import java.util.function.ToIntFunction;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class MergedAlleleList<A extends Allele> implements AlleleList<A> {
    private final AlleleList<A> left;
    private final AlleleList<A> right;
    private final AlleleList<A> result;
    private final IntList indexOfLeftAlleles;
    private final IntList indexOfRightAlleles;
    private List<IntList[]> indexesOfLeftGenotypes;
    private List<IntList[]> indexesOfRightGenotypes;


    private final GenotypeLikelihoodCalculators GT_LK_CALCULATORS = new GenotypeLikelihoodCalculators();

    private MergedAlleleList(final AlleleList<A> left, final AlleleList<A> right, final AlleleList<A> result, final IntList indexOfLeftAllelesOnResult,
                             final IntList indexOfRightAllelesOnResult) {
        this.left = left;
        this.right = right;
        this.result = result;
        this.indexOfLeftAlleles = indexOfLeftAllelesOnResult;
        this.indexOfRightAlleles = indexOfRightAllelesOnResult;
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
    public int indexOfAllele(Allele allele) {
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

    public int[] mapIntPerAlleleAttribute(final AlleleList<A> oldAlleleList, final int[] oldValues, final int missingValue) {
        Utils.nonNull(oldAlleleList);
        Utils.nonNull(oldValues);
        Utils.validate(oldValues.length == oldAlleleList.numberOfAlleles(), "");
        final IntList indexOfOldAlleles = calculateIndexesOfAlleles(oldAlleleList);
        final int[] newValues = new int[result.numberOfAlleles()];
        for (int i = 0; i < newValues.length; i++) {
            final int oldIndex = indexOfOldAlleles.getInt(i);
            newValues[i] = oldIndex < 0 ? missingValue : oldValues[oldIndex];
        }
        return newValues;
    }

    /**
     * Merge likelihoods.
     * @param oldAlleleList
     * @param ploidy
     * @param likelihoods
     * @return
     */
    public double[] mapGenotypeLikelihoods(final AlleleList<A> oldAlleleList, final int ploidy, final double[] likelihoods) {
        final List<IntList> indexOfGenotypes = calculateIndexesOfGenotypes(oldAlleleList, ploidy);
        final double[] result = new double[indexOfGenotypes.size()];
        final double minLikelihood = MathUtils.arrayMin(likelihoods);
        for (int i = 0; i < result.length; i++) {
            final IntList oldGenotypes = indexOfGenotypes.get(i);
            final int oldGenotypesSize = oldGenotypes.size();
            if (oldGenotypesSize == 0 || oldGenotypesSize == likelihoods.length) {
                result[i] = minLikelihood - MathUtils.log10(ploidy);
            } else if (oldGenotypes.size() == 1) {
                result[i] = likelihoods[oldGenotypes.getInt(0)];
            } else {
                int minLikelihoodIndex = oldGenotypes.get(0);
                double minLikelihoodValue = likelihoods[oldGenotypes.get(minLikelihoodIndex)];
                for (int j = 1; j < oldGenotypes.size(); j++) {
                    int candidateIndex = oldGenotypes.get(j);
                    if (likelihoods[candidateIndex] < minLikelihoodValue) {
                        minLikelihoodValue = likelihoods[candidateIndex];
                        minLikelihoodIndex = candidateIndex;
                    }
                }
                final GenotypeLikelihoodCalculator oldCalculator = GT_LK_CALCULATORS.getInstance(ploidy, oldAlleleList.numberOfAlleles());
                final GenotypeAlleleCounts example = oldCalculator.genotypeAlleleCountsAt(minLikelihoodIndex);
                int matched = 0;
                for (int k = 0; k < example.distinctAlleleCount(); k++) {
                    if (this.result.containsAllele(oldAlleleList.getAllele(example.alleleIndexAt(k)))) {
                        matched += example.alleleCountAt(k);
                    }
                }
                result[i] = minLikelihoodValue + MathUtils.log10(matched) - MathUtils.log10(ploidy);
            }
        }
        return result;

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


    private static <A extends Allele, B extends Allele> IntList calculateIndexesOfAlleles(final AlleleList<A> originalAlleles, final AlleleList<B> resultAlleles) {
        final int resultLength = originalAlleles.numberOfAlleles();
        if (resultLength == 0) {
            return IntLists.EMPTY_LIST;
        } else if (resultLength == 1) {
            final int index = resultAlleles.indexOfAllele(originalAlleles.getAllele(0));
            if (index > 0) {
                return IntLists.singleton(index);
            } else {
                return IntLists.singleton(resultAlleles.nonRefAlleleIndex());
            }
        } else {
            final IntList result = new IntArrayList(resultLength);
            final int nonRefIndex = originalAlleles.nonRefAlleleIndex();
            for (int i = 0; i < resultLength; i++) {
                final int index = originalAlleles.indexOfAllele(resultAlleles.getAllele(i));
                result.add(index >= 0 ? index : nonRefIndex);
            }
            return result;
        }
    }

    private <A extends Allele> IntList calculateIndexesOfAlleles(final AlleleList<A> originalAlleles) {
        if (originalAlleles.sameAlleles(left)) {
            return indexOfLeftAlleles;
        } else if (originalAlleles.sameAlleles(right)) {
            return indexOfRightAlleles;
        } else {
            return calculateIndexesOfAlleles(originalAlleles, result);
        }
    }

    private <A extends Allele> List<IntList> calculateIndexesOfGenotypes(final AlleleList<A> original, final int ploidy) {
        final IntList indexOfOriginalAlleles = calculateIndexesOfAlleles(original);

        final GenotypeLikelihoodCalculator resultCalculator = GT_LK_CALCULATORS.getInstance(ploidy, result.numberOfAlleles());
        final int resultLength = resultCalculator.genotypeCount();
        if (resultLength == 0) {
            return Collections.emptyList();
        } else {
            final List<IntList> result = new ArrayList<>(resultLength);
            final int[] alleleIndexBuffer = new int[ploidy];
            int[] alleleCountsByIndexBuffer = null;
            final GenotypeLikelihoodCalculator originalCalculator = GT_LK_CALCULATORS.getInstance(ploidy, original.numberOfAlleles());
            for (int i = 0; i < resultLength; i++) {
                final GenotypeAlleleCounts resultAlleleCounts = resultCalculator.genotypeAlleleCountsAt(i);
                resultAlleleCounts.copyAlleleIndexes(alleleIndexBuffer, 0);
                boolean notFound = false;
                for (int j = 0; j < ploidy; j++) {
                    final int index = alleleIndexBuffer[j] = indexOfOriginalAlleles.getInt(alleleIndexBuffer[j]);
                    notFound |= index == -1;
                }
                if (!notFound) {
                    result.add(IntLists.singleton(originalCalculator.allelesToIndex(alleleIndexBuffer)));
                } else {
                    alleleCountsByIndexBuffer = alleleCountsByIndexBuffer == null ? alleleCountsByIndexBuffer : new int[original.numberOfAlleles()];
                    final IntList resultList = new IntArrayList(20); // 20 will cover most small cases.
                    final int originalGenotypeCount = originalCalculator.genotypeCount();
                    for (final int alleleIndex : alleleIndexBuffer) {
                        if (alleleIndex >= 0) {
                            alleleCountsByIndexBuffer[alleleIndex]++;
                        }
                    }
                    for (int j = 0; j < originalGenotypeCount; j++) {
                        final GenotypeAlleleCounts originalGenotypeAlleleCounts = originalCalculator.genotypeAlleleCountsAt(j);
                        if (originalGenotypeAlleleCounts.greaterOrEqualAlleleCountsByIndex(alleleCountsByIndexBuffer)) {
                            resultList.add(j);
                        }
                    }
                    Arrays.fill(alleleCountsByIndexBuffer, 0);
                    result.add(resultList);
                }
            }
            return result;
        }
    }

}
