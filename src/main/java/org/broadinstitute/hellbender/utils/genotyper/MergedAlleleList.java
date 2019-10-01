package org.broadinstitute.hellbender.utils.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import it.unimi.dsi.fastutil.ints.AbstractIntList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntLists;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class MergedAlleleList<A extends Allele> implements AlleleList<A> {
    private final AlleleList<A> left;
    private final AlleleList<A> right;
    private final AlleleList<A> result;
    private final IntList indexOfLeftAlleles;
    private final IntList indexOfRightAlleles;
    private List<List<GenotypeMatch>> indexesOfLeftGenotypes;
    private List<List<GenotypeMatch>> indexesOfRightGenotypes;


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


    public <T> List<T> mapAlleleAnnotation(final AlleleList<?> alleles, final List<T> values, final VCFHeaderLineCount count, final T missingValue) {
        switch (count) {
            case INTEGER:
            case UNBOUNDED:
                return values;
            case A:
                return mapAlleleAnnotation(alleles, values, false, missingValue);
            case R:
                return mapAlleleAnnotation(alleles, values, true, missingValue);
            default:
                throw new IllegalArgumentException("unsupported annotation VCF count type: " + count);
        }
    }

    public int[] mapAlleleAnnotation(final AlleleList<?> alleles, final int[] values, final VCFHeaderLineCount count, final int missingValue) {
        final List<Integer> valueList = new AbstractList<Integer>() {
            @Override
            public Integer get(int index) {
                return values[index];
            }

            @Override
            public int size() {
                return values.length;
            }
        };
        final List<Integer> result = mapAlleleAnnotation(alleles, valueList, count, missingValue);
        return result.stream().mapToInt(Integer::intValue).toArray();
    }

    public <T> List<T> mapAlleleAnnotation(final AlleleList<?> oldAlleleList, final List<T> oldValues, final boolean includeReference, final T missingValue) {
        Utils.nonNull(oldAlleleList);
        Utils.nonNull(oldValues);
        int referenceIndex = oldAlleleList.indexOfUniqueReference();
        if (referenceIndex == -1 && oldAlleleList.indexOfReference() != -1) {
            referenceIndex = oldAlleleList.indexOfAllele(result.getReferenceAllele(), true, false);
        }
        if (includeReference) {
            Utils.validate(oldValues.size() == oldAlleleList.numberOfAlleles(), "bad input length");
            return unsafeMapAlleleAnnotation(oldAlleleList, oldValues, true, missingValue);
        } else if (referenceIndex == -1) {
            throw new IllegalArgumentException("cannot determine the reference allele");
        } else {
            if (oldValues.size() == oldAlleleList.numberOfAlleles() + 1) { // the input does not have a value for the ref.
                final List<T> paddedReferenceValue = new ArrayList<>(oldValues.size() + 1);
                paddedReferenceValue.addAll(oldValues.subList(0, referenceIndex));
                paddedReferenceValue.add(null);
                paddedReferenceValue.addAll(oldValues.subList(referenceIndex, oldValues.size()));
                return unsafeMapAlleleAnnotation(oldAlleleList, paddedReferenceValue, false, missingValue);
            } else if (oldValues.size() == oldAlleleList.numberOfAlleles()) { // the input contains a value for the ref.
                return unsafeMapAlleleAnnotation(oldAlleleList, oldValues, false, missingValue);
            } else {
                throw new IllegalArgumentException("wrong number of values");
            }
        }
    }

    private <T> List<T>  unsafeMapAlleleAnnotation(final AlleleList<?> oldAlleleList, final List<T> oldValues, final boolean includeReference, final T missingValue) {
        final IntList indexOfOldAlleles = result.alleleIndexMap(oldAlleleList, true, true);
        final int numberOfAlleles = numberOfAlleles();
        final List<T> result = new ArrayList<>(numberOfAlleles);
        for (int i = 0; i < numberOfAlleles; i++) {
            final int oldIndex = indexOfOldAlleles.getInt(i);
            result.add(oldIndex < 0 ? missingValue : oldValues.get(oldIndex));
        }
        return result;
    }

    /**
     * Maps likelihoods on a different allele list as they have been calculated on
     * the merged allele list.
     *
     * <p>
     *     The input source allele list could be one of the merged list or a different one.
     * </p>
     * <p>
     *     When merged list alleles cannot be mapped to the source allele one we try two strategies to
     *     rescue the affected likelihoods:
     *     <p>1. If the input source list has the special {@code <NON_REF>} allele we use it in-leu for
     *     all missing alleles in the merged list.
     *     </p>
     *     <p>2. If {@code <NON_REF>} is not present we calculate the genotype likelihoods based on
     *     the likelihoods of all the input genotypes that is shares one allele with using the following formula:
     *      <pre>
     *          lk = sum_i(lk[i] ^ 2 * #shared_alleles(gt[i], gt) / ploidy ) / sum_i(lk[i])
     *          i in [0,#source_gts] where #shared_alleels(gt, gt[i]) > 0
     *          lk[i] the ith source likelihood.
     *          #shared_alleles(a, b): number of alleles call shared between both genotypes.
     *      </pre>
     *     </p>
     *     <p>
     *         In the second approach we try to calculate the likelihood of a genotype as a weighted sum
     *         of the source genotypes where more likely have a greater weight. Then contribution of each
     *         of those source genotypes is "degraded" based on the number of non-shared allele calls.
     *     </p>
     *
     *     <p>
     *         Example merged list has A*,B,C whereas input source list has A*,B. So C is not present
     *         in the input:
     *
     *         lk_out(A*,A*) = lk_in(A*,A*) # no problem.
     *         lk_out(A*,B) = lk_in(A*,B) ...
     *         ...
     *         lk_out(A*,C)  = combine(lk_in(A*,A*), lk_in(A*,B)) using formula
     *                       = (lk_in(A*,A*)^2 * 0.5 + lk_in(A*,B)^2 * 0.5) /
     *                         (lk_in(A*,A*) + lk_in(A*,B))
     *         lk_out(B,C)   = (lk_in(B,B)^2 * 0.5 + lk_in(A*,B)^2 * 0.5) /
     *                               (lk_in(B,B) + lk_in(A*,B))
     *         lk_out(C,C)   = zero # as calculated below:
     *     </p>
     *
     *     <p>
     *         We set the lowest possible (zero) likelihood as the lowest input likelihood divide by 2 or
     *         the the lowest input likelihood divided by its ratio vs the second lowest whichever is less BUT
     *         never less than the ratio between the worst and the best likelihoods:
     *         So for input 0, 10, 20. It would be 30 (20 + (20 - 10)). For input 0, 10, 10 would be 13.01 (10 + 3.01 (=log10(2)*-10))
     *         for input 0, 0, 0 it would be 0 (using the third rule).
     *
     *     </p>
     *
     * </p>
     *
     * @param sourceAlleles the source allele list.
     * @param ploidy the ploidy for the sample involved.
     * @param sourceLikelihoods the source likelihoods.
     * @return never {@code null}, an array with the corresponding likelihoods for the
     *   merged list.
     */
    public double[] mapGenotypeLikelihoods(final AlleleList<A> sourceAlleles, final int ploidy, final double[] sourceLikelihoods) {
        final List<List<GenotypeMatch>> genotypeMap = calculateGenotypeMap(sourceAlleles, ploidy);
        final int resultLength = genotypeMap.size();
        final double[] resultLikelihoods = new double[resultLength];
        final double[] normalizedLikelihoods = normalizeLikelihoods(sourceLikelihoods, false);
        final double zeroLikelihood = calculateZeroLikelihood(normalizedLikelihoods);
        double[] log10SumBuffer = null;  // may not be needed so lazily instantiated.
        double[] log10SumBuffer2 = null; // may not be needed so lazily instantiated.
        for (int i = 0; i < resultLength; i++) {
            final List<GenotypeMatch> oldGenotypes = genotypeMap.get(i);
            final int oldGenotypesSize = oldGenotypes.size();
            if (oldGenotypesSize == 0) {
                resultLikelihoods[i] = zeroLikelihood;
            } else if (oldGenotypesSize == 1) {
                final GenotypeMatch match = oldGenotypes.get(0);
                resultLikelihoods[i] = normalizedLikelihoods[match.genotypeIndex] + match.log10Factor;
            } else {
                log10SumBuffer = log10SumBuffer == null ? new double[sourceLikelihoods.length] : log10SumBuffer;
                log10SumBuffer2 = log10SumBuffer2 == null ? new double[sourceLikelihoods.length] : log10SumBuffer2;
                for (int j = 0; j < oldGenotypesSize; j++) {
                    final GenotypeMatch match = oldGenotypes.get(j);
                    log10SumBuffer[j] = 2 * normalizedLikelihoods[match.genotypeIndex];
                    log10SumBuffer2[j] = normalizedLikelihoods[match.genotypeIndex];
                    log10SumBuffer[j] += match.log10Factor;
                }
                final double normalizer = MathUtils.log10SumLog10(log10SumBuffer2, 0, oldGenotypesSize);
                resultLikelihoods[i] = Math.max(MathUtils.log10SumLog10(log10SumBuffer, 0, oldGenotypesSize) - normalizer, zeroLikelihood);
            }
        }
        return normalizeLikelihoods(resultLikelihoods, true);
    }

    private double[] normalizeLikelihoods(final double[] likelihoods, final boolean canUseInput) {

        final int maxIndex = MathUtils.maxElementIndex(likelihoods);
        if (likelihoods[maxIndex] == 0.0) {
            return likelihoods;
        } else {
            final int length = likelihoods.length;
            final double max = likelihoods[maxIndex];
            final double[] result = canUseInput ? likelihoods : new double[length];
            for (int i = 0; i < length; i++) {
                result[i] = likelihoods[i] - max;
            }
            return result;
        }
    }

    private double calculateZeroLikelihood(final double[] likelihoods) {
        double worst = likelihoods[0];
        double secondWorst = worst;
        final int length = likelihoods.length;
        for (int i = 1; i < length; i++) {
            final double lk = likelihoods[i];
            if (lk < worst) {
                secondWorst = worst;
                worst = lk;
            } else if (lk < secondWorst) {
                secondWorst = lk;
            }
        }
        final double diff = secondWorst - worst;
        return diff < MathUtils.log10(2) ? worst - MathUtils.log10(2) : worst - diff;
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
        final int resultLength = resultAlleles.numberOfAlleles();
        if (resultLength == 0) {
            return IntLists.EMPTY_LIST;
        } else if (resultLength == 1) {
            return IntLists.singleton(originalAlleles.indexOfAllele(resultAlleles.getAllele(0), resultAlleles.getUniqueReferenceAllele(), true, true));
        } else {
            final Allele resultReference = resultAlleles.getUniqueReferenceAllele();
            final IntList result = new IntArrayList(resultLength);
            for (int i = 0; i < resultLength; i++) {
                result.add(originalAlleles.indexOfAllele(resultAlleles.getAllele(i), resultReference, true, true));
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

    private static class GenotypeMatch {
        private final int genotypeIndex;
        private final int matchedAlleles;
        private final double log10Factor;

        private GenotypeMatch(final int genotypeIndex, final int matchedAlleles, final double factor) {
           this.genotypeIndex = genotypeIndex;
           this.matchedAlleles = matchedAlleles;
           this.log10Factor = factor;
        }
    }

    private <A extends Allele> List<List<GenotypeMatch>> calculateGenotypeMap(final AlleleList<A> original, final int ploidy) {
        final IntList indexOfOriginalAlleles = calculateIndexesOfAlleles(original);

        final GenotypeLikelihoodCalculator resultCalculator = GT_LK_CALCULATORS.getInstance(ploidy, result.numberOfAlleles());
        final int resultLength = resultCalculator.genotypeCount();
        if (resultLength == 0) {
            return Collections.emptyList();
        } else {
            final List<List<GenotypeMatch>> result = new ArrayList<>(resultLength);
            final int[] alleleIndexBuffer = new int[ploidy];
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
                    result.add(Collections.singletonList(new GenotypeMatch(originalCalculator.allelesToIndex(alleleIndexBuffer), ploidy, 0)));
                } else {
                    final Int2IntMap allelesIndicesPresent = new Int2IntArrayMap(originalCalculator.alleleCount());
                    allelesIndicesPresent.defaultReturnValue(0);
                    for (int k = 0; k < ploidy; k++) {
                       if (alleleIndexBuffer[k] >= 0) {
                           allelesIndicesPresent.put(alleleIndexBuffer[k], allelesIndicesPresent.get(alleleIndexBuffer[k]) + 1);
                       }
                    }
                    if (allelesIndicesPresent.isEmpty()) {
                        result.add(Collections.emptyList());
                    } else {
                        final int[] uniqueAlleleIndicesPresent = allelesIndicesPresent.keySet().toIntArray();
                        final IntList intersectingGenotypeIndexes = originalCalculator.allelesContainingIndexes(uniqueAlleleIndicesPresent);
                        final List<GenotypeMatch> matches = new ArrayList<>(intersectingGenotypeIndexes.size());
                        for (final int index : intersectingGenotypeIndexes) {
                            int matched = 0;
                            final GenotypeAlleleCounts counts = originalCalculator.genotypeAlleleCountsAt(index);
                            for (int alleleRank = 0; alleleRank < counts.distinctAlleleCount(); alleleRank++) {
                                final int alleleIndex = counts.alleleIndexAt(alleleRank);
                                final int alleleCount = counts.alleleCountAt(alleleRank);
                                matched += Math.min(alleleCount, allelesIndicesPresent.get(alleleIndex));
                            }
                            matches.add(new GenotypeMatch(index, matched, Math.log10(matched) - Math.log10(ploidy)));
                        }
                        result.add(matches);
                    }
                }
            }
            return result;
        }
    }

}
