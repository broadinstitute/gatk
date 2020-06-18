package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.apache.commons.lang.mutable.MutableDouble;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.PriorityQueue;
import java.util.function.Consumer;

/**
 * Helper to calculate genotype likelihoods given a ploidy and an allele count (number of possible distinct alleles).
 */
public final class GenotypeLikelihoodCalculator implements Iterable<GenotypeAlleleCounts>{

    /**
     * Offset table for this calculator.
     *
     * <p>
     *     This is a shallow copy of {@link GenotypeLikelihoodCalculators#alleleFirstGenotypeOffsetByPloidy} when the calculator was created
     *     thus it follows the same format as that array. Please refer to its documentation.
     * </p>
     *
     * <p>You can assume that this offset table contain at least (probably more) the numbers corresponding to the allele count and ploidy for this calculator.
     * However since it might have more than that and so you must use {@link #alleleCount} and {@link #ploidy} when
     * iterating through this array rather that its length or the length of its components.</p>.
     */
    private final int[][] alleleFirstGenotypeOffsetByPloidy;

    /**
     * Genotype table for this calculator.
     *
     * <p>It is ensure that it contains all the genotypes for this calculator ploidy and allele count, maybe more. For
     * that reason you must use {@link #genotypeCount} when iterating through this array and not relay on its length.</p>
     */
    private final GenotypeAlleleCounts[] genotypeAlleleCounts;

    /**
     * Number of genotypes given this calculator {@link #ploidy} and {@link #alleleCount}.
     */
    private final int genotypeCount;

    private final int alleleCount;
    private final int ploidy;

    /**
     * Max-heap for integers used for this calculator internally.
     */
    private final PriorityQueue<Integer> alleleHeap;

    /**
     * Cache of the last genotype-allele-count requested using {@link #genotypeAlleleCountsAt(int)}, when it
     * goes beyond the maximum genotype-allele-count static capacity. Check on that method documentation for details.
     */
    private GenotypeAlleleCounts lastOverheadCounts;

    private static final int INITIAL_READ_CAPACITY = 10;
    /**
     * Indicates how many reads the calculator supports.
     *
     * This figure is increased dynamically by {@link #ensureReadCapacity(int) ensureReadCapacity}.
     */
    private int readCapacity = INITIAL_READ_CAPACITY;

    /**
     * Buffer field used as a temporary container when one value per read is stored.
     *
     * In the multiallelic calculation we accumulate the likelihood contribution of each read one allele at a time.  That is,
     * for genotype containing alleles A, B, C, we first fill the buffer with the likelihood contributions from allele A, then
     * we make a second pass and add the contributions from allele B, then allele C.  Traversing each allele row of the
     * likelihoods array in this manner is cache-friendly.
     */
    private double[] perReadBuffer = new double[INITIAL_READ_CAPACITY];

    /**
     * Creates a new calculator providing its ploidy and number of genotyping alleles.
     */
    protected GenotypeLikelihoodCalculator(final int ploidy, final int alleleCount,
                                           final int[][] alleleFirstGenotypeOffsetByPloidy,
                                           final GenotypeAlleleCounts[][] genotypeTableByPloidy) {
        Utils.validateArg(ploidy > 0, () -> "ploidy must be at least 1 but was " + ploidy);
        this.alleleFirstGenotypeOffsetByPloidy = alleleFirstGenotypeOffsetByPloidy;
        genotypeAlleleCounts = genotypeTableByPloidy[ploidy];
        this.alleleCount = alleleCount;
        this.ploidy = ploidy;
        genotypeCount = this.alleleFirstGenotypeOffsetByPloidy[ploidy][alleleCount];
        alleleHeap = new PriorityQueue<>(ploidy, Comparator.<Integer>naturalOrder().reversed());
    }

    /**
     * Makes sure that temporal arrays and matrices are prepared for a number of reads to process.
     * @param requestedCapacity number of read that need to be processed.
     */
    private void ensureReadCapacity(final int requestedCapacity) {
        if (readCapacity < requestedCapacity) {
            readCapacity = 2 * requestedCapacity;
            perReadBuffer = new double[readCapacity];
        }
    }

    /**
     * Give a list of alleles, returns the likelihood array index.
     * @param alleleIndices the indices of the alleles in the genotype, there should be as many repetition of an
     *                      index as copies of that allele in the genotype. Allele indices do not need to be sorted in
     *                      any particular way.
     *
     * @return never {@code null}.
     */
    public int allelesToIndex(final int... alleleIndices) {
        // Special case ploidy == 0.
        if (ploidy == 0) {
            return 0;
        }

        alleleHeap.clear();
        for (int i = 0; i < alleleIndices.length; i++) {
            alleleHeap.add(alleleIndices[i]);
        }
        return alleleHeapToIndex();
    }

    /**
     * Returns the number of possible genotypes given ploidy and the maximum allele index.
     * @return never {@code null}.
     */
    public int genotypeCount()  {
        return genotypeCount;
    }

    /**
     * Returns the genotype associated to a particular likelihood index.
     *
     * <p>If {@code index} is larger than {@link GenotypeLikelihoodCalculators#MAXIMUM_NUMBER_OF_CACHED_GENOTYPE_ALLELE_COUNTS_PER_CALCULATOR},
     *  this method will reconstruct that genotype-allele-count iteratively from the largest strongly referenced count available.
     *  or the last requested index genotype.
     *  </p>
     *
     * <p> Therefore if you are iterating through all genotype-allele-counts you should do sequentially and incrementally, to
     * avoid a large efficiency drop </p>.
     *
     * @param index query likelihood-index.
     * @return never {@code null}.
     */
    public GenotypeAlleleCounts genotypeAlleleCountsAt(final int index) {
        Utils.validateArg(index >= 0 && index < genotypeCount, () -> "invalid likelihood index: " + index + " >= " + genotypeCount
                    + " (genotype count for nalleles = " + alleleCount + " and ploidy " + ploidy);
        if (index < GenotypeLikelihoodCalculators.MAXIMUM_NUMBER_OF_CACHED_GENOTYPE_ALLELE_COUNTS_PER_CALCULATOR) {
            return genotypeAlleleCounts[index];
        } else if (lastOverheadCounts == null || lastOverheadCounts.index() > index) {
            final GenotypeAlleleCounts result = genotypeAlleleCounts[GenotypeLikelihoodCalculators.MAXIMUM_NUMBER_OF_CACHED_GENOTYPE_ALLELE_COUNTS_PER_CALCULATOR - 1].copy();
            result.increase(index - GenotypeLikelihoodCalculators.MAXIMUM_NUMBER_OF_CACHED_GENOTYPE_ALLELE_COUNTS_PER_CALCULATOR + 1);
            lastOverheadCounts = result;
            return result.copy();
        } else {
            lastOverheadCounts.increase(index - lastOverheadCounts.index());
            return lastOverheadCounts.copy();
        }
    }

    /**
     * Calculate the log10Likelihoods given the list of alleles and the likelihood map.
     *
     * @param log10Likelihoods the likelihood matrix all alleles vs all reads.
     *
     * @throws IllegalArgumentException if {@code alleleList} is {@code null} or {@code log10Likelihoods} is {@code null}
     *     or the alleleList size does not match the allele-count of this calculator, or there are missing allele vs
     *     read combinations in {@code log10Likelihoods}.
     *
     * @return never {@code null}.
     */
    public <EVIDENCE, A extends Allele> GenotypeLikelihoods genotypeLikelihoods(final LikelihoodMatrix<EVIDENCE, A> log10Likelihoods) {
        Utils.nonNull(log10Likelihoods);
        Utils.validateArg(log10Likelihoods.numberOfAlleles() == alleleCount, "mismatch between allele list and alleleCount");
        final int readCount = log10Likelihoods.evidenceCount();
        ensureReadCapacity(readCount);

        final double[][] log10LikelihoodsByAlleleAndRead = log10Likelihoods.copyAlleleLikelihoods();
        final boolean triallelicGenotypesPossible = alleleCount > 2 && ploidy > 2;

        // non-log space likelihoods for multiallelic computation
        final double[][] likelihoodsByAlleleAndRead = triallelicGenotypesPossible ? log10Likelihoods.copyAlleleLikelihoods() : null;
        final MutableDouble multiallelicNormalization = new MutableDouble(0);
        if (triallelicGenotypesPossible) {
            makeNormalizedNonLogLikelihoods(readCount, likelihoodsByAlleleAndRead, multiallelicNormalization);
        }

        final double[] result = new double[genotypeCount];

        Utils.stream(iterator()).forEach(alleleCounts -> {
            final int componentCount = alleleCounts.distinctAlleleCount();
            final int genotypeIndex = alleleCounts.index();
            if (componentCount == 1) {
                // homozygous case: log P(reads|AAAAA. . .) = sum_{reads} log P(read|A)
                final int allele = alleleCounts.alleleIndexAt(0);
                result[genotypeIndex] = MathUtils.sum(log10LikelihoodsByAlleleAndRead[allele]);
            } else if (componentCount == 2) {
                // biallelic het case: log P(reads | nA copies of A, nB copies of B) = sum_{reads} log[(nA * P(read | A) + nB * P(read | B))] -log(ploidy)
                final double[] logLks1 = log10LikelihoodsByAlleleAndRead[alleleCounts.alleleIndexAt(0)];
                final int freq1 = alleleCounts.alleleCountAt(0);
                final double log10Freq1 = MathUtils.log10(freq1);
                final double[] logLks2  = log10LikelihoodsByAlleleAndRead[alleleCounts.alleleIndexAt(1)];
                final double log10Freq2 = MathUtils.log10(ploidy - freq1);

                result[genotypeIndex] = new IndexRange(0, readCount).sum(r -> MathUtils.approximateLog10SumLog10(logLks1[r] + log10Freq1, logLks2[r] + log10Freq2))
                        - readCount * MathUtils.log10(ploidy);
            } else {
                // the multiallelic case is conceptually the same as the biallelic case but done in non-log space
                // We implement in a cache-friendly way by adding nA * P(read|A) to the per-read buffer for all reads (inner loop), and all alleles (outer loop)
                Arrays.fill(perReadBuffer,0, readCount, 0);
                alleleCounts.forEachAlleleIndexAndCount((a, f) -> new IndexRange(0, readCount).forEach(r -> perReadBuffer[r] += f * likelihoodsByAlleleAndRead[a][r]));
                result[genotypeIndex] = new IndexRange(0, readCount).sum(r -> FastMath.log10(perReadBuffer[r])) - readCount * MathUtils.log10(ploidy) + multiallelicNormalization.doubleValue();
            }
        });
        return GenotypeLikelihoods.fromLog10Likelihoods(result);
    }

    private void makeNormalizedNonLogLikelihoods(int readCount, double[][] likelihoodsByAlleleAndRead, final MutableDouble multiallelicNormalization) {
        // fill the per-read buffer with the maximum log-likelihood per read
        Arrays.fill(perReadBuffer, 0, readCount, Double.NEGATIVE_INFINITY);
        for (int a = 0; a < alleleCount; a++) {
            for (int r = 0; r < readCount; r++) {
                perReadBuffer[r] = FastMath.max(perReadBuffer[r], likelihoodsByAlleleAndRead[a][r]);
            }
        }

        // subtract these maxima
        for (int a = 0; a < alleleCount; a++) {
            for (int r = 0; r < readCount; r++) {
                likelihoodsByAlleleAndRead[a][r] -= perReadBuffer[r];
            }
        }
        // switch to non-log now that we have moved to numerically safe zero-max space
        new IndexRange(0, alleleCount).forEach(a -> MathUtils.applyToArrayInPlace(likelihoodsByAlleleAndRead[a], x -> Math.pow(10.0, x)));

        // the normalization is the sum of all the subtracted-off maximum log likelihoods
        multiallelicNormalization.setValue(MathUtils.sum(perReadBuffer, 0, readCount));
    }

    @Override
    public Iterator<GenotypeAlleleCounts> iterator() {
        return new Iterator<GenotypeAlleleCounts>() {
            private int genotypeIndex = 0;
            private final GenotypeAlleleCounts increasingAlleleCounts = genotypeCount < GenotypeLikelihoodCalculators.MAXIMUM_NUMBER_OF_CACHED_GENOTYPE_ALLELE_COUNTS_PER_CALCULATOR ?
                    null : genotypeAlleleCounts[GenotypeLikelihoodCalculators.MAXIMUM_NUMBER_OF_CACHED_GENOTYPE_ALLELE_COUNTS_PER_CALCULATOR - 1].copy();

            @Override
            public boolean hasNext() {
                return genotypeIndex < genotypeCount;
            }

            @Override
            public GenotypeAlleleCounts next() {
                if (genotypeIndex < GenotypeLikelihoodCalculators.MAXIMUM_NUMBER_OF_CACHED_GENOTYPE_ALLELE_COUNTS_PER_CALCULATOR) {
                    return genotypeAlleleCounts[genotypeIndex++];
                } else {
                    increasingAlleleCounts.increase();
                    genotypeIndex++;
                    return increasingAlleleCounts;
                }
            }
        };
    }

    /**
     * Returns the ploidy for this genotype likelihood calculator.
     * @return 0 or greater.
     */
    public int ploidy() {
        return ploidy;
    }

    /**
     * Returns the total number of alleles for this genotype calculator.
     * @return the number of alleles considered by this calculator.
     */
    public int alleleCount() {
        return alleleCount;
    }

    /**
     * Returns the likelihood index given the allele counts.
     *
     * @param alleleCountArray the query allele counts. This must follow the format returned by
     *  {@link GenotypeAlleleCounts#copyAlleleCounts} with 0 offset.
     *
     * @throws IllegalArgumentException if {@code alleleCountArray} is not a valid {@code allele count array}:
     *  <ul>
     *      <li>is {@code null},</li>
     *      <li>or its length is not even,</li>
     *      <li>or it contains any negatives,
     *      <li>or the count sum does not match the calculator ploidy,</li>
     *      <li>or any of the alleles therein is negative or greater than the maximum allele index.</li>
     *  </ul>
     *
     * @return 0 or greater but less than {@link #genotypeCount}.
     */
    public int alleleCountsToIndex(final int ... alleleCountArray) {
        Utils.nonNull(alleleCountArray, "the allele counts cannot be null");
        Utils.validateArg((alleleCountArray.length & 1) == 0, "the allele counts array cannot have odd length");
        alleleHeap.clear();
        for (int i = 0; i < alleleCountArray.length; i += 2) {
            final int index = alleleCountArray[i];
            final int count = alleleCountArray[i+1];
            Utils.validateArg(count >= 0, "no allele count can be less than 0");
            for (int j = 0; j < count; j++) {
                alleleHeap.add(index);
            }
        }
        return alleleHeapToIndex();
    }

    /**
     * Transforms the content of the heap into an index.
     *
     * <p>
     *     The heap contents are flushed as a result, so is left ready for another use.
     * </p>
     *
     * @return a valid likelihood index.
     */
    private int alleleHeapToIndex() {
        Utils.validateArg(alleleHeap.size() == ploidy, "the sum of allele counts must be equal to the ploidy of the calculator");
        Utils.validateArg(alleleHeap.peek() < alleleCount, () -> "invalid allele " + alleleHeap.peek() + " more than the maximum " + (alleleCount - 1));
        int result = 0;
        for (int p = ploidy; p > 0; p--) {
            final int allele = alleleHeap.remove();
            Utils.validateArg(allele >= 0, () -> "invalid allele " + allele + " must be equal or greater than 0 ");
            result += alleleFirstGenotypeOffsetByPloidy[p][allele];
        }
        return result;
    }

    /**
     * Composes a genotype index map given a allele index recoding.
     *
     * @param oldToNewAlleleIndexMap allele recoding. The ith entry indicates the index of the allele in original encoding
     *                               that corresponds to the ith allele index in the final encoding.
     *
     * @throws IllegalArgumentException if this calculator cannot handle the recoding provided. This is
     * the case when either {@code oldToNewAlleleIndexMap}'s length or any of its element (+ 1 as they are 0-based) is larger
     * this calculator's {@link #alleleCount()}. Also if any {@code oldToNewAllelesIndexMap} element is negative.
     *
     * @return never {@code null}.
     */
    public int[] genotypeIndexMap(final int[] oldToNewAlleleIndexMap, final GenotypeLikelihoodCalculators calculators) {
        Utils.nonNull(oldToNewAlleleIndexMap);
        final int resultAlleleCount = oldToNewAlleleIndexMap.length;
        Utils.validateArg(resultAlleleCount <= alleleCount, () -> "this calculator does not have enough capacity for handling "
                    + resultAlleleCount + " alleles ");
        final int resultLength = resultAlleleCount == alleleCount
                ? genotypeCount : calculators.genotypeCount(ploidy,resultAlleleCount);

        final int[] result = new int[resultLength];
        final int[] sortedAlleleCounts = new int[Math.max(ploidy, alleleCount) << 1];
        alleleHeap.clear();
        Utils.stream(iterator()).limit(resultLength).forEach(gac -> genotypeIndexMapPerGenotypeIndex(gac.index(), gac, oldToNewAlleleIndexMap, result, sortedAlleleCounts));
        return result;
    }

    /**
     * Performs the genotype mapping per new genotype index.
     *
     * @param newGenotypeIndex the target new genotype index.
     * @param alleleCounts tha correspond to {@code newGenotypeIndex}.
     * @param oldToNewAlleleIndexMap the allele mapping.
     * @param destination where to store the new genotype index mapping to old.
     * @param sortedAlleleCountsBuffer a buffer to re-use to get the genotype-allele-count's sorted allele counts.
     */
    private void genotypeIndexMapPerGenotypeIndex(final int newGenotypeIndex, final GenotypeAlleleCounts alleleCounts, final int[] oldToNewAlleleIndexMap, final int[] destination, final int[] sortedAlleleCountsBuffer) {
        final int distinctAlleleCount = alleleCounts.distinctAlleleCount();
        alleleCounts.copyAlleleCounts(sortedAlleleCountsBuffer,0);
        for (int j = 0, jj = 0; j < distinctAlleleCount; j++) {
            final int oldIndex = sortedAlleleCountsBuffer[jj++];
            final int repeats = sortedAlleleCountsBuffer[jj++];
            final int newIndex = oldToNewAlleleIndexMap[oldIndex];
            if (newIndex < 0 || newIndex >= alleleCount) {
                throw new IllegalArgumentException("found invalid new allele index (" + newIndex + ") for old index (" + oldIndex + ")");
            }
            for (int k = 0; k < repeats; k++) {
                alleleHeap.add(newIndex);
            }
        }
        final int genotypeIndex = alleleHeapToIndex(); // this cleans the heap for the next use.
        destination[newGenotypeIndex] = genotypeIndex;
    }

}
