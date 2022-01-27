package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.apache.commons.lang3.mutable.MutableDouble;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;

import java.util.Arrays;
import java.util.Iterator;

public class GenotypeLikelihoodCalculator implements Iterable<GenotypeAlleleCounts> {
    /**
     * Genotype table for this calculator.
     *
     * <p>It is ensure that it contains all the genotypes for this calculator ploidy and allele count, maybe more. For
     * that reason you must use {@link #genotypeCount} when iterating through this array and not relay on its length.</p>
     */
    private final GenotypeAlleleCounts[] genotypeAlleleCounts;

    final int genotypeCount;

    final int alleleCount;

    final int ploidy;

    final GenotypeIndexCalculator genotypeIndexCalculator;
    
    /**
     * Cache of the last genotype-allele-count requested using {@link #genotypeAlleleCountsAt(int)}, when it
     * goes beyond the maximum genotype-allele-count static capacity. Check on that method documentation for details.
     */
    private GenotypeAlleleCounts lastOverheadCounts;

    private static final int INITIAL_READ_CAPACITY = 10;

    /**
     * How many reads the calculator supports.
     *
     * This figure is increased dynamically by {@code ensureReadCapacity}.
     */
    private int readCapacity = INITIAL_READ_CAPACITY;

    /**
     * Buffer field used as a temporary container when one value per read is stored.
     *
     * In the multiallelic calculation we accumulate the likelihood contribution of each read one allele at a time.  That is,
     * for genotype containing alleles A, B, C, we first fill the buffer with the likelihood contributions from allele A, then
     * we make a second pass and add the contributions from allele B, then allele C.  Traversing all the reads in each
     * allele row of the likelihoods array in this manner is cache-friendly and makes an enormous difference in runtime.
     */
    private double[] perReadBuffer = new double[INITIAL_READ_CAPACITY];



    public GenotypeLikelihoodCalculator(final int ploidy, final int alleleCount,
                                        final GenotypeAlleleCounts[][] genotypeTableByPloidy) {
        genotypeAlleleCounts = genotypeTableByPloidy[ploidy];
        genotypeCount = (int) GenotypeLikelihoodCalculators.numberOfGenotpyes(ploidy, alleleCount);
        this.alleleCount = alleleCount;
        this.ploidy = ploidy;
        genotypeIndexCalculator = new GenotypeIndexCalculator(ploidy, alleleCount);
    }

    /**
     * Give a list of alleles, returns the likelihood array index.
     * @param alleles the indices of the alleles in the genotype, there should be as many repetition of an
     *                      index as copies of that allele in the genotype. Allele indices do not need to be sorted in
     *                      any particular way.
     *
     * @return never {@code null}.
     */
    public int allelesToIndex(final int... alleles) {
        return genotypeIndexCalculator.allelesToIndex(alleles);
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
     * <p>If {@code index} is larger than {@link GenotypeLikelihoodCalculators#MAXIMUM_CACHED_GENOTYPES_PER_CALCULATOR},
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
        if (index < GenotypeLikelihoodCalculators.MAXIMUM_CACHED_GENOTYPES_PER_CALCULATOR) {
            return genotypeAlleleCounts[index];
        } else if (lastOverheadCounts == null || lastOverheadCounts.index() > index) {
            final GenotypeAlleleCounts result = genotypeAlleleCounts[GenotypeLikelihoodCalculators.MAXIMUM_CACHED_GENOTYPES_PER_CALCULATOR - 1].copy();
            result.increase(index - GenotypeLikelihoodCalculators.MAXIMUM_CACHED_GENOTYPES_PER_CALCULATOR + 1);
            lastOverheadCounts = result;
            return result.copy();
        } else {
            lastOverheadCounts.increase(index - lastOverheadCounts.index());
            return lastOverheadCounts.copy();
        }
    }

    /**
     * Calculate the likelihoods given the list of alleles and the likelihood map.
     *
     * @param likelihoods the likelihood matrix all alleles vs all reads.
     *
     * @throws IllegalArgumentException if {@code alleleList} is {@code null} or {@code likelihoods} is {@code null}
     *     or the alleleList size does not match the allele-count of this calculator, or there are missing allele vs
     *     read combinations in {@code likelihoods}.
     *
     * @return never {@code null}.
     */
    public <EVIDENCE, A extends Allele> GenotypeLikelihoods genotypeLikelihoods(final LikelihoodMatrix<EVIDENCE, A> likelihoods) {
        final double[] readLikelihoodsByGenotypeIndex = getReadRawReadLikelihoodsByGenotypeIndex(likelihoods);
        return GenotypeLikelihoods.fromLog10Likelihoods(readLikelihoodsByGenotypeIndex);
    }

    /**
     * A helper method that actually does the matrix operations but returns the raw values.
     *
     * @param likelihoods   log 10 likelihood matrix indexed by allele, then read
     * @return the raw array (in log10 likelihoods space) of the GL for each genotype
     */
    <EVIDENCE, A extends Allele> double[] getReadRawReadLikelihoodsByGenotypeIndex(final LikelihoodMatrix<EVIDENCE, A> likelihoods) {
        Utils.nonNull(likelihoods);
        Utils.validateArg(likelihoods.numberOfAlleles() == alleleCount, "mismatch between allele list and alleleCount");
        final int readCount = likelihoods.evidenceCount();
        ensureReadCapacity(readCount);

        final double[][] likelihoodsByAlleleAndRead = likelihoods.asRealMatrix().getData();

        final boolean triallelicGenotypesPossible = alleleCount > 2 && ploidy > 2;

        // non-log space likelihoods for multiallelic computation
        // TODO: ditto
        final double[][] nonLogLikelihoodsByAlleleAndRead = triallelicGenotypesPossible ? likelihoods.asRealMatrix().getData() : null;
        final MutableDouble multiallelicNormalization = new MutableDouble(0);
        if (triallelicGenotypesPossible) {
            makeNormalizedNonLogLikelihoods(readCount, nonLogLikelihoodsByAlleleAndRead, multiallelicNormalization);
        }

        final double[] result = new double[genotypeCount];

        Utils.stream(iterator()).forEach(alleleCounts -> {
            final int componentCount = alleleCounts.distinctAlleleCount();
            final int genotypeIndex = alleleCounts.index();
            if (componentCount == 1) {
                // homozygous case: log P(reads|AAAAA. . .) = sum_{reads} log P(read|A)
                final int allele = alleleCounts.alleleIndexAt(0);
                result[genotypeIndex] = MathUtils.sum(likelihoodsByAlleleAndRead[allele]);
            } else if (componentCount == 2) {
                // biallelic het case: log P(reads | nA copies of A, nB copies of B) = sum_{reads} log[(nA * P(read | A) + nB * P(read | B))] -log(ploidy)
                final double[] logLks1 = likelihoodsByAlleleAndRead[alleleCounts.alleleIndexAt(0)];
                final int freq1 = alleleCounts.alleleCountAt(0);
                final double log10Freq1 = MathUtils.log10(freq1);
                final double[] logLks2  = likelihoodsByAlleleAndRead[alleleCounts.alleleIndexAt(1)];
                final double log10Freq2 = MathUtils.log10(ploidy - freq1);

                result[genotypeIndex] = new IndexRange(0, readCount).sum(r -> MathUtils.approximateLog10SumLog10(logLks1[r] + log10Freq1, logLks2[r] + log10Freq2))
                        - readCount * MathUtils.log10(ploidy);
            } else {
                // the multiallelic case is conceptually the same as the biallelic case but done in non-log space
                // We implement in a cache-friendly way by adding nA * P(read|A) to the per-read buffer for all reads (inner loop), and all alleles (outer loop)
                Arrays.fill(perReadBuffer,0, readCount, 0);
                alleleCounts.forEachAlleleIndexAndCount((a, f) -> new IndexRange(0, readCount).forEach(r -> perReadBuffer[r] += f * nonLogLikelihoodsByAlleleAndRead[a][r]));
                result[genotypeIndex] = new IndexRange(0, readCount).sum(r -> FastMath.log10(perReadBuffer[r])) - readCount * MathUtils.log10(ploidy) + multiallelicNormalization.doubleValue();
            }
        });
        return result;
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

    private void ensureReadCapacity(final int requestedCapacity) {
        if (readCapacity < requestedCapacity) {
            readCapacity = 2 * requestedCapacity;
            perReadBuffer = new double[readCapacity];
        }
    }


    // note that if the input has a high index that is not cached, it will be mutated in order to form the output
    private GenotypeAlleleCounts nextGenotypeAlleleCounts(final GenotypeAlleleCounts alleleCounts) {
        final int index = alleleCounts.index();
        if (index < (GenotypeLikelihoodCalculators.MAXIMUM_CACHED_GENOTYPES_PER_CALCULATOR - 1)) {
            return genotypeAlleleCounts[index + 1];
        } else if (index == GenotypeLikelihoodCalculators.MAXIMUM_CACHED_GENOTYPES_PER_CALCULATOR - 1) {
            return genotypeAlleleCounts[index].copy().increase();
        } else {
            return alleleCounts.increase();
        }
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
     * Returns the likelihood index given the allele counts in format (allele1, count1, allele2, count2. . . )
     *
     * @param alleleCountArray allele counts in the format returned by {@link GenotypeAlleleCounts#copyAlleleCounts}.
     */
    public int alleleCountsToIndex(final int ... alleleCountArray) {
        return genotypeIndexCalculator.alleleCountsToIndex(alleleCountArray);
    }

    /**
     * Composes a genotype index map given a allele index recoding such that result[i] is the index of the old
     * genotype corresponding to the ith new genotype.
     *
     * @param newToOldAlleleMap allele recoding such that newToOldAlleleMap[i] is the index of the old allele
     *                               corresponding to the ith new allele
     *
     * @throws IllegalArgumentException if this calculator cannot handle the recoding provided. This is
     * the case when either {@code newToOldAlleleMap}'s length or any of its element (+ 1 as they are 0-based) is larger
     * this calculator's {@link #alleleCount()}. Also if any {@code oldToNewAllelesIndexMap} element is negative.
     *
     * @return never {@code null}.
     */
    public int[] newToOldGenotypeMap(final int[] newToOldAlleleMap) {
        Utils.nonNull(newToOldAlleleMap);
        final int newAlleleCount = newToOldAlleleMap.length;
        Utils.validateArg(newAlleleCount <= alleleCount,
                () -> String.format("New allele count %d exceeds old allele count %d.", newAlleleCount, alleleCount));
        final int newGenotypeCount = newAlleleCount == alleleCount ? genotypeCount :
                GenotypeLikelihoodCalculators.genotypeCount(ploidy, newAlleleCount);

        final int[] result = new int[newGenotypeCount];
        GenotypeAlleleCounts newGAC = genotypeAlleleCounts[0];
        for (int i = 0; i < newGenotypeCount; i++) {
            result[i] = genotypeIndexCalculator.alleleCountsToIndex(newGAC, newToOldAlleleMap);

            if (i < newGenotypeCount - 1) {
                newGAC = nextGenotypeAlleleCounts(newGAC);
            }
        }
        return result;
    }

    @Override
    public Iterator<GenotypeAlleleCounts> iterator() {
        return new Iterator<GenotypeAlleleCounts>() {
            private int index = 0;
            private GenotypeAlleleCounts alleleCounts = genotypeAlleleCounts[0];

            @Override
            public boolean hasNext() {
                return index < genotypeCount;
            }

            @Override
            public GenotypeAlleleCounts next() {
                alleleCounts = index++ == 0 ? genotypeAlleleCounts[0] : nextGenotypeAlleleCounts(alleleCounts);
                return alleleCounts;
            }
        };
    }
}
