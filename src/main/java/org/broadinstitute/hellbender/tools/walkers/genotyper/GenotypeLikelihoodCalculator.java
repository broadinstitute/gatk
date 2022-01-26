package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.apache.commons.lang3.mutable.MutableInt;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;

import java.util.Iterator;

public class GenotypeLikelihoodCalculator implements Iterable<GenotypeAlleleCounts> {
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
    final int genotypeCount;
    /**
     * Number of genotyping alleles for this calculator.
     */
    final int alleleCount;
    /**
     * Ploidy for this calculator.
     */
    final int ploidy;

    final GenotypeIndexCalculator genotypeIndexCalculator;


    /**
     * Buffer used as a temporary container for likelihood components for genotypes stratified by reads.
     *
     * <p>
     *     It is indexed by genotype index and then by read index. The read capacity is increased as needed by calling
     *     {@link #ensureReadCapacity(int) ensureReadCapacity}.
     * </p>
     */
    final double[][] readLikelihoodsByGenotypeIndex;


    /**
     * Cache of the last genotype-allele-count requested using {@link #genotypeAlleleCountsAt(int)}, when it
     * goes beyond the maximum genotype-allele-count static capacity. Check on that method documentation for details.
     */
    private GenotypeAlleleCounts lastOverheadCounts;
    /**
     * Precomputed values of log_10(frequency) + log_10 Lk(read | allele)
     */
    double[] readAlleleLikelihoodByAlleleCount = null;

    /**
     * Indicates how many reads the calculator supports.
     *
     * <p>This figure is increased dynamically as per the
     * calculation request calling {@link #ensureReadCapacity(int) ensureReadCapacity}.<p/>
     */
    private int readCapacity = -1;


    public GenotypeLikelihoodCalculator(final int ploidy, final int alleleCount,
                                        final GenotypeAlleleCounts[][] genotypeTableByPloidy) {
        genotypeAlleleCounts = genotypeTableByPloidy[ploidy];
        genotypeCount = (int) GenotypeLikelihoodCalculators.numberOfGenotpyes(ploidy, alleleCount);
        this.alleleCount = alleleCount;
        this.ploidy = ploidy;
        genotypeIndexCalculator = new GenotypeIndexCalculator(ploidy, alleleCount);
        readLikelihoodsByGenotypeIndex = new double[genotypeCount][];
    }

    /**
     * Makes sure that temporal arrays and matrices are prepared for a number of reads to process.
     * @param requestedCapacity number of read that need to be processed.
     */
    public void ensureReadCapacity(final int requestedCapacity) {
        Utils.validateArg(requestedCapacity >= 0, "capacity may not be negative");
        if (readCapacity == -1) { // first time call.
            final int minimumCapacity = Math.max(requestedCapacity, 10); // Never go too small, 10 is the minimum.
            for (int i = 0; i < genotypeCount; i++) {
                readLikelihoodsByGenotypeIndex[i] = new double[minimumCapacity];
            }
            readCapacity = minimumCapacity;
        } else if (readCapacity < requestedCapacity) {
            final int doubleCapacity = (requestedCapacity << 1);
            for (int i = 0; i < genotypeCount; i++) {
                readLikelihoodsByGenotypeIndex[i] = new double[doubleCapacity];
            }
            readCapacity = doubleCapacity;
        }
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
     * @return the raw array (in log10 likelihoods space) of the GL for each genotype
     */
    <EVIDENCE, A extends Allele> double[] getReadRawReadLikelihoodsByGenotypeIndex(final LikelihoodMatrix<EVIDENCE, A> likelihoods) {
        Utils.nonNull(likelihoods);
        Utils.validateArg(likelihoods.numberOfAlleles() == alleleCount, "mismatch between allele list and alleleCount");
        final int readCount = likelihoods.evidenceCount();
        ensureReadCapacity(readCount);

        /**
         * Precompute a 1D flattened matrix corresponding to the 3D matrix:
         *
         * result[allele][frequency][read] = log_10(frequency) + log_10 Lk ( read_x | allele_y )
         *
         * That's right -- all this machinery just to precompute. . . addition.  If that strikes you not only as
         * premature optimization but as bald non-optimization, you're correct!
         */
        final double[] readLikelihoodComponentsByAlleleCount = null;    // null is a placeholder for later deletion

        final double[][] genotypeLikelihoodByRead = genotypeLikelihoodByRead(readLikelihoodComponentsByAlleleCount,readCount);
        return genotypeLikelihoods(genotypeLikelihoodByRead, readCount);
    }

    /**
     * Calculates the final genotype likelihood array out of the likelihoods for each genotype per read.
     *
     * @param readLikelihoodsByGenotypeIndex <i>[g][r]</i> likelihoods for each genotype <i>g</i> and <i>r</i>.
     * @param readCount number of reads in the input likelihood arrays in {@code genotypeLikelihoodByRead}.
     * @return never {@code null}, one position per genotype where the <i>i</i> entry is the likelihood of the ith
     *   genotype (0-based).
     */
    double[] genotypeLikelihoods(final double[][] readLikelihoodsByGenotypeIndex, final int readCount) {
        final double[] result = new double[genotypeCount];
        final double denominator = readCount * MathUtils.log10(ploidy);
        // instead of dividing each read likelihood by ploidy ( so subtract log10(ploidy) )
         // we multiply them all and the divide by ploidy^readCount (so substract readCount * log10(ploidy) )
        for (int g = 0; g < genotypeCount; g++) {
            result[g] = MathUtils.sum(readLikelihoodsByGenotypeIndex[g], 0, readCount) - denominator;
        }
        return result;
    }

    /**
     * Calculates the likelihood component of each read on each genotype.
     *
     * NOTE: this is not actually the read likelihood component for each genotype, it is the sum of the log read likelihoods components
     *       for each genotype without having been normalized by the the denominator of the ploidy, that happens in the final step
     *
     * @param readLikelihoodComponentsByAlleleCount precomputed matrix M such that
     *                                              result[allele][frequency][read] = log_10(frequency) + log_10 Lk ( read_x | allele_y )
     * @param readCount number of reads in {@code readLikelihoodComponentsByAlleleCount}.
     * @return never {@code null}.
     */
    protected double[][] genotypeLikelihoodByRead(final double[] readLikelihoodComponentsByAlleleCount, final int readCount) {

        // Here we don't use the convenience of {@link #genotypeAlleleCountsAt(int)} within the loop to spare instantiations of
        // GenotypeAlleleCounts class when we are dealing with many genotypes.
        GenotypeAlleleCounts alleleCounts = genotypeAlleleCounts[0];

        for (int genotypeIndex = 0; genotypeIndex < genotypeCount; genotypeIndex++) {
            final double[] readLikelihoods = this.readLikelihoodsByGenotypeIndex[genotypeIndex];
            final int componentCount = alleleCounts.distinctAlleleCount();
            switch (componentCount) {
                case 1: // for homozygous genotype with one allele with frequency = ploidy
                    // fill readLikelihoods[read] = log_10(ploidy) + log_10 Lk ( read | allele )
                case 2: // for biallelic genotype with alleles a1, a2 with frequencies f1, f2
                    // fill readLikelihoods[read] = log_10_sum_log_10(log_10(f1) + log_10 Lk ( read | a1 ), log_10(f2) + log_10 Lk ( read | a2 ))
                default: // general case with alleles a1, a2. . . with frequencies f1, f2. . .
                    // fill readLikelihoods[read] = log_10_sum_log_10(log_10(f1) + log_10 Lk ( read | a1 ), . . .)
            }
            if (genotypeIndex < genotypeCount - 1) {
                alleleCounts = nextGenotypeAlleleCounts(alleleCounts);
            }
        }
        return readLikelihoodsByGenotypeIndex;
    }

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
     * @param alleleCountArray the query allele counts. This must follow the format returned by
     *  {@link GenotypeAlleleCounts#copyAlleleCounts}.
     *
     * @throws IllegalArgumentException if {@code alleleCountArray} is null, has odd length, contains negative counts,
     * or has a total allele count different from the ploidy.

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
            private int genotypeIndex = 0;
            private final GenotypeAlleleCounts increasingAlleleCounts = genotypeCount < GenotypeLikelihoodCalculators.MAXIMUM_CACHED_GENOTYPES_PER_CALCULATOR ?
                    null : genotypeAlleleCounts[GenotypeLikelihoodCalculators.MAXIMUM_CACHED_GENOTYPES_PER_CALCULATOR - 1].copy();

            @Override
            public boolean hasNext() {
                return genotypeIndex < genotypeCount;
            }

            @Override
            public GenotypeAlleleCounts next() {
                if (genotypeIndex < GenotypeLikelihoodCalculators.MAXIMUM_CACHED_GENOTYPES_PER_CALCULATOR) {
                    return genotypeAlleleCounts[genotypeIndex++];
                } else {
                    genotypeIndex++;
                    return increasingAlleleCounts.increase();
                }
            }
        };
    }
}
