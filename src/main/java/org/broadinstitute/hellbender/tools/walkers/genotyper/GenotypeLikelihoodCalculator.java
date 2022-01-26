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
     * Precomputed values such that value[genotype g][read] is the log likelihood log_10 lk(g | r)
     *  = log_10_sum_log_10(log_10(f1) + log_10 Lk ( read | a1 ), . . .)
     *  = log_10(f1*lk(r|a1) + f2*lk(r|a2). . .)
     *  where the genotype g has alleles a1, a2. . . with frequencies f1, f2. . .
     */
    final double[][] readLikelihoodsByGenotypeIndex;


    /**
     * Cache of the last genotype-allele-count requested using {@link #genotypeAlleleCountsAt(int)}, when it
     * goes beyond the maximum genotype-allele-count static capacity. Check on that method documentation for details.
     */
    private GenotypeAlleleCounts lastOverheadCounts;

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
     * @return the raw array (in log10 likelihoods space) of the GL for each genotype
     */
    <EVIDENCE, A extends Allele> double[] getReadRawReadLikelihoodsByGenotypeIndex(final LikelihoodMatrix<EVIDENCE, A> likelihoods) {
        Utils.nonNull(likelihoods);
        Utils.validateArg(likelihoods.numberOfAlleles() == alleleCount, "mismatch between allele list and alleleCount");
        final int readCount = likelihoods.evidenceCount();


        // the result satisfies result[g] = sum_reads ( log_10_sum_log_10(log_10(f1) + log_10 Lk ( read | a1 ). . .) - num_reads*log_10 ploidy
        // where the genotype g has alleles a1, a2. . . with frequencies f1, f2. . .
        // please note that the above expression is the entirety of this class, obscured by a whole lot of optimization
        // tht the observant reader may notice is premature
        return genotypeLikelihoods(_____, readCount);
    }

    /**
     * Calculates the final genotype likelihood array out of the likelihoods for each genotype per read.
     *
     * @param readLikelihoodsByGenotypeIndex likelihoods log_10 lk(g | r) for genotypes g and reads r
     *                                       The values therein are entry[g][r] = log_10_sum_log_10(log_10(f1) + log_10 Lk ( read | a1 ), . . .)
     *                                       = log_10(f1*lk(r|a1) + f2*lk(r|a2). . .)
     *                                       where the genotype g has alleles a1, a2. . . with frequencies f1, f2. . .
     * @param readCount number of reads in the input likelihood arrays in {@code genotypeLikelihoodByRead}.
     * @return array result such that result[g] = sum_reads(log_10 lk(read | g)) - num_reads*log_10 ploidy
     */
    double[] genotypeLikelihoods(final double[][] readLikelihoodsByGenotypeIndex, final int readCount) {

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
