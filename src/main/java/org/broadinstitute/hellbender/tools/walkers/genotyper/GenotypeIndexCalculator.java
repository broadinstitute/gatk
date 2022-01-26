package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.apache.commons.lang3.mutable.MutableInt;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;

public class GenotypeIndexCalculator {
    private final int ploidy;
    private final int alleleCount;

    // instead of re-allocating, we re-use the same container
    private final int[] alleleContainer;



    public GenotypeIndexCalculator(final int ploidy, final int alleleCount) {
        this.ploidy = ploidy;
        this.alleleCount = alleleCount;
        alleleContainer = new int[ploidy];
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
        Utils.validateArg(alleles.length == ploidy,
                () -> String.format("Total allele count %d does not match ploidy %d.", alleles.length, ploidy));

        // Special case ploidy == 0.
        if (ploidy == 0) {
            return 0;
        }

        for (int n = 0; n < ploidy; n++) {
            alleleContainer[n] = alleles[n];
        }
        return calculateIndex(alleleContainer);
    }

    /**
     * Returns the genotype index given the allele counts in format (allele1, count1, allele2, count2. . . )
     *
     * @param alleleCountArray the query allele counts. This must follow the format returned by
     *  {@link GenotypeAlleleCounts#copyAlleleCounts}.
     *
     * @throws IllegalArgumentException if {@code alleleCountArray} is null, has odd length, contains negative counts,
     * or has a total allele count different from the ploidy.

     */
    public int alleleCountsToIndex(final int ... alleleCountArray) {
        Utils.nonNull(alleleCountArray, "the allele counts cannot be null");
        Utils.validateArg((alleleCountArray.length & 1) == 0, "the allele counts array cannot have odd length");
        int n = 0;
        int totalCount = 0;
        for (int i = 0; i < alleleCountArray.length; i += 2) {
            final int allele = alleleCountArray[i];
            final int count = alleleCountArray[i+1];
            totalCount += count;
            Utils.validateArg(count >= 0, "no allele count can be less than 0");
            for (int j = 0; j < count; j++) {
                alleleContainer[n++] = allele;
            }
        }
        final int totalCountLambda = totalCount; // variables must be final in Java lambda
        Utils.validateArg(totalCount == ploidy, () -> String.format("Total allele count %d does not match ploidy %d.", totalCountLambda, ploidy));
        return calculateIndex(alleleContainer);
    }

    /**
     * Calculate the "old" genotype index for the ploidy and allele count of this instance given a GenotypeAlleleCounts
     * object in some new basis of alleles and a int -> int map (in the form of an array) to translate from new allele
     * indices to the "old" allele indices of this instance.
     */
    public int alleleCountsToIndex(final GenotypeAlleleCounts newGAC, final int[] newToOldAlleleMap) {
        final MutableInt n = new MutableInt(0);
        newGAC.forEachAlleleIndexAndCount((newAllele, count) -> {
            final int oldAllele = newToOldAlleleMap[newAllele];
            new IndexRange(0, count).forEach(k -> alleleContainer[n.getAndIncrement()] = oldAllele);
        });

        return calculateIndex(alleleContainer);
    }

    /**
     * Example: suppose our genotype is ABC.  Then the index is the sum of (1) the number of ploidy 3 genotypes before
     * reaching C in the third position, (2) the number of ploidy 2 genotypes before reaching B in the 2nd position, and
     * (3) the number of ploidy 1 genotypes before reaching A in the 1st position.
     */
    private static int calculateIndex(final int[] alleles) {
        // traverse alleles from highest to lowest index
        sort(alleles);
        final int ploidy = alleles.length;
        return new IndexRange(0, ploidy).sumInt(n ->  {
            final int allele = alleles[ploidy - n - 1];
            return (int) GenotypeLikelihoodCalculators.numberOfGenotypesBeforeAllele(ploidy - n, allele);
        });
    }

    private static void sort(final int[] array) {
        final int len = array.length;
        if (len < 2) {
            return;
        } else if (len == 2) {
            if (array[1] < array[0]) {
                final int tmp = array[0];
                array[0] = array[1];
                array[1] = tmp;
            }
        } else {
            Arrays.sort(array);
        }
    }
}
