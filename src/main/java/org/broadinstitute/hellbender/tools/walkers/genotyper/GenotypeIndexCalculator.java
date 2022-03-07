package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;

/**
 * Utilities class for calculations involving the canonical enumeration of (unphased) genotypes.
 *
 * For diploid genotypes with alleles A, B, C. . . this ordering is AA, AB, BB, AC, BC, CC. . .
 *
 * For triploid genotypes it is AAA, AAB, ABB, BBB, AAC, ABC, BBC, ACC, BCC, CCC. . .
 *
 * Note that we may define the ordering recursively.  Letting g = {g_1,g_2,g_3. . .,g_N} and h = {h_1,h_2,h_3,h_N} = genotypes comprising
 * alleles g_1,g_2,g_3 and h_1,h_2,h_3, respectively:
 *    (i)   the order of haploid genotypes is simply the allele ordering
 *    (ii)  if g_N < h_N then g < h
 *    (iii) if g_N = h_N then the order is that of the first N-1 alleles
 *
 * Note also that whenever possible it is best to traverse all genotypes in the canonical order without the random index calculations
 * provided here.  However, when subsetting, reordering, merging, and adding alleles it is necessary to translate indices from
 * one basis of alleles to another.  In such cases efficient index calculations are important.
 */
public class GenotypeIndexCalculator {

    private GenotypeIndexCalculator() {}

    /**
     *     How many genotypes with given ploidy appear in the standard order before a given allele is reached.
     *
     *     For example, considering alleles A, B, C, D, etc ... (indexed 0, 1, 2, ... respectively):
     *     f(3,A) = f(3,0) = 0 as the first genotype AAA contains A.
     *     f(3,B) = f(3,1) = 1 as the second genotype AAB contains B.
     *     f(3,C) = f(3,2) = 4 as the first genotype that contains C, AAC follows: AAA AAB ABB BBB
     *     f(4,D) = f(4,3) = 15 as AAAD follows AAAA AAAB AABB ABBB BBBB AAAC AABC ABBC BBBC AACC ABCC BBCC ACCC BCCC CCCC
     *
     *     There is a simple closed-form expression for this.  Any genotype with ploidy p and a alleles can be encoded
     *     by p 'x's and a - 1 '/'s, where each x represents one allele count and each slash divides between consecutive alleles.
     *     For example, with p = 3 and a = 3 we have xxx// representing AAA, //xxx representing CCC, x/x/x representing ABC,
     *     and xx//x representing AAC.  It is easy to see that any such string corresponds to a genotype, and the number of such
     *     strings is given by the number of places to put the a-1 slashes within the p+a-1 total characters, which is
     *     simply the binomial coefficient (p+a-1)C(a-1).  Considering that allele indices are zero-based, we also have
     *     f(p,a) = (p+a-1)C(a-1).
     *
     *     See discussion at https://genome.sph.umich.edu/wiki/Relationship_between_Ploidy,_Alleles_and_Genotypes
     */
    public static long indexOfFirstGenotypeWithAllele(final int ploidy, final int allele) {
        return allele == 0 ? 0 : CombinatoricsUtils.binomialCoefficient(ploidy + allele - 1, allele - 1);
    }

    /**
     * Returns the number of possible genotypes given the ploidy and number of different alleles.
     * @param ploidy the requested ploidy.
     * @param alleleCount the requested number of alleles.
     *
     * @throws IllegalArgumentException if {@code ploidy} or {@code alleleCount} is negative or
     *                                      the number of genotypes is too large (more than {@link Integer#MAX_VALUE}).
     *
     * @return the number of genotypes given ploidy and allele count (0 or greater).
     */
    public static int genotypeCount(final int ploidy, final int alleleCount) {
        final long result = indexOfFirstGenotypeWithAllele(ploidy, alleleCount);
        Utils.validateArg(result != MathUtils.LONG_OVERFLOW && result < Integer.MAX_VALUE, () ->
                String.format("the number of genotypes is too large for ploidy %d and %d alleles: approx. %.0f", ploidy, alleleCount,
                        CombinatoricsUtils.binomialCoefficientDouble(ploidy + alleleCount - 1, alleleCount - 1)));
        return (int) result;
    }

    /**
     * Give a list of alleles, returns the likelihood array index.
     *
     * @param alleles the indices of the alleles in the genotype, there should be as many repetition of an
     *                      index as copies of that allele in the genotype. Allele indices do not need to be sorted in
     *                      any particular way.  For example, {A,A,B}, {A,B,A}, {B,A,A} are all valid inputs.
     *
     * @return never {@code null}.
     */
    public static int allelesToIndex(final int... alleles) {
        final int ploidy = alleles.length;
        return ploidy == 0 ? 0 : calculateIndex(Arrays.copyOf(alleles, ploidy));
    }

    /**
     * Returns the genotype index given the allele counts in format (allele1, count1, allele2, count2. . . )
     *
     * @param alleleCountArray the query allele counts.
     *
     * @throws IllegalArgumentException if {@code alleleCountArray} is null, has odd length, contains negative counts,
     * or has a total allele count different from the ploidy.
     */
    public static int alleleCountsToIndex(final int ... alleleCountArray) {
        Utils.nonNull(alleleCountArray, "the allele counts cannot be null");
        Utils.validateArg((alleleCountArray.length & 1) == 0, "the allele counts array cannot have odd length");
        int ploidy = 0;
        for (int i = 0; i < alleleCountArray.length; i += 2) {
            ploidy += alleleCountArray[i+1];
        }
        final int[] alleleContainer = new int[ploidy];


        int n = 0;
        for (int i = 0; i < alleleCountArray.length; i += 2) {
            final int allele = alleleCountArray[i];
            final int count = alleleCountArray[i+1];
            Utils.validateArg(count >= 0, "no allele count can be less than 0");
            for (int j = 0; j < count; j++, n++) {
                alleleContainer[n] = allele;
            }
        }
        return calculateIndex(alleleContainer);
    }

    /**
     * Calculate the "old" genotype index for the ploidy and allele count of this instance given a GenotypeAlleleCounts
     * object in some new basis of alleles and a int -> int map (in the form of an array) to translate from new allele
     * indices to the "old" allele indices of this instance.
     */
    public static int alleleCountsToIndex(final GenotypeAlleleCounts newGAC, final int[] newToOldAlleleMap) {
        final int[] alleleContainer = new int[newGAC.ploidy()];
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
        final int ploidy = alleles.length;

        // traverse alleles from highest to lowest index
        Arrays.sort(alleles);
        return new IndexRange(0, ploidy).sumInt(n ->  {
            final int allele = alleles[ploidy - n - 1];
            return (int) indexOfFirstGenotypeWithAllele(ploidy - n, allele);
        });
    }

    /**
     * Compute the maximally acceptable allele count (ref allele included) given the maximally acceptable genotype count.
     * @param ploidy            sample ploidy
     * @param maxGenotypeCount  maximum number of genotype count used to calculate upper bound on number of alleles given ploidy
     * @throws IllegalArgumentException if {@code ploidy} or {@code alleleCount} is negative.
     * @return                  the maximally acceptable allele count given ploidy and maximum number of genotypes acceptable
     */
    public static int computeMaxAcceptableAlleleCount(final int ploidy, final int maxGenotypeCount){
        Utils.validateArg(ploidy >= 0, () -> "negative ploidy " + ploidy);

        if (ploidy == 1) {
            return maxGenotypeCount;
        }
        final double logMaxGenotypeCount = FastMath.log(maxGenotypeCount);

        // Math explanation: genotype count is determined by ${P+A-1 \choose A-1}$, this leads to constraint
        // $\log(\frac{(P+A-1)!}{(A-1)!}) \le \log(P!G)$,
        // where $P$ is ploidy, $A$ is allele count, and $G$ is maxGenotypeCount
        // The upper and lower bounds of the left hand side of the constraint are $P \log(A-1+P)$ and $P \log(A)$
        // which require $A$ to be searched in interval $[exp{\log(P!G)/P} - (P-1), exp{\log(P!G)/P}]$
        // Denote $[10^{\log(P!G)/P}$ as $x$ in the code.

        final double x = FastMath.exp((CombinatoricsUtils.factorialLog(ploidy) + logMaxGenotypeCount)/ploidy );
        final int lower = (int)Math.floor(x) - ploidy - 1;
        final int upper = (int)Math.ceil(x);
        for(int a=upper; a>=lower; --a){// check one by one

            final double logGTCnt = CombinatoricsUtils.binomialCoefficientLog(ploidy+a-1, a-1);
            if(logMaxGenotypeCount >= logGTCnt) {
                return a;
            }
        }
        throw new GATKException("Code should never reach here.");
    }

    /**
     * Composes a genotype index map given a allele index recoding such that result[i] is the index of the old
     * genotype corresponding to the ith new genotype.
     *
     * @param newToOldAlleleMap allele recoding such that newToOldAlleleMap[i] is the index of the old allele
     *                               corresponding to the ith new allele
     *
     * @return never {@code null}.
     */
    public static int[] newToOldGenotypeMap(final int ploidy, final int[] newToOldAlleleMap) {
        Utils.nonNull(newToOldAlleleMap);
        final int newAlleleCount = newToOldAlleleMap.length;

        final int[] result = new int[genotypeCount(ploidy, newAlleleCount)];
        for (final GenotypeAlleleCounts newGAC : GenotypeAlleleCounts.iterable(ploidy, newAlleleCount)) {
            result[newGAC.index()] = alleleCountsToIndex(newGAC, newToOldAlleleMap);
        }

        return result;
    }
}
