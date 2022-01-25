package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Genotype likelihood calculator utility. This class is thread-safe since access to shared mutable state is
 * synchronized.
 *
 * <p>
 *     This class provide genotype likelihood calculators with any number of alleles able given an arbitrary ploidy and allele
 *     count (number of distinct alleles).
 * </p>
 */
public final class GenotypeLikelihoodCalculators {

    private static final Logger logger = LogManager.getLogger(GenotypeLikelihoodCalculators.class);

    /**
     * The current maximum ploidy supported by the tables.
     * <p>
     *     Its initial value indicates the initial capacity of the shared {@link #genotypeTableByPloidy}. Feel free
     *     to change it to anything reasonable that is non-negative.
     * </p>
     */
    private int maximumPloidy = 2; // its initial value is the initial capacity of the shared tables.

    /**
     * Maximum possible number of {@link GenotypeAlleleCounts} a {@link GenotypeLikelihoodCalculator} of fixed ploidy
     *  and allele count may cache in an array, allowing for fast random access of genotypes.
     *
     *  Iterating over all genotypes via the {@code increase} method is always fast, but accessing the allele counts
     *  for a particular genotype index, eg the one with the greatest likelihood, is slow without a cache.
     */
    public static final int MAXIMUM_CACHED_GENOTYPES_PER_CALCULATOR = 1000;

    /**
     * Mark to indicate genotype-count overflow due to a large number of allele and ploidy;
     */
    static final int GENOTYPE_COUNT_OVERFLOW = -1;

    /**
     * The current maximum allele index supported by the tables.
     * <p>
     *     Feel free to change this initial value to anything reasonable that is non-negative.
     * </p>
     */
    private int maximumAllele = 1; // its initial value is the initial capacity of the shared tables.


    /**
     * Shared table of genotypes given the ploidy sorted by their index in the likelihood array.
     *
     * <p>
     *  Its format is described in {@link #buildGenotypeAlleleCountsTable(int, int)}.
     * </p>
     */
    private GenotypeAlleleCounts[][] genotypeTableByPloidy =
            buildGenotypeAlleleCountsTable(maximumPloidy,maximumAllele);

    public GenotypeLikelihoodCalculators(){

    }

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
     */
    public static long numberOfGenotypesBeforeAllele(final int ploidy, final int allele) {
        return allele == 0 ? 0 : MathUtils.exactBinomialCoefficient(ploidy + allele - 1, allele - 1);
    }

    /**
     * Number of genotypes for a given ploidy and allele count.  Follows the same logic as
     * {@link GenotypeLikelihoodCalculator::numberOfGenotypesBeforeAllele}.  That is, the number of genotypes of ploidy
     * p and a alleles is equal to the number of genotypes that don't contain an imaginary a-th zero-based allele.
     *
     * We could inline this and just call the other function, but having a distinct name is much clearer.
     */
    public static long numberOfGenotpyes(final int ploidy, final int alleleCount) {
        return numberOfGenotypesBeforeAllele(ploidy, alleleCount);
    }

    /**
     * Composes a table with the lists of all possible genotype allele counts given the the ploidy and maximum allele index.
     * <p>
     *     The resulting matrix has at least as many rows as {@code maximumPloidy } + 1 as the first row with index 0 correspond
     *     to ploidy == 0. Each row array has as many positions as necessary to contain all possible genotype-allele-counts in increasing order.
     *     This quantity varies with the ploidy.
     * </p>
     *
     * <p>
     *     Therefore <code>result[3][4]</code> would contain the 5th genotype with ploidy 3, and <code>result[4].length</code>
     *     would be equal to the count of possible genotypes for ploidy 4.
     * </p>
     *
     * @param maximumPloidy maximum ploidy to use in queries to the resulting table.
     * @param maximumAllele maximum allele index to use in queries to the resulting table.
     *
     * @throws IllegalArgumentException if {@code maximumPloidy} or {@code maximumAllele} is negative, or {@code offsetTable} is {@code null},
     *   or it does not have the capacity to handle the requested maximum ploidy or allele index.
     *
     * @return never {@code null}.
     */
    private static GenotypeAlleleCounts[][] buildGenotypeAlleleCountsTable(final int maximumPloidy, final int maximumAllele) {
        checkPloidyAndMaximumAllele(maximumPloidy, maximumAllele);
        final int rowCount = maximumPloidy + 1;
        final GenotypeAlleleCounts[][] result = new GenotypeAlleleCounts[rowCount][]; // each row has a different number of columns.

        for (int ploidy = 0; ploidy <= maximumPloidy; ploidy++) {
            result[ploidy] = buildGenotypeAlleleCountsArray(ploidy, maximumAllele);
        }

        return result;
    }

    /**
     * Builds a genotype-allele-counts array given the genotype ploidy and how many genotype you need.
     * <p>
     *     The result is guarantee to have exactly {@code length} positions and the elements are sorted
     *     in agreement with the standard way to display genotypes following the VCF standard.
     * </p>
     *
     * <p> Notice that is possible to request ploidy ==0. In that case the resulting array will have repetitions
     * of the empty genotype allele count.
     * </p>
     *
     * <p>
     *     For example,
     *
     *     <pre>
     *         ploidy = 1, length = 5 : [ {A}, {B}, {C}, {D}, {E} ]
     *         ploidy = 2, length = 7 : [ {AA}, {AB}, {BB}, {AC}, {BC}, {CC}, {AD}
     *         ploidy = 3, length = 10 : [ {AAA}, {AAB}, {ABB}, {BBB}, {AAC}, {ABC}, {BBC}, {BCC}, {CCC}, {AAD} ]
     *     </pre>
     * </p>
     *
     * @param ploidy requested ploidy.
     * @param alleleCount number of different alleles that the genotype table must support.
     * @throws IllegalArgumentException if {@code ploidy} or {@code length} is negative.
     *
     * @return never {@code null}, follows the specification above.
     */
    private static GenotypeAlleleCounts[] buildGenotypeAlleleCountsArray(final int ploidy, final int alleleCount) {
        Utils.validateArg(ploidy >= 0, () -> "the requested ploidy cannot be negative: " + ploidy);
        Utils.validateArg(alleleCount >= 0, () -> "the requested maximum allele cannot be negative: " + alleleCount);
        final long length = numberOfGenotpyes(ploidy, alleleCount);
        final int strongRefLength = length == GENOTYPE_COUNT_OVERFLOW ? MAXIMUM_CACHED_GENOTYPES_PER_CALCULATOR : (int) Math.min(length, MAXIMUM_CACHED_GENOTYPES_PER_CALCULATOR);
        final GenotypeAlleleCounts[] result = new GenotypeAlleleCounts[strongRefLength];
        result[0] = GenotypeAlleleCounts.first(ploidy);
        for (int genotypeIndex = 1; genotypeIndex < strongRefLength; genotypeIndex++) {
            result[genotypeIndex] = result[genotypeIndex - 1].next();
        }
        return result;
    }


    /**
     * Returns an instance given its ploidy and the number of alleles.
     *
     * @param alleleCount the required allele-count.
     * @param ploidy the required ploidy-count.
     *
     * @throws IllegalArgumentException if either {@code ploidy} or {@code alleleCount} is negative, or the resulting number of genotypes is too large.
     *
     * @return never {@code null}.
     */
    public synchronized GenotypeLikelihoodCalculator getInstance(final int ploidy, final int alleleCount) {
        calculateGenotypeCountsUsingTablesAndValidate(ploidy, alleleCount);

        // At this point the tables must have at least the requested capacity, likely to be much more.
        return new GenotypeLikelihoodCalculator(ploidy, alleleCount, genotypeTableByPloidy);
    }

    /**
     * Calculate genotype counts using the tables and validate that there is no overflow
     */
    private synchronized void calculateGenotypeCountsUsingTablesAndValidate(final int ploidy, final int alleleCount) {
        checkPloidyAndMaximumAllele(ploidy, alleleCount);

        if (calculateGenotypeCountUsingTables(ploidy, alleleCount) == GENOTYPE_COUNT_OVERFLOW) {
            final double largeGenotypeCount = Math.pow(10, MathUtils.log10BinomialCoefficient(ploidy + alleleCount - 1, alleleCount - 1));
            throw new IllegalArgumentException(String.format("the number of genotypes is too large for ploidy %d and allele %d: approx. %.0f", ploidy, alleleCount, largeGenotypeCount));
        }
    }

    /**
     * Returns an instance of the DRAGEN genotypeLikelihoodCalculator given its ploidy and the number of alleles.
     *
     * @param alleleCount the required allele-count.
     * @param ploidy the required ploidy-count.
     *
     * @throws IllegalArgumentException if either {@code ploidy} or {@code alleleCount} is negative, or the resulting number of genotypes is too large.
     *
     * @return never {@code null}.
     */
    public synchronized GenotypeLikelihoodCalculatorDRAGEN getInstanceDRAGEN(final int ploidy, final int alleleCount) {
        Utils.validate(ploidy == 2, "DRAGEN genotyping mode currently only supports diploid samples");
        calculateGenotypeCountsUsingTablesAndValidate(ploidy, alleleCount);

        // At this point the tables must have at least the requested capacity, likely to be much more.
        return new GenotypeLikelihoodCalculatorDRAGEN(ploidy, alleleCount, genotypeTableByPloidy);
    }


    /**
     * Update of shared tables.
     *
     * @param requestedMaximumAllele the new requested maximum allele maximum.
     * @param requestedMaximumPloidy the new requested ploidy maximum.
     */
    private synchronized void ensureCapacity(final int requestedMaximumAllele, final int requestedMaximumPloidy) {
        if (requestedMaximumAllele <= maximumAllele && requestedMaximumPloidy <= maximumPloidy) {
            return;
        }

        maximumPloidy = Math.max(maximumPloidy, requestedMaximumPloidy);
        maximumAllele = Math.max(maximumAllele, requestedMaximumAllele);
        logger.debug("Expanding capacity ploidy:" + maximumPloidy + "->" + maximumPloidy + " allele:" +  maximumAllele +"->" + maximumAllele );
        genotypeTableByPloidy = buildGenotypeAlleleCountsTable(maximumPloidy, maximumAllele);
    }

    /**
     * Perform value checks on maximumPloidy and allele passed to diverse methods in this class.
     * <p>
     *     Throws an exception if there is any issues.
     * </p>
     *
     * @param ploidy the maximum ploidy value.
     * @param maximumAllele the maximum allele value.
     *
     * @throws IllegalArgumentException if either value is negative.
     */
    private static void checkPloidyAndMaximumAllele(final int ploidy, final int maximumAllele) {
        Utils.validateArg(ploidy >= 0, () -> "the ploidy provided cannot be negative: " + ploidy);
        Utils.validateArg(maximumAllele >= 0, () -> "the maximum allele index provided cannot be negative: " + maximumAllele);
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
    public int genotypeCount(final int ploidy, final int alleleCount) {

        final int result = calculateGenotypeCountUsingTables(ploidy, alleleCount);
        if (result == GENOTYPE_COUNT_OVERFLOW) {
            final double largeGenotypeCount = Math.pow(10, MathUtils.log10BinomialCoefficient(ploidy + alleleCount - 1, alleleCount - 1));
            throw new IllegalArgumentException(String.format("the number of genotypes is too large for ploidy %d and allele %d: approx. %.0f", ploidy, alleleCount, largeGenotypeCount));
        }
        return result;
    }

    /**
     * Compute the maximally acceptable allele count (ref allele included) given the maximally acceptable genotype count.
     * @param ploidy            sample ploidy
     * @param maxGenotypeCount  maximum number of genotype count used to calculate upper bound on number of alleles given ploidy
     * @throws IllegalArgumentException if {@code ploidy} or {@code alleleCount} is negative.
     * @return                  the maximally acceptable allele count given ploidy and maximum number of genotypes acceptable
     */
    public static int computeMaxAcceptableAlleleCount(final int ploidy, final int maxGenotypeCount){

        checkPloidyAndMaximumAllele(ploidy, ploidy); // a hack to check ploidy makes sense (could duplicate code but choice must be made)

        if (ploidy == 1) {
            return maxGenotypeCount;
        }
        final double log10MaxGenotypeCount = Math.log10(maxGenotypeCount);

        // Math explanation: genotype count is determined by ${P+A-1 \choose A-1}$, this leads to constraint
        // $\log(\frac{(P+A-1)!}{(A-1)!}) \le \log(P!G)$,
        // where $P$ is ploidy, $A$ is allele count, and $G$ is maxGenotypeCount
        // The upper and lower bounds of the left hand side of the constraint are $P \log(A-1+P)$ and $P \log(A)$
        // which require $A$ to be searched in interval $[10^{\log(P!G)/P} - (P-1), 10^{\log(P!G)/P}]$
        // Denote $[10^{\log(P!G)/P}$ as $x$ in the code.

        final double x = Math.pow(10, (MathUtils.log10Factorial(ploidy) + log10MaxGenotypeCount)/ploidy );
        final int lower = (int)Math.floor(x) - ploidy - 1;
        final int upper = (int)Math.ceil(x);
        for(int a=upper; a>=lower; --a){// check one by one

            final double log10GTCnt = MathUtils.log10BinomialCoefficient(ploidy+a-1, a-1);
            if(log10MaxGenotypeCount >= log10GTCnt) {
                return a;
            }
        }
        throw new GATKException("Code should never reach here.");
    }

    private synchronized int calculateGenotypeCountUsingTables(int ploidy, int alleleCount) {
        checkPloidyAndMaximumAllele(ploidy, alleleCount);
        ensureCapacity(alleleCount, ploidy);
        return (int) numberOfGenotpyes(ploidy, alleleCount);
    }
}