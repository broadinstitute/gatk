package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.Arrays;

/**
 * Genotype likelihood calculator utility.
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
     * Maximum possible number of genotypes that this calculator can handle.
     */
    public static final int MAXIMUM_STRONG_REF_GENOTYPE_PER_PLOIDY = 1000;

    /**
     * Mark to indicate genotype-count overflow due to a large number of allele and ploidy;
     */
    static final int GENOTYPE_COUNT_OVERFLOW = -1;

    /**
     * The current maximum allele index supported by the tables.
     * <p>
     *     Its initial value indicates the initial capacity of the shared {@link #alleleFirstGenotypeOffsetByPloidy} table.
     *     Feel free to change it to anything reasonable that is non-negative.
     * </p>
     */
    private int maximumAllele = 1; // its initial value is the initial capacity of the shared tables.

    /**
     * Shared copy of the offset table as described in {@link #buildGenotypeAlleleCountsTable(int, int, int[][])}.
     *
     * This reference holds the largest requested so far in terms of maximum-allele and maximum-ploidy.
     */
    private int[][] alleleFirstGenotypeOffsetByPloidy =
            buildAlleleFirstGenotypeOffsetTable(maximumPloidy, maximumAllele);


    /**
     * Shared table of genotypes give the ploidy sorted by their index in the likelihood array.
     *
     * <p>
     *  Its format is described in {@link #buildGenotypeAlleleCountsTable(int, int, int[][])}.
     * </p>
     */
    private GenotypeAlleleCounts[][] genotypeTableByPloidy =
            buildGenotypeAlleleCountsTable(maximumPloidy,maximumAllele,alleleFirstGenotypeOffsetByPloidy);

    public GenotypeLikelihoodCalculators(){

    }

    /**
     * Build the table with the genotype offsets based on ploidy and the maximum allele index with representation
     * in the genotype.
     * <p>
     * The result is a matrix containing the offset of the first genotype that contain a particular allele
     * stratified by ploidy.
     * <p>
     *     Row (first dimension) represent the ploidy, whereas
     *     the second dimension represents the allele.
     * </p>
     *
     * <p>
     *     Thus the value a position <i>[p][a]</i> indicates how many genotypes of ploidy <i>p</i> there are before the first
     *     one that contains allele <i>a</i>. <br/>
     *
     *     For example, considering ploidy 3 and alleles A, B, C, D, etc ... (indexed 0, 1, 2, ... respectively):
     *     <br/>
     *     [3][A] == [3][0] == 0 as the first genotype AAA contains A.
     *     <br/>
     *     [3][C] == [3][2] == 4 as the first genotype that contains C, AAC follows: AAA AAB ABB BBB
     *     <br/>
     *     [4][D] == [4][3] == 14  as the first genotype that contains D, AAAD follows: AAAA AAAB AABB ABBB BBBB AAAC
     *     AABC ABBC BBBC AACC ABCC BBCC ACCC BCCC CCCC.
     *
     * </p>
     *
     * <p>
     *     This value are calculated recursively as follows:
     * </p>
     * <pre>
     *
     *     Offset[p][a] := Offset[p-1][a] + Offset[p][a-1] when a > 0, p > 0
     *                     0                               when a == 0
     *                     1                               otherwise
     *
     *
     *         0 1 1  1  1  1   1 ...
     *         0 1 2  3  4  5   6 ...
     *         0 1 3  6 10 15  21 ...
     *         0 1 4 10 20 35  56 ...
     *         0 1 5 15 35 70 126 ...
     *         0 ..................
     * </pre>
     *
     * <p>
     *    Note: if someone can come with a close form computable 0(1) (respect to ploidy and allele count)
     *     please let the author know.
     * </p>
     *
     * <p>
     *     The matrix is guaranteed to have as many rows as indicated by {@code maximumPloidy} + 1; the first
     *     row refers to the special case of ploidy == 0, the second row to ploidy 1 and so forth. Thus the ploidy
     *     matches the index.
     * </p>
     * <p>
     *     The matrix is guaranteed to have as many columns as indicate by {@code maximumAllele} + 1. In this case however
     *     the first allele index 0 is a sense allele (typically the reference allele). The reason to have at least the total
     *     genotype count up to allele count {@link @alleleCapacity} that is equal to the offset of the first genotype
     *     of the following allele; thus we need an extra one.
     * </p>
     *
     * <p>
     *     Although it might seem non-sense to have genotypes of ploidy 0. The values in the first row are used when
     *     filling up values in row 1 and so forth so it is present for programmatic convenience.
     *     Offsets in this row are 0 for the first column and 1 for any others.
     * </p>
     *
     * @param maximumPloidy maximum supported ploidy.
     * @param maximumAllele maximum supported allele index.
     *
     * @throws IllegalArgumentException if {@code maximumPloidy} or {@code maximumAllele} is negative.
     *
     * @return never {@code null}, the matrix described with enough information to address
     *       problems concerning up to the requested maximum allele index and ploidy.
     */
    private static int[][] buildAlleleFirstGenotypeOffsetTable(final int maximumPloidy, final int maximumAllele) {
        checkPloidyAndMaximumAllele(maximumPloidy, maximumAllele);
        final int rowCount = maximumPloidy + 1;
        final int colCount = maximumAllele + 1;
        final int[][] result = new int[rowCount][colCount];

        // Ploidy 0 array must be { 0, 1, 1, ...., 1}
        Arrays.fill(result[0], 1, colCount, 1);
        // Now we take care of the rest of ploidies.
        // We leave the first allele offset to it correct value 0 by starting with allele := 1.
        for (int ploidy = 1; ploidy < rowCount; ploidy++) {
            for (int allele = 1; allele < colCount; allele++) {
                result[ploidy][allele] = result[ploidy][allele - 1] + result[ploidy - 1][allele];
                if (result[ploidy][allele] < result[ploidy][allele - 1]) {
                    result[ploidy][allele] = GENOTYPE_COUNT_OVERFLOW;
                }
            }
        }
        return result;
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
     * @param offsetTable an allele first genotype offset table as constructed using {@link #buildAlleleFirstGenotypeOffsetTable(int, int)}
     *                    that supports at least up to {@code maximumAllele} and {@code maximumPloidy}.
     *
     * @throws IllegalArgumentException if {@code maximumPloidy} or {@code maximumAllele} is negative, or {@code offsetTable} is {@code null},
     *   or it does not have the capacity to handle the requested maximum ploidy or allele index.
     *
     * @return never {@code null}.
     */
    private static GenotypeAlleleCounts[][] buildGenotypeAlleleCountsTable(final int maximumPloidy, final int maximumAllele, final int[][] offsetTable) {
        checkPloidyAndMaximumAllele(maximumPloidy, maximumAllele);
        checkOffsetTableCapacity(offsetTable,maximumPloidy,maximumAllele);
        final int rowCount = maximumPloidy + 1;
        final GenotypeAlleleCounts[][] result = new GenotypeAlleleCounts[rowCount][]; // each row has a different number of columns.

        for (int ploidy = 0; ploidy <= maximumPloidy; ploidy++) {
            result[ploidy] = buildGenotypeAlleleCountsArray(ploidy, maximumAllele, offsetTable);
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
     * @param genotypeOffsetTable table with the offset of the first genotype that contain an allele given
     *                            the ploidy and its index.
     *
     * @throws IllegalArgumentException if {@code ploidy} or {@code length} is negative.
     *
     * @return never {@code null}, follows the specification above.
     */
    private static GenotypeAlleleCounts[] buildGenotypeAlleleCountsArray(final int ploidy, final int alleleCount, final int[][] genotypeOffsetTable) {
        Utils.validateArg(ploidy >= 0, () -> "the requested ploidy cannot be negative: " + ploidy);
        Utils.validateArg(alleleCount >= 0, () -> "the requested maximum allele cannot be negative: " + alleleCount);
        final int length = genotypeOffsetTable[ploidy][alleleCount];
        final int strongRefLength = length == GENOTYPE_COUNT_OVERFLOW ? MAXIMUM_STRONG_REF_GENOTYPE_PER_PLOIDY : Math.min(length, MAXIMUM_STRONG_REF_GENOTYPE_PER_PLOIDY);
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
    public GenotypeLikelihoodCalculator getInstance(final int ploidy, final int alleleCount) {
        checkPloidyAndMaximumAllele(ploidy, alleleCount);

        if (calculateGenotypeCountUsingTables(ploidy, alleleCount) == GENOTYPE_COUNT_OVERFLOW) {
            final double largeGenotypeCount = Math.pow(10, MathUtils.log10BinomialCoefficient(ploidy + alleleCount - 1, alleleCount - 1));
            throw new IllegalArgumentException(String.format("the number of genotypes is too large for ploidy %d and allele %d: approx. %.0f", ploidy, alleleCount, largeGenotypeCount));
        }

        // At this point the tables must have at least the requested capacity, likely to be much more.
        return new GenotypeLikelihoodCalculator(ploidy, alleleCount, alleleFirstGenotypeOffsetByPloidy, genotypeTableByPloidy);
    }

    /**
     * Update of shared tables
     *
     * @param requestedMaximumAllele the new requested maximum allele maximum.
     * @param requestedMaximumPloidy the new requested ploidy maximum.
     */
    private void ensureCapacity(final int requestedMaximumAllele, final int requestedMaximumPloidy) {

        final boolean needsToExpandAlleleCapacity = requestedMaximumAllele > maximumAllele;
        final boolean needsToExpandPloidyCapacity = requestedMaximumPloidy > maximumPloidy;

        // Double check with the lock on to avoid double work.
        if (!needsToExpandAlleleCapacity && !needsToExpandPloidyCapacity) {
            return;
        }

        final int newMaximumPloidy = Math.max(maximumPloidy, requestedMaximumPloidy);
        final int newMaximumAllele = Math.max(maximumAllele, requestedMaximumAllele);

        logger.debug("Expanding capacity ploidy:" + maximumPloidy + "->" + newMaximumPloidy + " allele:" +  maximumAllele +"->" + newMaximumAllele );

        // Update tables first.
        alleleFirstGenotypeOffsetByPloidy = buildAlleleFirstGenotypeOffsetTable(newMaximumPloidy,newMaximumAllele);
        genotypeTableByPloidy = buildGenotypeAlleleCountsTable(newMaximumPloidy,newMaximumAllele,alleleFirstGenotypeOffsetByPloidy);

        if (needsToExpandAlleleCapacity) {
            maximumAllele = requestedMaximumAllele;
        }
        if (needsToExpandPloidyCapacity) {
            maximumPloidy = requestedMaximumPloidy;
        }
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

    private static void checkOffsetTableCapacity(final int[][] offsetTable, final int maximumPloidy, final int maximumAllele) {
        Utils.nonNull(offsetTable, "the allele first genotype offset table provided cannot be null");
        Utils.validateArg(offsetTable.length > maximumPloidy, () -> "the allele first genotype offset table provided does not have enough " +
                    "capacity for requested maximum ploidy: " + maximumPloidy);
        Utils.validateArg(offsetTable[0].length >= maximumAllele, () -> "the allele first genotype offset table provided does not have enough " +
                    "capacity for requested maximum allele index: " + maximumAllele);
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

    private int calculateGenotypeCountUsingTables(int ploidy, int alleleCount) {
        checkPloidyAndMaximumAllele(ploidy, alleleCount);
        if (ploidy > maximumPloidy || alleleCount > maximumAllele) {
            ensureCapacity(alleleCount, ploidy);
        }
        return alleleFirstGenotypeOffsetByPloidy[ploidy][alleleCount];
    }
}