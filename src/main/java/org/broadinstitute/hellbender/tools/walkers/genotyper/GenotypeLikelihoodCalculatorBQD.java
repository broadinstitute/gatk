package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;

/**
 * Helper to calculate genotype likelihoods given a ploidy and an allele count (number of possible distinct alleles).
 */
public final class GenotypeLikelihoodCalculatorBQD extends GenotypeLikelihoodCalculator {

    /**
     * Creates a new calculator providing its ploidy and number of genotyping alleles.
     */
    protected GenotypeLikelihoodCalculatorBQD(final int ploidy, final int alleleCount,
                                              final int[][] alleleFirstGenotypeOffsetByPloidy,
                                              final GenotypeAlleleCounts[][] genotypeTableByPloidy) {
        super(ploidy, alleleCount, alleleFirstGenotypeOffsetByPloidy, genotypeTableByPloidy);
        Utils.validateArg(ploidy > 0, () -> "ploidy must be at least 1 but was " + ploidy);
        // The number of possible components is limited by distinct allele count and ploidy.
    }


    /**
     * Calculates the likelihood component of each read on each genotype.
     *
     * @param readLikelihoodComponentsByAlleleCount [a][f][r] likelihood stratified by allele <i>a</i>, frequency in genotype <i>f</i> and
     *                                              read <i>r</i>.
     * @param readCount                             number of reads in {@code readLikelihoodComponentsByAlleleCount}.
     * @return never {@code null}.
     */
    protected double[][] genotypeLikelihoodByRead(final double[] readLikelihoodComponentsByAlleleCount, final int readCount) {

        // Here we don't use the convenience of {@link #genotypeAlleleCountsAt(int)} within the loop to spare instantiations of
        // GenotypeAlleleCounts class when we are dealing with many genotypes.
        GenotypeAlleleCounts alleleCounts = genotypeAlleleCounts[0];

//        for (int genotypeIndex = 0; genotypeIndex < genotypeCount; genotypeIndex++) {
//            final double[] readLikelihoods = this.readLikelihoodsByGenotypeIndex[genotypeIndex];
//            final int componentCount = alleleCounts.distinctAlleleCount();
//            switch (componentCount) {
//                case 1: //
//                    singleComponentGenotypeLikelihoodByRead(alleleCounts, readLikelihoods, readLikelihoodComponentsByAlleleCount, readCount);
//                    break;
//                case 2:
//                    twoComponentGenotypeLikelihoodByRead(alleleCounts,readLikelihoods,readLikelihoodComponentsByAlleleCount, readCount);
//                    break;
//                default:
//                   throw new GATKException.ShouldNeverReachHereException("BQD Genotyping model is only applied for Diploid Samples ");
//            }
//            if (genotypeIndex < genotypeCount - 1) {
//                alleleCounts = nextGenotypeAlleleCounts(alleleCounts);
//            }
//        }
        return readLikelihoodsByGenotypeIndex;
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
//    public <EVIDENCE, A extends Allele> GenotypeLikelihoods genotypeLikelihoods(final LikelihoodMatrix<EVIDENCE, A> likelihoods, final List< StrandForwardSorted,) {
//        Utils.nonNull(likelihoods);
//        Utils.validateArg(likelihoods.numberOfAlleles() == alleleCount, "mismatch between allele list and alleleCount");
//        final int readCount = likelihoods.evidenceCount();
//        ensureReadCapacity(readCount);
//
//        /// [x][y][z] = z * LnLk(Read_x | Allele_y)
//        final double[] readLikelihoodComponentsByAlleleCount
//                = readLikelihoodComponentsByAlleleCount(likelihoods);
//        final double[][] genotypeLikelihoodByRead = genotypeLikelihoodByRead(readLikelihoodComponentsByAlleleCount,readCount);
//        final double[] readLikelihoodsByGenotypeIndex = genotypeLikelihoods(genotypeLikelihoodByRead, readCount);
//        return GenotypeLikelihoods.fromLog10Likelihoods(readLikelihoodsByGenotypeIndex);
//    }
}
