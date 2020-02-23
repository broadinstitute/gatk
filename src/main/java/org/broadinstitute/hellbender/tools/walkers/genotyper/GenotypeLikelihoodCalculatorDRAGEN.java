package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;

/**
 * Helper to calculate genotype likelihoods for DRAGEN advanced genotyping models (BQD - Base Quality Dropout, and FRD - Foreign Reads Detection).
 *
 * This object is simply a thin wrapper on top of a regular GenotypeLikelihoods object with some extra logic for handling new impouts to the genotyper:
 *  - both BQD and FRD rely on per-read per-genotype scores as would be computed for the standard genotyper, rather than pay the cost of recomputing these
 *    for each of the 3 independent models this GenotypeLikelihoodCalculator simply makes the computation once and relies on the fact that the underlying
 *    readLikelihoodsByGenotypeIndex is still populated from the previous call. To this end strict object equality tests have been implemented to ensure
 *    that the cache is populated with the correct likelihoods before running either of the advanced models.
 */
public final class GenotypeLikelihoodCalculatorDRAGEN extends GenotypeLikelihoodCalculator {
    // measure of the likelihood of an error base occurring if we have triggered a base quality dropout
    static final double BQD_FIXED_DEFAULT_ALPHA = 0.5;

    // Cache for enforcing the strictness of using the filled array with the correct likelihoods object
    private LikelihoodMatrix<?, ?> cachedLikelihoodsObject = null;

    /**
     * Creates a new calculator providing its ploidy and number of genotyping alleles.
     */
    protected GenotypeLikelihoodCalculatorDRAGEN(final int ploidy, final int alleleCount,
                                                 final int[][] alleleFirstGenotypeOffsetByPloidy,
                                                 final GenotypeAlleleCounts[][] genotypeTableByPloidy) {
        super(ploidy, alleleCount, alleleFirstGenotypeOffsetByPloidy, genotypeTableByPloidy);
        Utils.validateArg(ploidy > 0, () -> "ploidy must be at least 1 but was " + ploidy);
        // The number of possible components is limited by distinct allele count and ploidy.
    }

    /**
     *
     * @param sampleLikelihoods
     * @param strandForward
     * @param strandReverse
     * @param forwardHomopolymerAdjustment
     * @param reverseHomopolymerAdjustment
     * @return
     */
    public <A extends Allele> double[] calculateBQDLikelihoods(final LikelihoodMatrix<GATKRead, A> sampleLikelihoods,
                                                               final List<Pair<Pair<GATKRead,Integer>, Integer>> strandForward,
                                                               final List<Pair<Pair<GATKRead,Integer>, Integer>> strandReverse,
                                                               final double forwardHomopolymerAdjustment,
                                                               final double reverseHomopolymerAdjustment) {
        //TODO put a very stringent check that we have not invalidated the cache because we will be relying on it to get home in the storm
        // First we invalidate the cache
        Utils.validate(sampleLikelihoods == cachedLikelihoodsObject, "There was a mismatch between the sample stored by the genotyper and the one requesed for BQD, this will result in invalid genotyping");

        //TODO omptize me
        // TODO is this actually refAllele?
        final Allele refAllele = sampleLikelihoods.getAllele(0);

        //Determine the size of an allele page for the readsLikelihoodsByAlleleFrequency table
        final int readCount = sampleLikelihoods.evidenceCount();
        final int alleleDataSize = readCount * (ploidy + 1);

        for(int gtAlleleIndex = 0; gtAlleleIndex < sampleLikelihoods.numberOfAlleles(); gtAlleleIndex++) {
            //TODO how in the hell do i actaully calculate this offset correctly... blehhhhhhh
            double[] readLikelihoodsForGT = readLikelihoodsByGenotypeIndex[getGTMatrixOffsetForAlleleIndex(gtAlleleIndex)];

            for(int errorAlleleIndex = 0; errorAlleleIndex < sampleLikelihoods.numberOfAlleles(); errorAlleleIndex++) {
                // We only want to make calls on SNPs for now
                if (sampleLikelihoods.getAllele(gtAlleleIndex) == sampleLikelihoods.getAllele(errorAlleleIndex) ||
                        sampleLikelihoods.getAllele(gtAlleleIndex).length() != refAllele.length() ||
                        sampleLikelihoods.getAllele(errorAlleleIndex).length() != refAllele.length()) {
                    continue;
                }

                // This selects the index for the reads page in the table corresponding to the likelihoods of the read given allele frequency of 1
                int offsetForReadLikelihoodGivenAlleleIndex = alleleDataSize * errorAlleleIndex + readCount;

                // Forwards strand computations
                double[] prob_read_given_error_mode_f = new double[strandForward.size()];
                for (int i = 0; i < prob_read_given_error_mode_f.length; i++) {
                    int readIndex = strandForward.get(i).getRight();
                    if (readIndex == -1) {
                        prob_read_given_error_mode_f[i] = 0.0;
                    } else {
                        // P(R_i | E) = alpha * P(R_i | ErrorAllele) + (1-alpha) * P(R | Gt)
                        //TODO the first one is in logspace the second is probably also but must be clear
                        prob_read_given_error_mode_f[i] =
                                (BQD_FIXED_DEFAULT_ALPHA * readAlleleLikelihoodByAlleleCount[offsetForReadLikelihoodGivenAlleleIndex + readIndex])
                                + (1 - BQD_FIXED_DEFAULT_ALPHA) * readLikelihoodsForGT[readIndex];
                    }
                }
            }
        }
    }

    /**
     *
     * @param sampleLikelihoods
     * @param strandForward
     * @param strandReverse
     * @param forwardHomopolymerAdjustment
     * @param reverseHomopolymerAdjustment
     * @return
//     */
//    public <A extends Allele> double[] calculateFRDLikelihoods(LikelihoodMatrix<GATKRead, A> sampleLikelihoods, List<Pair<GATKRead, Integer>> strandForward, List<Pair<GATKRead, Integer>> strandReverse, double forwardHomopolymerAdjustment, double reverseHomopolymerAdjustment) {
//        //TODO put a very stringent check that we have not invalidated the cache because we will be relying on it to get home in the storm
//        // First we invalidate the cache
//        Utils.validate(sampleLikelihoods == cachedLikelihoodsObject, "There was a mismatch between the sample stored by the genotyper and the one requesed for BQD, this will result in invalid genotyping");
//
//    }


    /**
     * See {@link GenotypeLikelihoodCalculator#genotypeLikelihoods}. This wrapper just enforces that the likelihoods object is recorded in the cache.
     *
     * @return never {@code null}.
     */
    public <EVIDENCE, A extends Allele> GenotypeLikelihoods genotypeLikelihoods(final LikelihoodMatrix<EVIDENCE, A> likelihoods) {
        cachedLikelihoodsObject = null;
        GenotypeLikelihoods output = super.genotypeLikelihoods(likelihoods);
        cachedLikelihoodsObject = likelihoods;
        return output;
    }

    /**
     * See {@link GenotypeLikelihoodCalculator#getReadRawReadLikelihoodsByGenotypeIndex}. This wrapper just enforces that the likelihoods object is recorded in the cache.
     *
     * @return never {@code null}.
     */
    public <EVIDENCE, A extends Allele> double[] rawGenotypeLikelihoods(final LikelihoodMatrix<EVIDENCE, A> likelihoods) {
        cachedLikelihoodsObject = null;
        double[] output = super.getReadRawReadLikelihoodsByGenotypeIndex(likelihoods);
        cachedLikelihoodsObject = likelihoods;
        return output;
    }
}
