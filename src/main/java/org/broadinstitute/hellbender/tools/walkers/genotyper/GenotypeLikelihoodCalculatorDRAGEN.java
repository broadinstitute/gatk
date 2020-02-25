package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.GenotypeUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Arrays;
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

    // PhredScaled adjustment applied to the BQD score
    static final double PHRED_SCALED_ADJUSTMENT_FOR_BQ_SCORE = 2.5;

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
                                                               final double reverseHomopolymerAdjustment,
                                                               final GenotypeLikelihoodCalculators calculators) {
        //TODO put a very stringent check that we have not invalidated the cache because we will be relying on it to get home in the storm
        // First we invalidate the cache
        Utils.validate(sampleLikelihoods == cachedLikelihoodsObject, "There was a mismatch between the sample stored by the genotyper and the one requesed for BQD, this will result in invalid genotyping");
        final double[] outputArray = new double[genotypeCount];
        Arrays.fill(outputArray, Double.POSITIVE_INFINITY);

        //TODO omptize me
        // TODO is this actually refAllele?
        final Allele refAllele = sampleLikelihoods.getAllele(0);

        //Determine the size of an allele page for the readsLikelihoodsByAlleleFrequency table
        final int readCount = sampleLikelihoods.evidenceCount();
        final int alleleDataSize = readCount * (ploidy + 1);

        for(int gtAlleleIndex = 0; gtAlleleIndex < sampleLikelihoods.numberOfAlleles(); gtAlleleIndex++) {
            //TODO how in the hell do i actaully calculate this offset correctly... blehhhhhhh
            //This is crufty, it just so happens that the index of the homozygous genotype corresponds to the maximum genotype count per field.
            //This should be pulled off as a calculator in some genotyping class.
            final int indexForGT = calculators.genotypeCount(ploidy, gtAlleleIndex) - 1;
            double[] readLikelihoodsForGT = readLikelihoodsByGenotypeIndex[indexForGT];

            for(int errorAlleleIndex = 0; errorAlleleIndex < sampleLikelihoods.numberOfAlleles(); errorAlleleIndex++) {
                // We only want to make calls on SNPs for now
                if (sampleLikelihoods.getAllele(gtAlleleIndex) == sampleLikelihoods.getAllele(errorAlleleIndex) ||
                        sampleLikelihoods.getAllele(gtAlleleIndex).length() != refAllele.length() ||
                        sampleLikelihoods.getAllele(errorAlleleIndex).length() != refAllele.length()) {
                    continue;
                }

                // This selects the index for the reads page in the table corresponding to the likelihoods of the read given allele frequency of 1
                final int offsetForReadLikelihoodGivenAlleleIndex = alleleDataSize * errorAlleleIndex + readCount;

                // BQD scores by strand
                double minScoreFoundForwardsStrand = computeBADModelForStrandData(strandForward, forwardHomopolymerAdjustment, readLikelihoodsForGT, offsetForReadLikelihoodGivenAlleleIndex);
                double minScoreFoundReverseStrand = computeBADModelForStrandData(strandReverse, reverseHomopolymerAdjustment, readLikelihoodsForGT, offsetForReadLikelihoodGivenAlleleIndex);

                double probability = Math.pow(10.0, (-1.0 * (minScoreFoundForwardsStrand + minScoreFoundReverseStrand)/10.0));
                //////
                //TODO this is improtant, this is where the piror gets applied
                //////
                //double prior_adjusted_probability = probability * Math.pow(10.0, (-1.0 * (SNP_HOM_PRIOR)/10.0));
                // TODO figure out how to mark this one down for debuggng in the future
                outputArray[indexForGT] = Math.min(outputArray[indexForGT], probability);
            }
        }
        return outputArray;
    }

    /**
     * Helper function that actually manages the math for BQD;
     *
     * This method works by combining the computed genotype scores for reads with the raw allele likelihoods scores for th
     *
     * @param positionSortedReads  Reads pairs objects (Pair<Pair<read,basequalityforSNP>, sampleReadIndex>) objects sorted in the correct order for partitioning.
     *                             This means that the "error" reads in the partition are sorted by read cycle first in the provided list
     * @param homopolymerAdjustment  Penalty to be applied to reads based on the homopolymer run (this should be precomputed for the ref site in quesiton)
     * @param readLikelihoodsForGT  The array corresponding to the log_10 genotype scores for the genotype in question
     * @param offsetForReadLikelihoodGivenAlleleIndex
     * @return phred scale liklihood for a BQD error mode for reads in the given direction according to the offsets requested
     */
    private double computeBADModelForStrandData(final List<Pair<Pair<GATKRead, Integer>, Integer>> positionSortedReads,
                                        final double homopolymerAdjustment,final  double[] readLikelihoodsForGT,
                                        final int offsetForReadLikelihoodGivenAlleleIndex) {
        if (positionSortedReads.isEmpty()) {
            return 0.0; // TODO check up on this
        }

        // Forwards strand tables
        final double[] cumulative_prob_read_given_error_mode_f = new double[positionSortedReads.size() + 1];
        final double[] cumulative_mean_base_quality_phred_adjusted = new double[positionSortedReads.size() + 1];
        final double[] cumulative_homozygous_genotype_score = new double[positionSortedReads.size() + 1];

        int totalBaseQuality = 0;
        // Iterate over the reads and populate the cumulative arrays
        for (int i = 1; i < cumulative_prob_read_given_error_mode_f.length; i++) {
            int readIndex = positionSortedReads.get(i - 1).getRight();

            // Populate the error probability array
            //?????????????:?:?????> tODO check this one for what the scores look like, we may need to pull out of log space
            cumulative_prob_read_given_error_mode_f[i] = -10.0 * ((cumulative_prob_read_given_error_mode_f[i - 1])/-10.0 +
                    (readIndex == -1 ? 0.0 :
                            (BQD_FIXED_DEFAULT_ALPHA * readAlleleLikelihoodByAlleleCount[offsetForReadLikelihoodGivenAlleleIndex + readIndex])
                            + (1 - BQD_FIXED_DEFAULT_ALPHA) * readLikelihoodsForGT[readIndex]));

            // Populate the mean base quality array
            totalBaseQuality += positionSortedReads.get(i).getLeft().getRight();
            cumulative_mean_base_quality_phred_adjusted[i] = Math.max(0, ((1.0 * totalBaseQuality / i) * PHRED_SCALED_ADJUSTMENT_FOR_BQ_SCORE) - homopolymerAdjustment);

            // Calculate the cumulative genotype score
            //?????????????:?:?????> tODO check this one for what the scores look like, we may need to pull out of log space
            cumulative_homozygous_genotype_score[i] = cumulative_homozygous_genotype_score[i - 1] + -10 * readLikelihoodsForGT[readIndex];
        }

        // Now we find the best partitioning N for the forwards evaluation of the data
        double minScoreFound = Double.POSITIVE_INFINITY;
        for (int n = 0; n < cumulative_mean_base_quality_phred_adjusted.length; n++) {
            minScoreFound = Math.min(minScoreFound,
                    cumulative_mean_base_quality_phred_adjusted[n] + cumulative_prob_read_given_error_mode_f[n] + cumulative_homozygous_genotype_score[n]);
        }

        //TODO this will be where i put my debuglogging (if i had any)
        return minScoreFound;
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
    public <A extends Allele> double[] calculateFRDLikelihoods(LikelihoodMatrix<GATKRead, A> sampleLikelihoods, List<Pair<GATKRead, Integer>> strandForward, List<Pair<GATKRead, Integer>> strandReverse, double forwardHomopolymerAdjustment, double reverseHomopolymerAdjustment) {
        //TODO put a very stringent check that we have not invalidated the cache because we will be relying on it to get home in the storm
        // First we invalidate the cache
        Utils.validate(sampleLikelihoods == cachedLikelihoodsObject, "There was a mismatch between the sample stored by the genotyper and the one requesed for BQD, this will result in invalid genotyping");
        final double[] outputArray = new double[genotypeCount];
        Arrays.fill(outputArray, Double.POSITIVE_INFINITY);

//        for(int gtAlleleIndex = 0; gtAlleleIndex < sampleLikelihoods.numberOfAlleles(); gtAlleleIndex++) {
//            //TODO how in the hell do i actaully calculate this offset correctly... blehhhhhhh
//            //This is crufty, it just so happens that the index of the homozygous genotype corresponds to the maximum genotype count per field.
//            //This should be pulled off as a calculator in some genotyping class.
//            final int indexForGT = calculators.genotypeCount(ploidy, gtAlleleIndex) - 1;
//            double[] readLikelihoodsForGT = readLikelihoodsByGenotypeIndex[indexForGT];
//
//            for (int foreignAlleleIndex = 0; foreignAlleleIndex < sampleLikelihoods.numberOfAlleles(); foreignAlleleIndex++) {
//
//
//            }
//        }

        for (genotypeIndex )

    }


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
