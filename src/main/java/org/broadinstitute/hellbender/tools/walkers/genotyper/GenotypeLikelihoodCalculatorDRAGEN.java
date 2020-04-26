package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.apache.commons.collections.ListUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

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

    private final double cachedLog10Alpha;
    private final double cachedLog10AlphaInverse;

    /**
     * Creates a new calculator providing its ploidy and number of genotyping alleles.
     */
    protected GenotypeLikelihoodCalculatorDRAGEN(final int ploidy, final int alleleCount,
                                                 final int[][] alleleFirstGenotypeOffsetByPloidy,
                                                 final GenotypeAlleleCounts[][] genotypeTableByPloidy) {
        super(ploidy, alleleCount, alleleFirstGenotypeOffsetByPloidy, genotypeTableByPloidy);
        Utils.validateArg(ploidy > 0, () -> "ploidy must be at least 1 but was " + ploidy);
        // The number of possible components is limited by distinct allele count and ploidy.
        cachedLog10Alpha = Math.log10(BQD_FIXED_DEFAULT_ALPHA);
        cachedLog10AlphaInverse = Math.log10(1 - BQD_FIXED_DEFAULT_ALPHA);
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
                                                               final List<DRAGENBQDGenotypesModel.DragenReadContainer> strandForward,
                                                               final List<DRAGENBQDGenotypesModel.DragenReadContainer> strandReverse,
                                                               final double forwardHomopolymerAdjustment,
                                                               final double reverseHomopolymerAdjustment,
                                                               final GenotypeLikelihoodCalculators calculators) {
        //TODO put a very stringent check that we have not invalidated the cache because we will be relying on it to get home in the storm
        // First we invalidate the cache
        Utils.validate(sampleLikelihoods == cachedLikelihoodsObject, "There was a mismatch between the sample stored by the genotyper and the one requesed for BQD, this will result in invalid genotyping");
        final double[] outputArray = new double[genotypeCount];
        Arrays.fill(outputArray, Double.NEGATIVE_INFINITY);


        final Allele refAllele = sampleLikelihoods.getAllele(0);

        //Determine the size of an allele page for the readsLikelihoodsByAlleleFrequency table
        final int readCount = sampleLikelihoods.evidenceCount();
        final int alleleDataSize = readCount * (ploidy + 1);

        for(int gtAlleleIndex = 0; gtAlleleIndex < sampleLikelihoods.numberOfAlleles(); gtAlleleIndex++) {
            //TODO how in the hell do i actaully calculate this offset correctly... blehhhhhhh
            //This is crufty, it just so happens that the index of the homozygous genotype corresponds to the maximum genotype count per field.
            //This should be pulled off as a calculator in some genotyping class.
            final int indexForGT = calculators.genotypeCount(ploidy, gtAlleleIndex + 1) - 1;
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
                double minScoreFoundForwardsStrand = computeBQDModelForStrandData(strandForward, forwardHomopolymerAdjustment, readLikelihoodsForGT, offsetForReadLikelihoodGivenAlleleIndex, true);
                double minScoreFoundReverseStrand = computeBQDModelForStrandData(strandReverse, reverseHomopolymerAdjustment, readLikelihoodsForGT, offsetForReadLikelihoodGivenAlleleIndex, false);

                double modelScoreInLog10 = (minScoreFoundForwardsStrand + minScoreFoundReverseStrand)/-10.0;
                //////
                // NOTE we have not applied the prior here, this is because that gets applied downstream to each genotype in the array.
                // since the prior is applied evenly to the defualt and error models this should not change this selection here.
                outputArray[indexForGT] = Math.max(outputArray[indexForGT], modelScoreInLog10);
            }
        }
        return outputArray;
    }

    /**
     * Helper function that actually manages the math for BQD;
     *
     * This method works by combining the computed genotype scores for reads with the raw allele likelihoods scores for th
     *
     * @param positionSortedReads  Reads pairs objects (Pair<Pair<read,readBaseOffset>, sampleReadIndex>) objects sorted in the correct order for partitioning.
     *                             This means that the "error" reads in the partition are sorted by read cycle first in the provided list
     * @param homopolymerAdjustment  Penalty to be applied to reads based on the homopolymer run (this should be precomputed for the ref site in quesiton)
     * @param readLikelihoodsForGT  The array corresponding to the log_10 genotype scores for the genotype in question
     * @param offsetForReadLikelihoodGivenAlleleIndex
     * @return phred scale liklihood for a BQD error mode for reads in the given direction according to the offsets requested
     */
    private double computeBQDModelForStrandData(final List<DRAGENBQDGenotypesModel.DragenReadContainer> positionSortedReads,
                                                final double homopolymerAdjustment, final double[] readLikelihoodsForGT,
                                                final int offsetForReadLikelihoodGivenAlleleIndex, final boolean forwards) {
        if (positionSortedReads.isEmpty()) {
            return 0.0; // TODO check up on this
        }

        // Forwards strand tables (all in phred space for the sake of conveient debugging with provided scripts)
        final double[] cumulative_p_R_for_E = new double[positionSortedReads.size() + 1];
        final double[] cumulative_mean_base_quality_phred_adjusted = new double[positionSortedReads.size() + 1];
        final double[] cumulative_P_GT = new double[positionSortedReads.size() + 1];

        double totalBaseQuality = 0;
        int baseQualityDenominator = 0; // We track this seperately because not every read overlaps the SNP in quesiton due to padding.
        // Iterate over the reads and populate the cumulative arrays
        for (int i = 1; i < cumulative_p_R_for_E.length; i++) {
            final DRAGENBQDGenotypesModel.DragenReadContainer container = positionSortedReads.get(i - 1);
            final int readIndex = container.getIndexInLikelihoodsObject();

            // Retrieve the homozygous genotype score and the error allele scores for the read in question
            final double homozygousGenotypeContribution;
            final double errorAlleleContribution;
            if (readIndex != -1) {
                homozygousGenotypeContribution = readLikelihoodsForGT[readIndex] - MathUtils.log10(2);
                errorAlleleContribution = readAlleleLikelihoodByAlleleCount[offsetForReadLikelihoodGivenAlleleIndex + readIndex];
            } else {
                // If read index == -1 then we are evaluating a read that was rejected by the HMM and therefore doesn't have genotype scores
                homozygousGenotypeContribution = 0;
                errorAlleleContribution = 0;
            }

//            System.out.println("read index:"+readIndex+" ihomGT:"+-10*homozygousGenotypeContribution+"  errorAlleleContribution:"+-10*errorAlleleContribution + " difference in phred: "+ -10*(errorAlleleContribution - homozygousGenotypeContribution));

            // Populate the error probability array in phred space
            // Calculation: Alpha * P(r|E_allele) + (1 - Alpha) * P(r | G_homozygousGT))
            double phredContributionForRead = (homozygousGenotypeContribution==0 && errorAlleleContribution==0) ? 0 : -10 *
                    MathUtils.approximateLog10SumLog10(errorAlleleContribution + cachedLog10Alpha,
                                                       homozygousGenotypeContribution + cachedLog10AlphaInverse);
            cumulative_p_R_for_E[i] = cumulative_p_R_for_E[i-1] + phredContributionForRead;

            // Populate the cumulative genotype contribution array with the score for this read
            // Calculation: (P(r | G_A1) + P(r | G_A2)) / 2
            cumulative_P_GT[i] = cumulative_P_GT[i - 1] + -10 * homozygousGenotypeContribution;

            // Populate the mean base quality array
            if (container.hasValidBaseQuality()) {
                totalBaseQuality += container.getBaseQuality();
                baseQualityDenominator++;
            }
            cumulative_mean_base_quality_phred_adjusted[i] = Math.max(0,
                    ((totalBaseQuality / (baseQualityDenominator == 0 ? 1 : baseQualityDenominator)) * PHRED_SCALED_ADJUSTMENT_FOR_BQ_SCORE) - homopolymerAdjustment);
        }

        // Now we find the best partitioning N for the forwards evaluation of the data
        double minScoreFound = Double.POSITIVE_INFINITY;
        int nIndexUsed = 0;
        double lastProbE=0;
        double lastGQQual=0;
        for (int n = 0; n < cumulative_mean_base_quality_phred_adjusted.length; n++) {
            final double bqdScore = cumulative_mean_base_quality_phred_adjusted[n] + cumulative_p_R_for_E[n] + (cumulative_P_GT[cumulative_P_GT.length-1] - cumulative_P_GT[n]);
//            System.out.println(String.format("n=%d: %.2f, cum_phred_bq=%.2f, cum_prob_r_Error=%.2f, prob_G_remaining=%.2f",
//                    n, bqdScore, cumulative_mean_base_quality_phred_adjusted[n], cumulative_p_R_for_E[n],
//                    (cumulative_P_GT[cumulative_P_GT.length-1] - cumulative_P_GT[n])));
            if (minScoreFound > bqdScore) {
                minScoreFound = bqdScore;
                nIndexUsed = n;
            }
        }

        // Debug output for the genotyper to see into the calculation itself
//        System.out.println(String.format("n%d=%2d, best_phred_score =%5.2f q_mean=%5.2f, alpha=%4.2f, Ph(E)=%4.2f;  ", forwards?1:0,
//                nIndexUsed, minScoreFound, cumulative_mean_base_quality_phred_adjusted[nIndexUsed],
//                BQD_FIXED_DEFAULT_ALPHA, cumulative_p_R_for_E[nIndexUsed]));

        //TODO this will be where i put my debuglogging (if i had any)
        return minScoreFound;
    }

    /**
     *
     * @param sampleLikelihoods
     * @param strandForward
     * @param strandReverse
     * @param forwardHomopolymerAdjustment
     * @param snipAprioriHet
     * @return
     */
    //TODO this WILLLL need to handle indel priors
    public <A extends Allele> double[] calculateFRDLikelihoods(LikelihoodMatrix<GATKRead, A> sampleLikelihoods,
                                                               final List<DRAGENBQDGenotypesModel.DragenReadContainer> readContainers,
                                                               final double snipAprioriHet, final double indelAprioriHet, final GenotypeLikelihoodCalculators calculators) {
        Utils.validate(sampleLikelihoods == cachedLikelihoodsObject, "There was a mismatch between the sample stored by the genotyper and the one requesed for BQD, this will result in invalid genotyping");
        final double[] outputArray = new double[genotypeCount];
        Arrays.fill(outputArray, Double.NEGATIVE_INFINITY);

        final Allele refAllele = sampleLikelihoods.getAllele(0);

        //Determine the size of an allele page for the readsLikelihoodsByAlleleFrequency table
        final int readCount = sampleLikelihoods.evidenceCount();
        final int alleleDataSize = readCount * (ploidy + 1);

        for(int fAlleleIndex = 0; fAlleleIndex < sampleLikelihoods.numberOfAlleles(); fAlleleIndex++) {
            // ignore symbolic alleles
            if (sampleLikelihoods.getAllele(fAlleleIndex).isSymbolic() ) {
                continue;
            }
            final boolean isIndel = sampleLikelihoods.getAllele(fAlleleIndex).length() != refAllele.length();
            final int offsetForReadLikelihoodGivenAlleleIndex = alleleDataSize * fAlleleIndex + readCount;

            // Here we generate a set of the critical log10(P(F)) values that we will iterate over
            final Set<Double> criticalThresholdsForward = new HashSet<>();
            final Set<Double> criticalThresholdsReverse = new HashSet<>();
            final Set<Double> criticalThresholds = new HashSet<>();
            computeCriticalValues(criticalThresholdsForward, criticalThresholdsReverse, criticalThresholds, readContainers, fAlleleIndex == 0 ? 0 : (isIndel? indelAprioriHet : snipAprioriHet)/-10); // simplified in line with DRAGEN, uses 1 alleledist for both snp and indels
            final List<Double> criticalThresholdsSortedForward = criticalThresholdsForward.stream().sorted(Double::compareTo).collect(Collectors.toList());
            final List<Double> criticalThresholdsSortedReverse = criticalThresholdsReverse.stream().sorted(Double::compareTo).collect(Collectors.toList());
            final List<Double> criticalThresholdsSorted = criticalThresholds.stream().sorted(Double::compareTo).collect(Collectors.toList());
            System.out.println("fIndex: "+fAlleleIndex+" criticalValues: \n"+criticalThresholds.stream().map(d -> Double.toString(d)).collect(Collectors.joining("\n")));

            // iterate over all of the homozygous genotypes for the given allele
            for(int gtAlleleIndex = 0; gtAlleleIndex < sampleLikelihoods.numberOfAlleles(); gtAlleleIndex++) {
                // Skip over the allele corresponding to the "foreign" allele
                if (gtAlleleIndex == fAlleleIndex) {
                    continue;
                }

                //This is crufty, it just so happens that the index of the homozygous genotype corresponds to the maximum genotype count per field.
                //This should be pulled off as a calculator in some genotyping class.
                final int indexForGT = calculators.genotypeCount(ploidy, gtAlleleIndex + 1) - 1;
                double[] readLikelihoodsForGT = readLikelihoodsByGenotypeIndex[indexForGT];

                // TODO this will have to account for the indel priors later in its existance
                final double gtAllelePrior = 0.0;

                System.out.println("\nForwards Strands: ");
                final double maxLog10FForwardsStrand = computeFRDModelForStrandData(readContainers, c -> !c.isReverseStrand() , offsetForReadLikelihoodGivenAlleleIndex, readLikelihoodsForGT, criticalThresholdsSortedForward, gtAllelePrior);
                System.out.println("\nReverse Strands: ");
                final double maxLog10FReverseStrand = computeFRDModelForStrandData(readContainers, c -> c.isReverseStrand(), offsetForReadLikelihoodGivenAlleleIndex, readLikelihoodsForGT, criticalThresholdsSortedReverse, gtAllelePrior);
                System.out.println("\nBoth Strands: ");
                final double maxLog10FBothStrands = computeFRDModelForStrandData(readContainers, c -> true, offsetForReadLikelihoodGivenAlleleIndex, readLikelihoodsForGT, criticalThresholdsSorted, gtAllelePrior);

//                System.out.println("gtAlleleIndex : "+gtAlleleIndex+ " fAlleleIndex: "+fAlleleIndex +" forwards: "+maxLog10FForwardsStrand+" reverse: "+maxLog10FReverseStrand+" both: "+maxLog10FBothStrands);
                outputArray[indexForGT] = Math.max(outputArray[indexForGT], maxLog10FForwardsStrand);
                outputArray[indexForGT] = Math.max(outputArray[indexForGT], maxLog10FReverseStrand);
                outputArray[indexForGT] = Math.max(outputArray[indexForGT], maxLog10FBothStrands);
            }

        }

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

//        for (genotypeIndex )
        return outputArray;
    }


    /**
     *
      * @param positionSortedReads
     * @param predicate
     * @param offsetForReadLikelihoodGivenAlleleIndex
     * @param readLikelihoodsForGT
     * @param criticalThresholdsSorted
     * @param gtAllelePrior
     * @return log10 score for frd model given the data and the predicate
     */
    private double computeFRDModelForStrandData(final List<DRAGENBQDGenotypesModel.DragenReadContainer> positionSortedReads, final Predicate<DRAGENBQDGenotypesModel.DragenReadContainer> predicate,
                                                final int offsetForReadLikelihoodGivenAlleleIndex, final double[] readLikelihoodsForGT, final List<Double> criticalThresholdsSorted, final double gtAllelePrior) {
        if (positionSortedReads.isEmpty()) {
            return Double.NEGATIVE_INFINITY; // TODO check up on this
        }
        System.out.println("GT Allele Prior: "+gtAllelePrior);

        int counter = 0;
        double maxLpspi = Double.NEGATIVE_INFINITY;

        for (final Double lpf : criticalThresholdsSorted) {
            double f_ratio = 0.0;
            double f_denom = 0.0;
            double localMaxLpspi = Double.NEGATIVE_INFINITY;

            // iterate over the
            for (final DRAGENBQDGenotypesModel.DragenReadContainer container : positionSortedReads) {
                final int readIndex = container.getIndexInLikelihoodsObject();

                // Keep track of the aggregate support for the foreign allele
                if (predicate.test(container)) {
                    double LPd_r_F = container.getPhredPFValue() - 0.0000001 <= lpf ?
                            readAlleleLikelihoodByAlleleCount[offsetForReadLikelihoodGivenAlleleIndex + readIndex] :
                            Double.NEGATIVE_INFINITY;
                    double lp_r_GT = readLikelihoodsForGT[readIndex] - MathUtils.log10(2);

                    f_ratio += Math.pow(10, LPd_r_F - MathUtils.approximateLog10SumLog10(LPd_r_F, lp_r_GT));
                    f_denom++;
                }
            }

            // Don't learn the beta but approximate it based on the read support for the alt
            double beta = Math.min(f_ratio / f_denom, 0.5);
            double log10beta = Math.log10(beta);
            double log10betaDelta = Math.log10(1.0 - beta);

            double LP_R_GF = 0.0;

            // iterate over the containers again using the approximated beta constraint
            for (final DRAGENBQDGenotypesModel.DragenReadContainer container : positionSortedReads) {
                final int readIndex = container.getIndexInLikelihoodsObject();

                double lp_r_GT = readLikelihoodsForGT[readIndex] - MathUtils.log10(2);

                // COMPUTE THE MODEL FOR THE STRAND IN QUESTION
                if (predicate.test(container)) {
                    double LPd_r_F = container.getPhredPFValue() - 0.0000001 <= lpf ?
                            readAlleleLikelihoodByAlleleCount[offsetForReadLikelihoodGivenAlleleIndex + readIndex] :
                            Double.NEGATIVE_INFINITY;

                    LP_R_GF += MathUtils.approximateLog10SumLog10(log10beta + LPd_r_F, log10betaDelta + lp_r_GT);
                } else {
                    LP_R_GF += lp_r_GT;
                }
            }
            // Allele prior for error allele, plus posterior for foreign event, plus model posterior
            double LPsi = gtAllelePrior + lpf + LP_R_GF;
            localMaxLpspi = Math.max(localMaxLpspi, LPsi);

            System.out.println("beta: "+beta+" localMaxLpspi: " + localMaxLpspi + " for lpf: "+lpf+" with LP_R_GF: "+LP_R_GF+" index: "+counter++);
            maxLpspi = Math.max(maxLpspi, localMaxLpspi);
        }

        return maxLpspi;
    }


    // helper method to populate the reads containers properly with their critical values and store them in the provided set
    // NOTE: this has the side effect of setting the DragenReadContainer setPhredPFValue() values for the reads for the given set of alleles
    private void computeCriticalValues(final Set<Double> criticalThresholdsForwards, final Set<Double> criticalThresholdsReverse, final Set<Double> criticalThresholdsTotal, final List<DRAGENBQDGenotypesModel.DragenReadContainer> container, final double log10MapqPriorAdjustment) {
        for (int i = 0; i < container.size(); i++) {
            final DRAGENBQDGenotypesModel.DragenReadContainer readContainer = container.get(i);
            final double log10CriticalValue = readContainer.getMappingQuality() / -10.0 + log10MapqPriorAdjustment;
            readContainer.setPhredPFValue(log10CriticalValue);
            // Split the critical thresholds up by their applicable strands in order to avoid repeated work
            if (readContainer.isReverseStrand()) {
                criticalThresholdsReverse.add(log10CriticalValue);
            } else {
                criticalThresholdsForwards.add(log10CriticalValue);
            }
            criticalThresholdsTotal.add(log10CriticalValue);
        }
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
