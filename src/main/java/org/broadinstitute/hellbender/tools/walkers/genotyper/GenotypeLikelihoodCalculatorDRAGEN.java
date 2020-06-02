package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;

/**
 * Helper to calculate genotype likelihoods for DRAGEN advanced genotyping models (BQD - Base Quality Dropout, and FRD - Foreign Reads Detection).
 *
 * This object is simply a thin wrapper on top of a regular GenotypeLikelihoods object with some extra logic for handling new inputs to the genotyper:
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

    // Debug output stream to be managed by the HaplotypeCallerEngine
    private PrintStream genotyperDebugStream;

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
     * Calculate the BQD model outputs to the likelihoods array.
     * This method handles splitting the model by strand and selecting the best scoring parameters across the two for return in the likelihoods array.
     *
     * BQD needs to see reads that have been disqualified in {@link org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods#filterPoorlyModeledEvidence(ToDoubleFunction)} as
     * well as reads that only overlap the variant in question in their low quality ends. Reads in the former category do not have their hmm scores
     * accounted for in the genotyping model, whereas reads in the later category do. All reads, (disqualified, low quality ends, and all others)
     * are sorted by the cycle-count of the SNP being genotyped, and the average base qualities for the SNP base are computed across partitions
     * of the reads in aggregate.
     *
     * NOTES:
     * - The model will not handle indel alleles
     * - The model currently does not support mixed-allele mode (modifying 0/1 GTs in addition to 0/0 GTs)
     *
     * @param sampleLikelihoods allele liklihoods containing data for reads
     * @param strandForward list of reads in the forwards orientation overlapping the site
     * @param strandReverse list of reads in the reverse orientation overlapping the site
     * @param paddedReference reference bases (with padding) used for calculating homopolymer adjustemnt
     * @param offsetForRefIntoEvent offset of the variant into the reference event
     * @param calculators likelihoods calculators object pre-filled with scores
     * @return An array corresponding to the likelihoods array score for BQD, with Double.NEGATIVE_INFINITY filling all mixed allele/indel allelse
     */
    public <A extends Allele> double[] calculateBQDLikelihoods(final LikelihoodMatrix<GATKRead, A> sampleLikelihoods,
                                                               final List<DRAGENGenotypesModel.DragenReadContainer> strandForward,
                                                               final List<DRAGENGenotypesModel.DragenReadContainer> strandReverse,
                                                               final byte[] paddedReference,
                                                               final int offsetForRefIntoEvent,
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
                // TODO super validate this
                byte baseOfErrorAllele = sampleLikelihoods.getAllele(errorAlleleIndex).getBases()[0];

                double forwardHomopolymerAdjustment = FRDBQDUtils.computeForwardHomopolymerAdjustment(paddedReference, offsetForRefIntoEvent, baseOfErrorAllele);
                double reverseHomopolymerAdjustment = FRDBQDUtils.computeReverseHomopolymerAdjustment(paddedReference, offsetForRefIntoEvent, baseOfErrorAllele);

                // This selects the index for the reads page in the table corresponding to the likelihoods of the read given allele frequency of 1
                final int offsetForReadLikelihoodGivenAlleleIndex = alleleDataSize * errorAlleleIndex + readCount;

                // BQD scores by strand
                double minScoreFoundForwardsStrand = computeBQDModelForStrandData(strandForward, forwardHomopolymerAdjustment, readLikelihoodsForGT, offsetForReadLikelihoodGivenAlleleIndex, true, errorAlleleIndex);
                double minScoreFoundReverseStrand = computeBQDModelForStrandData(strandReverse, reverseHomopolymerAdjustment, readLikelihoodsForGT, offsetForReadLikelihoodGivenAlleleIndex, false, errorAlleleIndex);

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
     * This method works by combining the computed genotype scores for reads with the raw allele likelihoods scores for each evidence
     *
     * @param positionSortedReads  Reads pairs objects (Pair<Pair<read,readBaseOffset>, sampleReadIndex>) objects sorted in the correct order for partitioning.
     *                             This means that the "error" reads in the partition are sorted by read cycle first in the provided list
     * @param homopolymerAdjustment  Penalty to be applied to reads based on the homopolymer run (this should be precomputed for the ref site in quesiton)
     * @param readLikelihoodsForGT  The array corresponding to the log_10 genotype scores for the genotype in question
     * @param offsetForReadLikelihoodGivenAlleleIndex
     * @return phred scale liklihood for a BQD error mode for reads in the given direction according to the offsets requested
     */
    private double computeBQDModelForStrandData(final List<DRAGENGenotypesModel.DragenReadContainer> positionSortedReads,
                                                final double homopolymerAdjustment, final double[] readLikelihoodsForGT,
                                                final int offsetForReadLikelihoodGivenAlleleIndex, final boolean forwards, final int errorAlleleIndex) {
        if (positionSortedReads.isEmpty()) {
            return 0.0; // TODO check up on this
        }
        if (genotyperDebugStream != null) {
            genotyperDebugStream.println("errorAllele index: " + errorAlleleIndex + " theta: " + (forwards ? "1" : "2") + " homopolymerAdjustment: " + homopolymerAdjustment);
        }

        // Forwards strand tables (all in phred space for the sake of conveient debugging with provided scripts)
        final double[] cumulative_p_R_for_E = new double[positionSortedReads.size() + 1];
        final double[] cumulative_mean_base_quality_phred_adjusted = new double[positionSortedReads.size() + 1];
        final double[] cumulative_P_GT = new double[positionSortedReads.size() + 1];

        double totalBaseQuality = 0;
        int baseQualityDenominator = 0; // We track this seperately because not every read overlaps the SNP in quesiton due to padding.
        // Iterate over the reads and populate the cumulative arrays
        for (int i = 1; i < cumulative_p_R_for_E.length; i++) {
            final DRAGENGenotypesModel.DragenReadContainer container = positionSortedReads.get(i - 1);
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

////            System.out.println("read index:"+readIndex+" ihomGT:"+-10*homozygousGenotypeContribution+"  errorAlleleContribution:"+-10*errorAlleleContribution + " difference in phred: "+ -10*(errorAlleleContribution - homozygousGenotypeContribution));

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
        if (genotyperDebugStream != null) {
            genotyperDebugStream.println(String.format("theta=%d n%d=%2d, best_phred_score =%5.2f q_mean=%5.2f, alpha=%4.2f, Ph(E)=%4.2f;  ", forwards ? 1 : 0,
                    (forwards ? 1 : 2), nIndexUsed, minScoreFound, cumulative_mean_base_quality_phred_adjusted[nIndexUsed],
                    BQD_FIXED_DEFAULT_ALPHA, cumulative_p_R_for_E[nIndexUsed]));
        }

        //TODO this will be where i put my debuglogging (if i had any)
        return minScoreFound;
    }

    /**
     * Calculate the BQD model outputs to the likelihoods array.
     *
     * This method is responsible for computing critical phred-mapping quality adjustments for the entire pool of reads (Disqualified reads,
     * reads only overlapping in low quality ends, and otherwise) and selecting true-allele/error-allele combinations as well as strand model
     * combinations (all forward reads/ all reverse reads/ all reads) and calling {@link #computeFRDModelForStrandData} for each of these
     * combinations selecting the best scoring columns in the final likelihoods array output.
     *
     * Like BQD this model genotypes with all reads that overlap the site in either their accepted bases or low quality ends, but it does
     * not include disqualified reads for genotyping. All reads are used for computing the critical values for the mapping quality cutoffs.
     *
     * NOTES:
     * - The model currently does not support mixed-allele mode (modifying 0/1 GTs in addition to 0/0 GTs)
     * - The model will not treat symbolic alleles specially, always treating them as indels. This might or might not be the best way to handle them.
     *
     * @param sampleLikelihoods the likelihoods object with allele liklihoods for the reads to be genotyped
     * @param ploidyModelLikelihoods standard genotyping model allele likelihoods (to be used for maxEffectiveDepthAdjustment)
     * @param readContainers reads (both forwards and reverse orientation as well as disqualified reads) overlapping the site in question
     * @param snipAprioriHet prior for heterozygus SNP allele
     * @param indelAprioriHet prior for heterozygus indel alleles based on the STRE tables if present
     * @param maxEffectiveDepthForHetAdjustment maxEffectiveDepthAdjustment used to reduce the effect of FRD at high depth sites (0 means no adjustment)
     * @param calculators DRAGENlikelihoodsCalculator object to manage GT array math
     * @return a likelihoods array corrsponding to the log10 likelihoods scores for the best combination of model parameters for each Genotype (Double.NEGATIVE_INFINITY for Genotypes not considered)
     */
    public <A extends Allele> double[] calculateFRDLikelihoods(final LikelihoodMatrix<GATKRead, A> sampleLikelihoods, final double[] ploidyModelLikelihoods,
                                                               final List<DRAGENGenotypesModel.DragenReadContainer> readContainers,
                                                               final double snipAprioriHet, final double indelAprioriHet, final int maxEffectiveDepthForHetAdjustment,
                                                               final GenotypeLikelihoodCalculators calculators) {
        Utils.validate(sampleLikelihoods == cachedLikelihoodsObject, "There was a mismatch between the sample stored by the genotyper and the one requested for BQD, this will result in invalid genotyping");
        final double[] outputArray = new double[genotypeCount];
        Arrays.fill(outputArray, Double.NEGATIVE_INFINITY);

        final Allele refAllele = sampleLikelihoods.getAllele(0);

        //Determine the size of an allele page for the readsLikelihoodsByAlleleFrequency table
        final int readCount = sampleLikelihoods.evidenceCount();
        final int alleleDataSize = readCount * (ploidy + 1);

        for(int fAlleleIndex = 0; fAlleleIndex < sampleLikelihoods.numberOfAlleles(); fAlleleIndex++) {
            // ignore symbolic alleles
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
            if (genotyperDebugStream != null) {
                genotyperDebugStream.println("fIndex: " + fAlleleIndex + " criticalValues: \n" + criticalThresholds.stream().map(d -> Double.toString(d)).collect(Collectors.joining("\n")));
            }
            // iterate over all of the homozygous genotypes for the given allele
            for(int gtAlleleIndex = 0; gtAlleleIndex < sampleLikelihoods.numberOfAlleles(); gtAlleleIndex++) {
                // Skip over the allele corresponding to the "foreign" allele
                if (gtAlleleIndex == fAlleleIndex) {
                    continue;
                }
                // For right now we allow symbolic alleles, but this might be subject to change
//                if (sampleLikelihoods.getAllele(fAlleleIndex).isSymbolic() ) {
//                    continue;
//                }

                //This is crufty, it just so happens that the index of the homozygous genotype corresponds to the maximum genotype count per field.
                //This should be pulled off as a calculator in some genotyping class.
                final int indexForGT = calculators.genotypeCount(ploidy, gtAlleleIndex + 1) - 1;
                double[] readLikelihoodsForGT = readLikelihoodsByGenotypeIndex[indexForGT];

                // TODO restore the critical thresholds
                if (genotyperDebugStream != null) {
                    genotyperDebugStream.println("indexForGT "+indexForGT+ " ooffsetForReadLikelihoodGivenAlleleIndex ="+offsetForReadLikelihoodGivenAlleleIndex);
                    genotyperDebugStream.println("\nForwards Strands: ");
                }
                final double[] maxLog10FForwardsStrand = computeFRDModelForStrandData(readContainers, c -> !c.isReverseStrand() , offsetForReadLikelihoodGivenAlleleIndex, readLikelihoodsForGT, criticalThresholdsSorted);
                if (genotyperDebugStream != null) {  genotyperDebugStream.println("\nReverse Strands: ");}
                final double[] maxLog10FReverseStrand = computeFRDModelForStrandData(readContainers, c -> c.isReverseStrand(), offsetForReadLikelihoodGivenAlleleIndex, readLikelihoodsForGT, criticalThresholdsSorted);
                if (genotyperDebugStream != null) {  genotyperDebugStream.println("\nBoth Strands: ");}
                final double[] maxLog10FBothStrands = computeFRDModelForStrandData(readContainers, c -> true, offsetForReadLikelihoodGivenAlleleIndex, readLikelihoodsForGT, criticalThresholdsSorted);

                if (genotyperDebugStream != null) {
                    genotyperDebugStream.println("gtAlleleIndex : "+gtAlleleIndex+ " fAlleleIndex: "+fAlleleIndex +" forwards: "+maxLog10FForwardsStrand+" reverse: "+maxLog10FReverseStrand+" both: "+maxLog10FBothStrands);

                }
                double[] localBestModel = maxLog10FForwardsStrand;
                if (localBestModel[0] < maxLog10FReverseStrand[0]) {
                    localBestModel = maxLog10FReverseStrand;
                }
                if (localBestModel[0] < maxLog10FBothStrands[0]) {
                    localBestModel = maxLog10FBothStrands;
                }

                // Handle max effective depth adjustment if specified
                if (maxEffectiveDepthForHetAdjustment > 0) {
                    // Use the index corresponding the mixture of F and
                    double localBestModelScore = localBestModel[0] - localBestModel[1];
                    int closestGTAlleleIndex = allelesToIndex(gtAlleleIndex, fAlleleIndex);
                    double log10LikelihoodsForPloyidyModel = ploidyModelLikelihoods[closestGTAlleleIndex] - MathUtils.log10(2);
                    int depthForGenotyping = sampleLikelihoods.evidenceCount();
                    double adjustedBestModel = log10LikelihoodsForPloyidyModel + ((localBestModelScore - log10LikelihoodsForPloyidyModel)
                            * ((Math.min(depthForGenotyping, maxEffectiveDepthForHetAdjustment) * 1.0) / depthForGenotyping));
                    outputArray[indexForGT] = Math.max(outputArray[indexForGT], adjustedBestModel + localBestModel[1]);

                    if (genotyperDebugStream != null) {
                        genotyperDebugStream.println("best FRD likelihoods: "+localBestModelScore+" P(F) score used: "+localBestModel[1]+"  use MaxEffectiveDepth: "+maxEffectiveDepthForHetAdjustment);
                        genotyperDebugStream.println("Using array index "+closestGTAlleleIndex+" for mixture gt with likelihood of "+log10LikelihoodsForPloyidyModel+" adjusted based on depth: "+depthForGenotyping);
                        genotyperDebugStream.println("p_rG_adj : "+adjustedBestModel);
                    }
                } else {
                    outputArray[indexForGT] = Math.max(outputArray[indexForGT], localBestModel[0]);
                }

            }

        }


        return outputArray;
    }


    /**
     * @param positionSortedReads read containers to use for genotyping
     * @param predicate predicate used to select the correct orientation combination for reads when genotyping
     * @param offsetForReadLikelihoodGivenAlleleIndex offset corresponding to the Error Allele in the reads likelihoods object array
     * @param readLikelihoodsForGT reads likelihoods for Genotype array table
     * @param criticalThresholdsSorted critical thresholds to use for this orientation combination
     * @return two doubles, index 0 is the frd score and the second is log p(F()) score used to adjust the score
     */
    private double[] computeFRDModelForStrandData(final List<DRAGENGenotypesModel.DragenReadContainer> positionSortedReads, final Predicate<DRAGENGenotypesModel.DragenReadContainer> predicate,
                                                  final int offsetForReadLikelihoodGivenAlleleIndex, final double[] readLikelihoodsForGT, final List<Double> criticalThresholdsSorted) {
        if (positionSortedReads.isEmpty()) {
            return new double[]{Double.NEGATIVE_INFINITY, 0};
        }

        int counter = 0;
        double maxLpspi = Double.NEGATIVE_INFINITY;
        double lpfApplied = 0;

        for (final Double logProbFAllele : criticalThresholdsSorted) {
            double f_ratio = 0.0;
            double f_denom = 0.0;
            double localMaxLpspi = Double.NEGATIVE_INFINITY;

            // iterate over the
            for (final DRAGENGenotypesModel.DragenReadContainer container : positionSortedReads) {
                // Ignore reads that were disqualified by the HMM
                if (container.wasFilteredByHMM()) {
                    continue;
                }

                final int readIndex = container.getIndexInLikelihoodsObject();

                // Keep track of the aggregate support for the foreign allele
                if (predicate.test(container)) {
                    // Only include reads with mapping quality adjustment < the critical threshold being used (i.e. exclude reads with MQ > than the threshold)
                    double LPd_r_F = container.getPhredPFValue() + 0.0000001 <= logProbFAllele ?
                            Double.NEGATIVE_INFINITY :
                            readAlleleLikelihoodByAlleleCount[offsetForReadLikelihoodGivenAlleleIndex + readIndex];
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
            for (final DRAGENGenotypesModel.DragenReadContainer container : positionSortedReads) {
                // Ignore reads that were disqualified by the HMM
                if (container.wasFilteredByHMM()) {
                    continue;
                }
                final int readIndex = container.getIndexInLikelihoodsObject();

                double lp_r_GT = readLikelihoodsForGT[readIndex] - MathUtils.log10(2);

                // COMPUTE THE MODEL FOR THE STRAND IN QUESTION
                if (predicate.test(container)) {
                    double LPd_r_F = container.getPhredPFValue() + 0.0000001 <= logProbFAllele ?
                            Double.NEGATIVE_INFINITY :
                            readAlleleLikelihoodByAlleleCount[offsetForReadLikelihoodGivenAlleleIndex + readIndex];

                    LP_R_GF += MathUtils.approximateLog10SumLog10(log10beta + LPd_r_F, log10betaDelta + lp_r_GT);
                } else {
                    LP_R_GF += lp_r_GT;
                }
            }
            // Allele prior for error allele, plus posterior for foreign event, plus model posterior
            double LPsi = logProbFAllele + LP_R_GF; // NOTE unlike DRAGEN we apply the prior to the combined likelihoods array after the fact so gtAllelePrior is not included at this stage
            localMaxLpspi = Math.max(localMaxLpspi, LPsi);

            if (genotyperDebugStream != null) {
                genotyperDebugStream.println("beta: "+beta+" localMaxLpspi: " + localMaxLpspi + " for lpf: "+logProbFAllele+" with LP_R_GF: "+LP_R_GF+" index: "+counter++);
            }
            if (localMaxLpspi > maxLpspi) {
                maxLpspi = Math.max(maxLpspi, localMaxLpspi);
                lpfApplied = logProbFAllele;
            }
        }

        //TODO javaize this
        // TODO soon should not need to use the LPF applied here...
        return new double[]{maxLpspi, lpfApplied};
    }


    // TODO: for reviewer... this code is meant to handle the performance regression brought about by FRD, essentially in DRAGEN for every
    // TODO  combination of strandedness and mapping quality a computation is done, this is vaguly unnessisary, essentially this leads to
    // TODO  a bias towards strandedness/allele combinations that have very few reads supporting them by viritue of the fact that a different
    // TODO  strand/allele combination had at least one low mapping quality read. I have reverted the change right now and consequently this genotyping
    // TODO  is taking somewhere in the order of ~5-6% runtime on the profiler whereas otherwise it could correspond to much less at the expense
    // TODO  of not matching DRAGEN properly.
    // helper method to populate the reads containers properly with their critical values and store them in the provided set
    // NOTE: this has the side effect of setting the DragenReadContainer setPhredPFValue() values for the reads for the given set of alleles
    private void computeCriticalValues(final Set<Double> criticalThresholdsForwards, final Set<Double> criticalThresholdsReverse, final Set<Double> criticalThresholdsTotal, final List<DRAGENGenotypesModel.DragenReadContainer> container, final double log10MapqPriorAdjustment) {
        for (int i = 0; i < container.size(); i++) {
            final DRAGENGenotypesModel.DragenReadContainer readContainer = container.get(i);
            final double log10CriticalValue = readContainer.getPhredScaledMappingQuality() / -10.0 + log10MapqPriorAdjustment;
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

    // Add debug out stream to to output genotyper debug information if necessary
    public void addGenotyperDebugOutputStream(final PrintStream genotyperDebugStream) {
        this.genotyperDebugStream = genotyperDebugStream;
    }
}
