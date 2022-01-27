package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerGenotypingDebugger;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.GATKRead;

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
    static final double BQD_FIXED_ERROR_RATE = 0.5;

    // PhredScaled adjustment applied to the BQD score (this controls the weight of the base quality prior term in the BQD calculation)
    static final double PHRED_SCALED_ADJUSTMENT_FOR_BQ_SCORE = 2.5;

    private final double cachedLog10ErrorRate;
    private final double cachedLog10NonErrorRate;

    /**
     * Creates a new calculator providing its ploidy and number of genotyping alleles.
     */
    protected GenotypeLikelihoodCalculatorDRAGEN(final int ploidy, final int alleleCount,
                                                 final GenotypeAlleleCounts[][] genotypeTableByPloidy) {
        super(ploidy, alleleCount, genotypeTableByPloidy);
        Utils.validateArg(ploidy > 0, () -> "ploidy must be at least 1 but was " + ploidy);
        // The number of possible components is limited by distinct allele count and ploidy.
        cachedLog10ErrorRate = Math.log10(BQD_FIXED_ERROR_RATE);
        cachedLog10NonErrorRate = Math.log10(1 - BQD_FIXED_ERROR_RATE);
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
     * @param sampleLikelihoods allele likelihoods containing data for reads
     * @param strandForward list of reads in the forwards orientation overlapping the site
     * @param strandReverse list of reads in the reverse orientation overlapping the site
     * @param paddedReference reference bases (with padding) used for calculating homopolymer adjustemnt
     * @param offsetForRefIntoEvent offset of the variant into the reference event
     * @return An array corresponding to the likelihoods array score for BQD, with Double.NEGATIVE_INFINITY filling all mixed allele/indel allelse
     */
    public <A extends Allele> double[] calculateBQDLikelihoods(final LikelihoodMatrix<GATKRead, A> sampleLikelihoods,
                                                               final List<DRAGENGenotypesModel.DragenReadContainer> strandForward,
                                                               final List<DRAGENGenotypesModel.DragenReadContainer> strandReverse,
                                                               final byte[] paddedReference,
                                                               final int offsetForRefIntoEvent) {
        final double[] outputArray = new double[genotypeCount];
        Arrays.fill(outputArray, Double.NEGATIVE_INFINITY);

        final Allele refAllele = sampleLikelihoods.getAllele(0);

        for (int gtAlleleIndex = 0; gtAlleleIndex < sampleLikelihoods.numberOfAlleles(); gtAlleleIndex++) {
            final int indexForGT = genotypeIndexCalculator.alleleCountsToIndex(gtAlleleIndex, ploidy);

            for (int errorAlleleIndex = 0; errorAlleleIndex < sampleLikelihoods.numberOfAlleles(); errorAlleleIndex++) {
                // We only want to make calls on SNPs for now
                if (sampleLikelihoods.getAllele(gtAlleleIndex) == sampleLikelihoods.getAllele(errorAlleleIndex) ||
                        sampleLikelihoods.getAllele(gtAlleleIndex).length() != refAllele.length() ||
                        sampleLikelihoods.getAllele(errorAlleleIndex).length() != refAllele.length()) {
                    continue;
                }
                // TODO super validate this
                final byte baseOfErrorAllele = sampleLikelihoods.getAllele(errorAlleleIndex).getBases()[0];

                final double forwardHomopolymerAdjustment = FRDBQDUtils.computeForwardHomopolymerAdjustment(paddedReference, offsetForRefIntoEvent, baseOfErrorAllele);
                final double reverseHomopolymerAdjustment = FRDBQDUtils.computeReverseHomopolymerAdjustment(paddedReference, offsetForRefIntoEvent, baseOfErrorAllele);

                // BQD scores by strand
                final double minScoreFoundForwardsStrand = computeBQDModelForStrandData(sampleLikelihoods, strandForward, forwardHomopolymerAdjustment, true, gtAlleleIndex, errorAlleleIndex);
                final double minScoreFoundReverseStrand = computeBQDModelForStrandData(sampleLikelihoods, strandReverse, reverseHomopolymerAdjustment, false, gtAlleleIndex, errorAlleleIndex);

                final double modelScoreInLog10 = (minScoreFoundForwardsStrand + minScoreFoundReverseStrand) * -0.1;
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
     * @param sampleLikelihoods allele likelihoods containing data for reads
     * @param positionSortedReads  Reads pairs objects (Pair<Pair<read,readBaseOffset>, sampleReadIndex>) objects sorted in the correct order for partitioning.
     *                             This means that the "error" reads in the partition are sorted by read cycle first in the provided list
     * @param homopolymerAdjustment  Penalty to be applied to reads based on the homopolymer run (this should be precomputed for the ref site in quesiton)
     * @return phred scale likelihood for a BQD error mode for reads in the given direction according to the offsets requested
     */
    private <A extends Allele> double computeBQDModelForStrandData(final LikelihoodMatrix<GATKRead, A> sampleLikelihoods,
                                                                   final List<DRAGENGenotypesModel.DragenReadContainer> positionSortedReads,
                                                                   final double homopolymerAdjustment,
                                                                   final boolean forwards, final int homozygousAlleleIndex, final int errorAlleleIndex) {
        // If no reads are found for a particular strand direction return no adjusted likelihoods score for those (non-existent) reads
        if (positionSortedReads.isEmpty()) {
            return 0.0;
        }
        if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
            HaplotypeCallerGenotypingDebugger.println("errorAllele index: " + errorAlleleIndex + " theta: " + (forwards ? "1" : "2") + " homopolymerAdjustment: " + homopolymerAdjustment);
        }

        // Forwards strand tables (all in phred space for the sake of conveient debugging with provided scripts)
        final int evidenceSize = positionSortedReads.size();
        final double[] cumulativeProbReadForErrorAllele = new double[evidenceSize + 1];
        final double[] cumulativeMeanBaseQualityPhredAdjusted = new double[evidenceSize + 1];
        final double[] cumulativeProbGenotype = new double[evidenceSize + 1];

        double totalBaseQuality = 0;
        int baseQualityDenominator = 0; // We track this separately because not every read overlaps the SNP in question due to padding.
        // Iterate over the reads and populate the cumulative arrays
        for (int i = 1; i < cumulativeProbReadForErrorAllele.length; i++) {
            final DRAGENGenotypesModel.DragenReadContainer container = positionSortedReads.get(i - 1);
            final int readIndex = container.getIndexInLikelihoodsObject();

            // Retrieve the homozygous genotype score and the error allele scores for the read in question
            final double homozygousGenotypeContribution;
            final double errorAlleleContribution;
            if (readIndex != -1) {
                homozygousGenotypeContribution = sampleLikelihoods.get(homozygousAlleleIndex, readIndex);
                errorAlleleContribution = sampleLikelihoods.get(errorAlleleIndex, readIndex);
            } else {
                // If read index == -1 then we are evaluating a read that was rejected by the HMM and therefore doesn't have genotype scores
                homozygousGenotypeContribution = 0;
                errorAlleleContribution = 0;
            }

            // Populate the error probability array in phred space
            // Calculation: Alpha * P(r|E_allele) + (1 - Alpha) * P(r | G_homozygousGT))
            double phredContributionForRead = (homozygousGenotypeContribution==0 && errorAlleleContribution==0) ? 0 : -10 *
                    MathUtils.approximateLog10SumLog10(errorAlleleContribution + cachedLog10ErrorRate,
                                                       homozygousGenotypeContribution + cachedLog10NonErrorRate);
            cumulativeProbReadForErrorAllele[i] = cumulativeProbReadForErrorAllele[i-1] + phredContributionForRead;

            // Populate the cumulative genotype contribution array with the score for this read
            // Calculation: (P(r | G_A1) + P(r | G_A2)) / 2
            cumulativeProbGenotype[i] = cumulativeProbGenotype[i - 1] + -10 * homozygousGenotypeContribution;

            // Populate the mean base quality array
            if (container.hasValidBaseQuality()) {
                totalBaseQuality += container.getBaseQuality();
                baseQualityDenominator++;
            }
            cumulativeMeanBaseQualityPhredAdjusted[i] = Math.max(0,
                    ((totalBaseQuality / (baseQualityDenominator == 0 ? 1 : baseQualityDenominator)) * PHRED_SCALED_ADJUSTMENT_FOR_BQ_SCORE) - homopolymerAdjustment);
        }

        // Now we find the best partitioning N for the forwards evaluation of the data
        double minScoreFound = Double.POSITIVE_INFINITY;
        int nIndexUsed = 0;
        for (int n = 0; n < cumulativeMeanBaseQualityPhredAdjusted.length; n++) {
            final double bqdScore = cumulativeMeanBaseQualityPhredAdjusted[n] + cumulativeProbReadForErrorAllele[n] + (cumulativeProbGenotype[cumulativeProbGenotype.length-1] - cumulativeProbGenotype[n]);
            if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
                HaplotypeCallerGenotypingDebugger.println(String.format("n=%d: %.2f, cum_phred_bq=%.2f, cum_prob_r_Error=%.2f, prob_G_remaining=%.2f",
                        n, bqdScore, cumulativeMeanBaseQualityPhredAdjusted[n], cumulativeProbReadForErrorAllele[n],
                        (cumulativeProbGenotype[cumulativeProbGenotype.length - 1] - cumulativeProbGenotype[n])));
            }
            if (minScoreFound > bqdScore) {
                minScoreFound = bqdScore;
                nIndexUsed = n;
            }
        }

        // Debug output for the genotyper to see into the calculation itself
        if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
            HaplotypeCallerGenotypingDebugger.println(String.format("theta=%d n%d=%2d, best_phred_score =%5.2f q_mean=%5.2f, alpha=%4.2f, Ph(E)=%4.2f;  ", forwards ? 1 : 0,
                    (forwards ? 1 : 2), nIndexUsed, minScoreFound, cumulativeMeanBaseQualityPhredAdjusted[nIndexUsed],
                    BQD_FIXED_ERROR_RATE, cumulativeProbReadForErrorAllele[nIndexUsed]));
        }

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
     * @param sampleLikelihoods the likelihoods object with allele likelihoods for the reads to be genotyped
     * @param ploidyModelLikelihoods standard genotyping model allele likelihoods (to be used for maxEffectiveDepthAdjustment)
     * @param readContainers reads (both forwards and reverse orientation as well as disqualified reads) overlapping the site in question
     * @param snipAprioriHet prior for heterozygus SNP allele
     * @param indelAprioriHet prior for heterozygus indel alleles based on the STRE tables if present
     * @param maxEffectiveDepthForHetAdjustment maxEffectiveDepthAdjustment used to reduce the effect of FRD at high depth sites (0 means no adjustment)
     * @return a likelihoods array corrsponding to the log10 likelihoods scores for the best combination of model parameters for each Genotype (Double.NEGATIVE_INFINITY for Genotypes not considered)
     */
    public <A extends Allele> double[] calculateFRDLikelihoods(final LikelihoodMatrix<GATKRead, A> sampleLikelihoods, final double[] ploidyModelLikelihoods,
                                                               final List<DRAGENGenotypesModel.DragenReadContainer> readContainers,
                                                               final double snipAprioriHet, final double indelAprioriHet, final int maxEffectiveDepthForHetAdjustment) {
        final double[] outputArray = new double[genotypeCount];
        Arrays.fill(outputArray, Double.NEGATIVE_INFINITY);

        final Allele refAllele = sampleLikelihoods.getAllele(0);

        for (int fAlleleIndex = 0; fAlleleIndex < sampleLikelihoods.numberOfAlleles(); fAlleleIndex++) {
            // ignore symbolic alleles
            final boolean isIndel = sampleLikelihoods.getAllele(fAlleleIndex).length() != refAllele.length();

            // Here we generate a set of the critical log10(P(F)) values that we will iterate over
            final FRDCriticalThresholds thresholds = computeCriticalValues(readContainers, fAlleleIndex == 0 ? 0 : (isIndel? indelAprioriHet : snipAprioriHet) * -0.1); // simplified in line with DRAGEN, uses 1 alleledist for both snp and indels

            if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
                HaplotypeCallerGenotypingDebugger.println("fIndex: " + fAlleleIndex + " criticalValues: \n" + thresholds.getCriticalThresholdsTotal().stream().map(d -> Double.toString(d)).collect(Collectors.joining("\n")));
            }
            // iterate over all of the homozygous genotypes for the given allele
            for (int gtAlleleIndex = 0; gtAlleleIndex < sampleLikelihoods.numberOfAlleles(); gtAlleleIndex++) {
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
                final int indexForGT = genotypeIndexCalculator.alleleCountsToIndex(gtAlleleIndex, ploidy);

                // TODO restore the critical thresholds
                if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
                    HaplotypeCallerGenotypingDebugger.println("indexForGT "+indexForGT);
                    HaplotypeCallerGenotypingDebugger.println("\nForwards Strands: ");
                }
                final double[] maxLog10FForwardsStrand = computeFRDModelForStrandData(sampleLikelihoods, gtAlleleIndex, fAlleleIndex, readContainers,
                        c -> !c.isReverseStrand() , thresholds.getCriticalThresholdsTotal());
                if (HaplotypeCallerGenotypingDebugger.isEnabled()) {  HaplotypeCallerGenotypingDebugger.println("\nReverse Strands: ");}
                final double[] maxLog10FReverseStrand = computeFRDModelForStrandData(sampleLikelihoods, gtAlleleIndex, fAlleleIndex, readContainers,
                        c -> c.isReverseStrand(), thresholds.getCriticalThresholdsTotal());
                if (HaplotypeCallerGenotypingDebugger.isEnabled()) {  HaplotypeCallerGenotypingDebugger.println("\nBoth Strands: ");}
                final double[] maxLog10FBothStrands = computeFRDModelForStrandData(sampleLikelihoods, gtAlleleIndex, fAlleleIndex, readContainers,
                        c -> true, thresholds.getCriticalThresholdsTotal());

                if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
                    HaplotypeCallerGenotypingDebugger.println("gtAlleleIndex : "+gtAlleleIndex+ " fAlleleIndex: "+fAlleleIndex +" forwards: "+maxLog10FForwardsStrand+" reverse: "+maxLog10FReverseStrand+" both: "+maxLog10FBothStrands);
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
                    final double localBestModelScore = localBestModel[0] - localBestModel[1];
                    final int closestGTAlleleIndex = allelesToIndex(gtAlleleIndex, fAlleleIndex);
                    final double log10LikelihoodsForPloyidyModel = ploidyModelLikelihoods[closestGTAlleleIndex] - -MathUtils.LOG10_ONE_HALF;
                    final int depthForGenotyping = sampleLikelihoods.evidenceCount();
                    final double adjustedBestModel = log10LikelihoodsForPloyidyModel + ((localBestModelScore - log10LikelihoodsForPloyidyModel)
                            * ((Math.min(depthForGenotyping, maxEffectiveDepthForHetAdjustment) * 1.0) / depthForGenotyping));
                    outputArray[indexForGT] = Math.max(outputArray[indexForGT], adjustedBestModel + localBestModel[1]);

                    if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
                        HaplotypeCallerGenotypingDebugger.println("best FRD likelihoods: "+localBestModelScore+" P(F) score used: "+localBestModel[1]+"  use MaxEffectiveDepth: "+maxEffectiveDepthForHetAdjustment);
                        HaplotypeCallerGenotypingDebugger.println("Using array index "+closestGTAlleleIndex+" for mixture gt with likelihood of "+log10LikelihoodsForPloyidyModel+" adjusted based on depth: "+depthForGenotyping);
                        HaplotypeCallerGenotypingDebugger.println("p_rG_adj : "+adjustedBestModel);
                    }
                } else {
                    outputArray[indexForGT] = Math.max(outputArray[indexForGT], localBestModel[0]);
                }

            }

        }


        return outputArray;
    }


    /**
     * @param sampleLikelihoods the likelihoods object with allele likelihoods for the reads to be genotyped
     * @param homozygousAlleleIndex index of allele in homzygous genotype whose likelihood is to be adjusted
     * @param fAlleleIndex index of foreign allele within likelihoods matrix
     * @param positionSortedReads read containers to use for genotyping
     * @param predicate predicate used to select the correct orientation combination for reads when genotyping
     * @param criticalThresholdsSorted critical thresholds to use for this orientation combination
     * @return two doubles, index 0 is the frd score and the second is log p(F()) score used to adjust the score
     */
    private <A extends Allele> double[] computeFRDModelForStrandData(final LikelihoodMatrix<GATKRead, A> sampleLikelihoods,
                                                                     final int homozygousAlleleIndex, final int fAlleleIndex,
                                                                     final List<DRAGENGenotypesModel.DragenReadContainer> positionSortedReads,
                                                                     final Predicate<DRAGENGenotypesModel.DragenReadContainer> predicate,
                                                                     final List<Double> criticalThresholdsSorted) {
        if (positionSortedReads.isEmpty()) {
            return new double[]{Double.NEGATIVE_INFINITY, 0};
        }

        int counter = 0;
        double maxLpspi = Double.NEGATIVE_INFINITY;
        double lpfApplied = 0;

        for (final Double logProbFAllele : criticalThresholdsSorted) {
            double fAlleleProbRatio = 0.0;
            double fAlleleProbDenom = 0.0;
            double localMaxLpspi = Double.NEGATIVE_INFINITY;

            // iterate over the reads to compute the foreign allele alpha to use for genotyping with FRD
            for (final DRAGENGenotypesModel.DragenReadContainer container : positionSortedReads) {
                // Ignore reads that were disqualified by the HMM
                if (container.wasFilteredByHMM()) {
                    continue;
                }

                final int readIndex = container.getIndexInLikelihoodsObject();

                // Keep track of the aggregate support for the foreign allele
                if (predicate.test(container)) {
                    // Only include reads with mapping quality adjustment < the critical threshold being used (i.e. exclude reads with MQ > than the threshold)
                    final double LPd_r_F = container.getPhredPFValue() + 0.0000001 <= logProbFAllele ?
                            Double.NEGATIVE_INFINITY :
                            sampleLikelihoods.get(fAlleleIndex, readIndex);
                    final double lp_r_GT = sampleLikelihoods.get(homozygousAlleleIndex, readIndex);

                    fAlleleProbRatio += Math.pow(10, LPd_r_F - MathUtils.approximateLog10SumLog10(LPd_r_F, lp_r_GT));
                    fAlleleProbDenom++;
                }
            }

            // Don't learn the beta but approximate it based on the read support for the alt
            final double foreignAlleleLikelihood = Math.min(fAlleleProbRatio / fAlleleProbDenom, 0.5);
            final double log10ForeignAlleleLikelihood = Math.log10(foreignAlleleLikelihood);
            final double log10NotForeignAlleleLikelihood = Math.log10(1.0 - foreignAlleleLikelihood);
            double cumulativeLog10LikelihoodOfForeignReadHypothesis = 0.0; // LP_R_GF

            // iterate over the containers again using the approximated beta constraint
            for (final DRAGENGenotypesModel.DragenReadContainer container : positionSortedReads) {
                // Ignore reads that were disqualified by the HMM
                if (container.wasFilteredByHMM()) {
                    continue;
                }
                final int readIndex = container.getIndexInLikelihoodsObject();

                final double log10LikelihoodReadForGenotype = sampleLikelihoods.get(homozygousAlleleIndex, readIndex);

                // COMPUTE THE MODEL FOR THE STRAND IN QUESTION
                if (predicate.test(container)) {
                    final double log10LikelihoodOfForeignAlleleGivenLPFCutoff = container.getPhredPFValue() + 0.0000001 <= logProbFAllele ?
                            Double.NEGATIVE_INFINITY :
                            sampleLikelihoods.get(fAlleleIndex, readIndex);;

                    cumulativeLog10LikelihoodOfForeignReadHypothesis += MathUtils.approximateLog10SumLog10(log10ForeignAlleleLikelihood + log10LikelihoodOfForeignAlleleGivenLPFCutoff, log10NotForeignAlleleLikelihood + log10LikelihoodReadForGenotype);
                } else {
                    cumulativeLog10LikelihoodOfForeignReadHypothesis += log10LikelihoodReadForGenotype;
                }
            }
            // Allele prior for error allele, plus posterior for foreign event, plus model posterior
            double LPsi = logProbFAllele + cumulativeLog10LikelihoodOfForeignReadHypothesis; // NOTE unlike DRAGEN we apply the prior to the combined likelihoods array after the fact so gtAllelePrior is not included at this stage
            localMaxLpspi = Math.max(localMaxLpspi, LPsi);

            if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
                HaplotypeCallerGenotypingDebugger.println("beta: "+foreignAlleleLikelihood+" localMaxLpspi: " + localMaxLpspi + " for lpf: "+logProbFAllele+" with LP_R_GF: "+cumulativeLog10LikelihoodOfForeignReadHypothesis+" index: "+counter++);
            }
            if (localMaxLpspi > maxLpspi) {
                maxLpspi = Math.max(maxLpspi, localMaxLpspi);
                lpfApplied = logProbFAllele;
            }
        }

        // TODO soon should not need to use the LPF applied here...
        return new double[]{maxLpspi, lpfApplied};
    }


    // TODO: for reviewer... this code is meant to handle the performance regression brought about by FRD. Essentially in DRAGEN for every
    // TODO  combination of strandedness and mapping quality a computation is done, this is vaguely unnecessary. Because this leads to
    // TODO  a bias towards strandedness/allele combinations that have very few reads supporting them by virtue of the fact that a different
    // TODO  strand/allele combination had at least one low mapping quality read. I have reverted the change right now and consequently this genotyping
    // TODO  is taking somewhere in the order of ~5-6% runtime on the profiler whereas otherwise it could correspond to much less at the expense
    // TODO  of not matching DRAGEN properly.
    // helper method to populate the reads containers properly with their critical values and store them in the provided set
    // NOTE: this has the side effect of setting the DragenReadContainer setPhredPFValue() values for the reads for the given set of alleles
    private FRDCriticalThresholds computeCriticalValues(final List<DRAGENGenotypesModel.DragenReadContainer> container, final double log10MapqPriorAdjustment) {
        final Set<Double> criticalThresholdsForwards = new HashSet<>();
        final Set<Double> criticalThresholdsReverse = new HashSet<>();
        final Set<Double> criticalThresholdsTotal = new HashSet<>();

        for (int i = 0; i < container.size(); i++) {
            final DRAGENGenotypesModel.DragenReadContainer readContainer = container.get(i);
            final double log10CriticalValue = readContainer.getPhredScaledMappingQuality() * -0.1 + log10MapqPriorAdjustment;
            readContainer.setPhredPFValue(log10CriticalValue);
            // Split the critical thresholds up by their applicable strands in order to avoid repeated work
//            if (readContainer.isReverseStrand()) {
//                criticalThresholdsReverse.add(log10CriticalValue);
//            } else {
//                criticalThresholdsForwards.add(log10CriticalValue);
//            }
            criticalThresholdsTotal.add(log10CriticalValue);
        }

        return new FRDCriticalThresholds(criticalThresholdsForwards, criticalThresholdsReverse, criticalThresholdsTotal);
    }


    /**
     * Helper class for storing FRD sorted and de-duplicated critical thresholds generated from reads to be accessed by subsequent calls.
     */
    private class FRDCriticalThresholds {
        private final List<Double> criticalThresholdsForwards;
        private final List<Double> criticalThresholdsReverse;
        private final List<Double> criticalThresholdsTotal;

        private FRDCriticalThresholds(final Set<Double> criticalThresholdsForwards, final Set<Double> criticalThresholdsReverse, final Set<Double> criticalThresholdsTotal) {
            this.criticalThresholdsForwards = criticalThresholdsForwards.stream().sorted(Double::compareTo).collect(Collectors.toList());
            this.criticalThresholdsReverse = criticalThresholdsReverse.stream().sorted(Double::compareTo).collect(Collectors.toList());
            this.criticalThresholdsTotal = criticalThresholdsTotal.stream().sorted(Double::compareTo).collect(Collectors.toList());
        }

        public List<Double> getCriticalThresholdsTotal() {
            return criticalThresholdsTotal;
        }


        // TODO while these are currently not being used a potential optimization to FRD would involve feeding the FRD computations only thresholds
        // TODO that are relevant/possible to be selected for a given strandedness. DRAGEN does not do this optimization and I have opted to not implement
        // TODO it here since it is complicated by the fact that critical thresholds that are not present in reads for a given strand might end up being
        // TODO selected if it is lower than a threshold where there is a change.
        public List<Double> getCriticalThresholdsForwards() {
            return criticalThresholdsForwards;
        }
        public List<Double> getCriticalThresholdsReverse() {
            return criticalThresholdsReverse;
        }
    }
}
