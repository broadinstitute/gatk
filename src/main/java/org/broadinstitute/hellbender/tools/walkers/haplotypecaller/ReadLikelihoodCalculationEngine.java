package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.function.ToDoubleFunction;

/**
 * Common interface for assembly-haplotype vs reads likelihood engines.
 */
public interface ReadLikelihoodCalculationEngine extends AutoCloseable {
    final int MAX_STR_UNIT_LENGTH = 8;
    final int MAX_REPEAT_LENGTH   = 20;

    enum Implementation {
        /**
         * Classic full pair-hmm all haplotypes vs all reads.
         */
        PairHMM,

        FlowBased,

        FlowBasedHMM
    }


    /**
     * Calculates the likelihood of reads across many samples evaluated against haplotypes resulting from the
     * active region assembly process.
     *
     * @param assemblyResultSet the input assembly results.
     * @param samples the list of targeted samples.
     * @param perSampleReadList the input read sets stratified per sample.
     *
     * @throws IllegalArgumentException if any parameter is {@code null}.
     *
     * @return never {@code null}, and with at least one entry for input sample (keys in {@code perSampleReadList}.
     *    The value maps can be potentially empty though.
     */
    public default AlleleLikelihoods<GATKRead, Haplotype> computeReadLikelihoods(AssemblyResultSet assemblyResultSet, SampleList samples,
                                                                         Map<String, List<GATKRead>> perSampleReadList,
                                                                         boolean filterPoorly) {
        Utils.nonNull(assemblyResultSet, "assemblyResultSet is null");
        final List<Haplotype> haplotypeList = assemblyResultSet.getHaplotypeList();
        final SAMFileHeader hdr= assemblyResultSet.getRegionForGenotyping().getHeader();

        return computeReadLikelihoods(haplotypeList, hdr, samples, perSampleReadList, filterPoorly);
    }

    /**
     * Implementation of computeReadLikelihoods that is intended to be implemented by implementers. Calculates the likelihood of the reads
     * across many samples evaluated against haplotypes resulting from the active region assembly prcess.
     *
     * @param haplotypeList List of haplotypes over which to to score likelihoods.
     * @param hdr header to extract read metadata for calling from.
     * @param samples the list of targeted samples.
     * @param perSampleReadList input read sets stratified per sample.
     *
     * @throws IllegalArgumentException if any parameter is {@code null}.
     *
     * @return never {@code null}, and with at least one entry for input sample (keys in {@code perSampleReadList}.
     *    The value maps can be potentially empty though.
     */
    public AlleleLikelihoods<GATKRead, Haplotype> computeReadLikelihoods(final List<Haplotype> haplotypeList,
                                                                         final SAMFileHeader hdr,
                                                                         final SampleList samples,
                                                                         final Map<String, List<GATKRead>> perSampleReadList, final boolean filterPoorly);

    /**
     * This method must be called when the client is done with likelihood calculations.
     * It closes any open resources.
     */
    @Override
    public void close();


    /**
     * Static methods related to read disqualification to be used by implementing classes of this interface
     *
     * @param result                        AlleleLikelihoods object from which reads will be removed
     * @param dynamicDisqualification       whether to apply DRAGEN-GATK dynamic disqualification threshold to better handle low BQ reads
     * @param expectedErrorRatePerBase      amount of uncertainty to tolerate per read-base (in log scale)
     * @param readDisqualificationScale     constant used to scale the dynamic read disqualificaiton threshold
     */
    default void filterPoorlyModeledEvidence(final AlleleLikelihoods<GATKRead, Haplotype> result, final boolean dynamicDisqualification, final double expectedErrorRatePerBase, final double readDisqualificationScale) {
        if (dynamicDisqualification) {
            result.filterPoorlyModeledEvidence(daynamicLog10MinLiklihoodModel(readDisqualificationScale, log10MinTrueLikelihood(expectedErrorRatePerBase, false)));
        } else {
            result.filterPoorlyModeledEvidence(log10MinTrueLikelihood(expectedErrorRatePerBase, true));
        }
    }

    static ToDoubleFunction<GATKRead> daynamicLog10MinLiklihoodModel(final double dynamicRadQualConstant, final ToDoubleFunction<GATKRead> log10MinTrueLikelihood) {
        return read -> {
            final double dynamicThreshold = calculateLog10DynamicReadQualThreshold(read, dynamicRadQualConstant);
            final double log10MaxLikelihoodForTrueAllele = log10MinTrueLikelihood.applyAsDouble(read);
            if (dynamicThreshold < log10MaxLikelihoodForTrueAllele ) {
                if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
                    HaplotypeCallerGenotypingDebugger.println("For read "+ read.getName() + " replacing old threshold ("+log10MaxLikelihoodForTrueAllele+") with new threshold: "+dynamicThreshold);
                }
                return dynamicThreshold;
            } else {
                return log10MaxLikelihoodForTrueAllele;
            }
        };
    }

    /**
     * Standard read disqulification threshold algorithm
     */
    default ToDoubleFunction<GATKRead> log10MinTrueLikelihood(final double maximumErrorPerBase, final boolean capLikelihoods) {
        return read -> {
            // TODO this might be replaced by an explicit calculation
            final int qualifiedReadLength = read.getTransientAttribute(PairHMMLikelihoodCalculationEngine.HMM_BASE_QUALITIES_TAG) != null ? ((byte[])read.getTransientAttribute(PairHMMLikelihoodCalculationEngine.HMM_BASE_QUALITIES_TAG)).length : read.getLength();
            final double maxErrorsForRead = capLikelihoods ? Math.min(2.0, Math.ceil(qualifiedReadLength * maximumErrorPerBase)) : Math.ceil(qualifiedReadLength * maximumErrorPerBase);
            final double log10QualPerBase = -4.0;
            return maxErrorsForRead * log10QualPerBase;
        };
    }

    /**
     * Dynamic likelihood DRAGEN-GATK threshold.
     */
    static double calculateLog10DynamicReadQualThreshold(final GATKRead read, final double dynamicReadQualConstant) {
        double sumMean = 0;
        double sumVariance = 0;

        final byte[] baseQualities = read.getOptionalTransientAttribute(PairHMMLikelihoodCalculationEngine.HMM_BASE_QUALITIES_TAG, byte[].class)
                .orElseGet(read::getBaseQualities);

        for (final int qualByte : baseQualities) {
            final int bq = 0xFF & qualByte; // making sure that larger BQ are not casted into negatives.
            // bound the base qualities for lookup between 1 and 40
            final int entryIndex = bq <= 1 ? 0 : Math.min(MAXIMUM_DYNAMIC_QUAL_THRESHOLD_ENTRY_BASEQ, bq) - 1;
            final int meanOffset = entryIndex * DYNAMIC_QUAL_THRESHOLD_TABLE_ENTRY_LENGTH + DYNAMIC_QUAL_THRESHOLD_TABLE_ENTRY_MEAN_OFFSET;
            final int varOffset = meanOffset + 1;
            sumMean +=      dynamicReadQualThreshLookupTable[meanOffset];
            sumVariance +=  dynamicReadQualThreshLookupTable[varOffset];
        }

        final double threshold = sumMean + dynamicReadQualConstant * Math.sqrt(sumVariance);
        return QualityUtils.qualToErrorProbLog10(threshold); // = threshold * -.1;
    }

    // TODO i don't like having a lookup table be static like this, i would prefer this be computed at initialization (with the default values being saved as a test)
    // table used for disqualifying reads for genotyping
    // Format for each row of table: baseQ, mean, variance
    // Actual threshold is calculated over the length of the read as:
    // sum(means) + K * sqrt(sum(variances))
    static double dynamicReadQualThreshLookupTable[] = {
            //baseQ,mean,variance
            1,  5.996842844, 0.196616587, 2,  5.870018422, 1.388545569, 3,  5.401558531, 5.641990128,
            4,  4.818940919, 10.33176216, 5,  4.218758304, 14.25799688, 6,  3.646319832, 17.02880749,
            7,  3.122346753, 18.64537883, 8,  2.654731979, 19.27521677, 9,  2.244479156, 19.13584613,
            10, 1.88893867,  18.43922003, 11, 1.583645342, 17.36842261, 12, 1.3233807, 16.07088712,
            13, 1.102785365, 14.65952563, 14, 0.916703025, 13.21718577, 15, 0.760361881, 11.80207947,
            16, 0.629457387, 10.45304833, 17, 0.520175654, 9.194183767, 18, 0.42918208,  8.038657241,
            19, 0.353590663, 6.991779595, 20, 0.290923699, 6.053379213, 21, 0.23906788,  5.219610436,
            22, 0.196230431, 4.484302033, 23, 0.160897421, 3.839943445, 24, 0.131795374, 3.27839108,
            25, 0.1078567,   2.791361596, 26, 0.088189063, 2.370765375, 27, 0.072048567, 2.008921719,
            28, 0.058816518, 1.698687797, 29, 0.047979438, 1.433525748, 30, 0.039111985, 1.207526336,
            31, 0.031862437, 1.015402928, 32, 0.025940415, 0.852465956, 33, 0.021106532, 0.714585285,
            34, 0.017163711, 0.598145851, 35, 0.013949904, 0.500000349, 36, 0.011332027, 0.41742159,
            37, 0.009200898, 0.348056286, 38, 0.007467036, 0.289881373, 39, 0.006057179, 0.241163527,
            40, 0.004911394, 0.200422214};

    static final int MAXIMUM_DYNAMIC_QUAL_THRESHOLD_ENTRY_BASEQ = 40;
    static final int DYNAMIC_QUAL_THRESHOLD_TABLE_ENTRY_LENGTH = 3;
    static final int DYNAMIC_QUAL_THRESHOLD_TABLE_ENTRY_MEAN_OFFSET = 1;


    /**
     * The expected rate of random sequencing errors for a read originating from its true haplotype.
     *
     * For example, if this is 0.01, then we'd expect 1 error per 100 bp.
     */
    double DEFAULT_EXPECTED_ERROR_RATE_PER_BASE = 0.02;

    @VisibleForTesting
    static Pair<byte[], Integer> findTandemRepeatUnits(byte[] readBases, int offset) {
        int maxBW = 0;
        byte[] bestBWRepeatUnit = {readBases[offset]};
        for (int str = 1; str <= MAX_STR_UNIT_LENGTH; str++) {
            // fix repeat unit length
            //edge case: if candidate tandem repeat unit falls beyond edge of read, skip
            if (offset+1-str < 0) {
                break;
            }

            // get backward repeat unit and # repeats
            maxBW = GATKVariantContextUtils.findNumberOfRepetitions(readBases, offset - str + 1,  str , readBases, 0, offset + 1, false);
            if (maxBW > 1) {
                bestBWRepeatUnit = Arrays.copyOfRange(readBases, offset - str + 1, offset + 1);
                break;
            }
        }
        byte[] bestRepeatUnit = bestBWRepeatUnit;
        int maxRL = maxBW;

        if (offset < readBases.length-1) {
            byte[] bestFWRepeatUnit = {readBases[offset+1]};
            int maxFW = 0;

            for (int str = 1; str <= MAX_STR_UNIT_LENGTH; str++) {
                // fix repeat unit length
                //edge case: if candidate tandem repeat unit falls beyond edge of read, skip
                if (offset+str+1 > readBases.length) {
                    break;
                }

                // get forward repeat unit and # repeats
                maxFW = GATKVariantContextUtils.findNumberOfRepetitions(readBases, offset + 1, str, readBases, offset + 1, readBases.length-offset -1, true);
                if (maxFW > 1) {
                    bestFWRepeatUnit = Arrays.copyOfRange(readBases, offset + 1, offset+str+1);
                    break;
                }
            }
            // if FW repeat unit = BW repeat unit it means we're in the middle of a tandem repeat - add FW and BW components
            if (Arrays.equals(bestFWRepeatUnit, bestBWRepeatUnit)) {
                maxRL = maxBW + maxFW;
                bestRepeatUnit = bestFWRepeatUnit; // arbitrary
            } else {
                // tandem repeat starting forward from current offset.
                // It could be the case that best BW unit was different from FW unit, but that BW still contains FW unit.
                // For example, TTCTT(C) CCC - at (C) place, best BW unit is (TTC)2, best FW unit is (C)3.
                // but correct representation at that place might be (C)4.
                // Hence, if the FW and BW units don't match, check if BW unit can still be a part of FW unit and add
                // representations to total
                final byte[] testString = Arrays.copyOfRange(readBases, 0, offset + 1);
                maxBW = GATKVariantContextUtils.findNumberOfRepetitions(bestFWRepeatUnit, testString, false);
                maxRL = maxFW + maxBW;
                bestRepeatUnit = bestFWRepeatUnit;
            }
        }

        if(maxRL > MAX_REPEAT_LENGTH) {
            maxRL = MAX_REPEAT_LENGTH;
        }
        return Pair.of(bestRepeatUnit, maxRL);
    }


}
