package org.broadinstitute.hellbender.utils.pairhmm;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;

public class DragstrCasesSamplerArgumentCollection {

    @Argument(doc = "Possible Gap-Penalty values for the DRAGstr model parameter esimation. " +
            "These are expressed in Phred scaled values with the following format: start:step:end. For example the default '10:1.0:50' indicate the sequence starting at 10 finishing at 50 sampled at 1.0 intervals.",
            optional = true,
            fullName = DragstrConstants.GP_VALUES_ARGUMENT_FULL_NAME)
    public DoubleSequence phredGpValues = DragstrConstants.DEFAULT_PHRED_GP_VALUES;

    @Argument(doc = "Possible a-priori probabilities for the heterozygous indel call for the DRAGstr model parameter esimation. " +
            "These are expressed in Phred scaled values with the following format: start:step:end. For example the default '10:1.0:50' indicate the sequence starting at 10 finishing at 50 sampled at 1.0 intervals.",
            optional = true,
            fullName = DragstrConstants.API_VALUES_ARGUMENT_FULL_NAME)
    public DoubleSequence phredApiValues = DragstrConstants.DEFAULT_PHRED_API_VALUES;

    @Argument(doc = "Possible a-priori probabilities for the heterozygous indel call for the DRAGstr model parameter esimation. " +
            "These are expressed in Phred scaled values with the following format: start:step:end. For example the default '10:1.0:50' indicate the sequence starting at 10 finishing at 50 sampled at 1.0 intervals.",
            optional = true,
            fullName = DragstrConstants.HET_TO_HOM_RATIO_FULL_NAME,
            minValue = 0.0,
            maxValue = Double.MAX_VALUE)
    public double hetToHomRatio = DragstrConstants.DEFAULT_HET_TO_HOM_RATIO;

    @Argument(doc = "Minimum number of sites for a repeat count and period length pair. We will combine pairs that have a smaller number of such cases (same period but +/- 1 repeat count)",
            optional = true,
            fullName = DragstrConstants.MIN_LOCI_COUNT_FULL_NAME,
            minValue = 1.0)
    public int minLociCount = DragstrConstants.DEFAULT_MIN_LOCI_COUNT;

    @Argument(doc = "<Not quite understand this one yet>",
            optional = true,
            fullName = DragstrConstants.API_MONO_THRESHOLD_FULL_NAME,
            minValue = 0.0)
    public int apiMonothresh = DragstrConstants.DEFAULT_API_MONO_THRESHOLD;

    @Argument(doc = "Minimum GOP value that is going to be ever estimated for any period and repeat count",
            optional = true,
            fullName = DragstrConstants.MIN_GOP_FULL_NAME,
            minValue = 0.0)
    public double minGOP = DragstrConstants.DEFAULT_MIN_GOP;

    @Argument(doc = "Maximum GOP value that is going to be ever estimated for any period and repeat count",
            optional = true,
            fullName = DragstrConstants.MAX_GOP_FULL_NAME,
            minValue = 0.0)
    public double maxGOP = DragstrConstants.DEFAULT_MAX_GOP;

    private void validate() {
        if (phredGpValues != DragstrConstants.DEFAULT_PHRED_GP_VALUES) {
            for (final double d : phredGpValues.toDoubleArray()) {
                if (!Double.isFinite(d) || d < 0) {
                    throw new CommandLineException.BadArgumentValue(DragstrConstants.GP_VALUES_ARGUMENT_FULL_NAME, "Not a valid Phred value: " + d);
                }
            }
        } else if (phredApiValues != DragstrConstants.DEFAULT_PHRED_API_VALUES) {
            for (final double d : phredApiValues.toDoubleArray()) {
                if (!Double.isFinite(d) || d < 0) {
                    throw new CommandLineException.BadArgumentValue(DragstrConstants.API_VALUES_ARGUMENT_FULL_NAME, "Not a valid Phred value: " + d);
                }
            }
        } else if (!Double.isFinite(hetToHomRatio) || hetToHomRatio <= 0) {
            throw new CommandLineException.BadArgumentValue(DragstrConstants.HET_TO_HOM_RATIO_FULL_NAME, "must be finite and greater than 0 but found " + hetToHomRatio);
        } else if (minLociCount < 1) {
            throw new CommandLineException.BadArgumentValue(DragstrConstants.MIN_LOCI_COUNT_FULL_NAME, "must be greater than 0 but found " + minLociCount);
        } else if (minGOP > maxGOP) {
            throw new CommandLineException.BadArgumentValue(DragstrConstants.MIN_GOP_FULL_NAME, "the min GOP " + minGOP + " cannot be larger than the max " + maxGOP);
        }
    }

    @Argument(fullName = DragstrConstants.MAX_PERIOD_ARGUMENT_FULL_NAME,
              doc = "maximum period considered in STR analysis",
              optional = true)
    public int maxPeriod = 8;

    @Argument(fullName = DragstrConstants.MAX_REPEATS_ARGUMENT_FULL_NAME,
              doc = "maximum repeat count considered in STR analysis",
              optional = true)
    public int maxRepeats = 20;

    @Argument(fullName = DragstrConstants.MIN_DEPTH_ARGUMENT_FULL_NAME,
              shortName = DragstrConstants.MIN_DEPTH_ARGUMENT_SHORT_NAME,
              doc = "Minimum coverage to consider a locus for sampling",
              optional = true)
    public int minDepth = DragstrConstants.DEFAULT_MIN_DEPTH;

    @Argument(fullName = DragstrConstants.USE_ALL_EVIDENCE_FULL_NAME,
            doc = "do not sample, simply use all cases available (slow for not downsampled STR loci sets)",
            optional = true)
    public boolean useAllEvidence = false;

    @Argument(fullName = DragstrConstants.RANDOM_SEED_ARGUMENT_FULL_NAME,
            doc = "random number generator base seed",
            optional = true)
    public int randomSeed = DragstrConstants.DEFAULT_RANDOM_SEED;

    @Argument(fullName = DragstrConstants.SAMPLING_PADDING_ARGUMENT_FULL_NAME,
            doc = "bases on either side of the repeat that are included in the STR pileup",
            optional = true)
    public int pileupPadding = DragstrConstants.DEFAULT_SAMPLING_PADDING;

    @Argument(fullName = DragstrConstants.SAMPLING_MIN_MQ_ARGUMENT_FULL_NAME,
            doc = "the minimum read mapping quality allowed in sampled loci. Any read with a lower MQ will result in discarding that locus",
            optional = true)
    public int samplingMinMQ = DragstrConstants.DEFAULT_SAMPLING_MIN_MQ;

    @Argument(fullName = DragstrConstants.SAMPLING_MAX_MQ_ARGUMENT_FULL_NAME,
            doc = "the maximum number of sites to sample per period and repeat count combination",
            optional = true)
    public int maxCount = DragstrConstants.DEFAULT_SAMPLING_MAX_COUNT;

    @Argument(fullName = DragstrConstants.SAMPLING_MIN_BQ_THRESHOLD_ARGUMENT_FULL_NAME,
            doc = "Base quality thershold for base-call. Reads that contain a number bases that map on the STR will not qualify for sampling",
            optional = true)
    public int baseQualThreshold = DragstrConstants.DEFAULT_SAMPLING_BASE_QUAL_THRESHOLD;

    @Argument(fullName = org.broadinstitute.hellbender.utils.pairhmm.DragstrConstants.SAMPLING_MAX_BQ_EXCEPTIONS_ALLOWED_ARGUMENT_FULL_NAME,
        doc ="Maximum number of STR overlapping base call with low quality allowed for any read to be considered for further analysis",
        optional =true)
    public int baseQualExceptionsAllowed = DragstrConstants.DEFAULT_SAMPLING_MAX_BQ_EXCEPTIONS_ALLOWED;

    @Argument(fullName = org.broadinstitute.hellbender.utils.pairhmm.DragstrConstants.SAMPLING_MAX_CASE_COUNT_ARGUMENT_FULL_NAME,
        doc ="maximum number of sites sampled for each combination of period and repeat count",
        optional =true)
    public int maximumNumberOfCases = DragstrConstants.DEFAULT_SAMPLING_MAX_COUNT;

    @Argument(fullName = org.broadinstitute.hellbender.utils.pairhmm.DragstrConstants.SAMPLING_MIN_NON_REF_CONTAING_CASE_COUNT_ARGUMENT_FULL_NAME,
        doc ="targeted minimum number of sites sampled that contain non-ref reads for each combination of period and repeat count",
        optional =true)
    public int targetMinimumNonRefCases =DragstrConstants.DEFAULT_SAMPLING_MIN_NON_REF_CONTAINING_CASE_COUNT;

    @Argument(fullName = org.broadinstitute.hellbender.utils.pairhmm.DragstrConstants.SAMPLING_MIN_CASE_COUNT_ARGUMENT_FULL_NAME,
        doc ="minimum number of sites sampled for each combination of period and repeat count",
        optional =true)
    public int minimumNumberOfCases =  DragstrConstants.DEFAULT_SAMPLING_MIN_CASE_COUNT;
}