package org.broadinstitute.hellbender.utils.pairhmm;


import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;

public class DragstrModelEstimatorArgumentCollection {

    public static final String GP_VALUES_ARGUMENT_FULL_NAME = "gp-values";
    public static final String API_VALUES_ARGUMENT_FULL_NAME = "api-values";
    public static final String HET_TO_HOM_RATIO_FULL_NAME = "het-to-hom-ratio";
    public static final String MIN_LOCI_COUNT_FULL_NAME = "min-loci-count";
    public static final String API_MONO_THRESHOD_FULL_NAME = "api-mono-threshold";
    public static final String MIN_GOP_FULL_NAME = "min-gop";
    public static final String MAX_GOP_FULL_NAME = "max-gop";
    public static final DoubleSequence DEFAULT_PHRED_GP_VALUES = new DoubleSequence("10:1.0:50");
    public static final DoubleSequence DEFAULT_PHRED_API_VALUES = new DoubleSequence("0:1.0:40");
    public static final double DEFAULT_HET_TO_HOM_RATIO = 2.0;
    public static final int DEFAULT_MIN_LOCI_COUNT = 50;
    public static final int DEFAULT_API_MONO_THRESHOLD = 3;
    public static final double DEFAULT_MIN_GOP = 10;
    public static final double DEFAULT_MAX_GOP = 50;

    @Argument(doc = "Possible Gap-Penalty values for the DRAGstr model parameter esimation. " +
            "These are expressed in Phred scaled values with the following format: start:step:end. For example the default '10:1.0:50' indicate the sequence starting at 10 finishing at 50 sampled at 1.0 intervals.",
             optional = true,
             fullName = GP_VALUES_ARGUMENT_FULL_NAME)
    public DoubleSequence phredGpValues = DEFAULT_PHRED_GP_VALUES;

    @Argument(doc = "Possible a-priori probabilities for the heterozygous indel call for the DRAGstr model parameter esimation. " +
            "These are expressed in Phred scaled values with the following format: start:step:end. For example the default '10:1.0:50' indicate the sequence starting at 10 finishing at 50 sampled at 1.0 intervals.",
            optional = true,
            fullName = API_VALUES_ARGUMENT_FULL_NAME)
    public DoubleSequence phredApiValues = DEFAULT_PHRED_API_VALUES;

    @Argument(doc = "Possible a-priori probabilities for the heterozygous indel call for the DRAGstr model parameter esimation. " +
            "These are expressed in Phred scaled values with the following format: start:step:end. For example the default '10:1.0:50' indicate the sequence starting at 10 finishing at 50 sampled at 1.0 intervals.",
            optional = true,
            fullName = HET_TO_HOM_RATIO_FULL_NAME,
            minValue = 0.0,
            maxValue = Double.MAX_VALUE)
    public double hetToHomRatio = DEFAULT_HET_TO_HOM_RATIO;

    @Argument(doc = "Minimum number of sites for a repeat count and period length pair. We will combine pairs that have a smaller number of such cases (same period but +/- 1 repeat count)",
              optional = true,
              fullName = MIN_LOCI_COUNT_FULL_NAME,
              minValue = 1.0)
    public int minLociCount = DEFAULT_MIN_LOCI_COUNT;

    @Argument(doc = "<Not quite understand this one yet>",
              optional = true,
              fullName = API_MONO_THRESHOD_FULL_NAME,
              minValue = 0.0)
    public int apiMonothresh = DEFAULT_API_MONO_THRESHOLD;

    @Argument(doc = "Minimum GOP value that is going to be ever estimated for any period and repeat count",
              optional = true,
              fullName = MIN_GOP_FULL_NAME,
              minValue = 0.0)
    public double minGOP = DEFAULT_MIN_GOP;

    @Argument(doc = "Maximum GOP value that is going to be ever estimated for any period and repeat count",
            optional = true,
            fullName = MAX_GOP_FULL_NAME,
            minValue = 0.0)
    public double maxGOP = DEFAULT_MAX_GOP;

    @Argument(fullName = "dont-adjust-gop", optional = true)
    public boolean dontPostAdjustmentOfGOP = false;

    private void validate() {
        if (phredGpValues != DEFAULT_PHRED_GP_VALUES) {
            for (final double d : phredGpValues.toDoubleArray()) {
                if (!Double.isFinite(d) || d < 0) {
                    throw new CommandLineException.BadArgumentValue(GP_VALUES_ARGUMENT_FULL_NAME, "Not a valid Phred value: " + d);
                }
            }
        } else if (phredApiValues != DEFAULT_PHRED_API_VALUES) {
            for (final double d : phredApiValues.toDoubleArray()) {
                if (!Double.isFinite(d) || d < 0) {
                    throw new CommandLineException.BadArgumentValue(API_VALUES_ARGUMENT_FULL_NAME, "Not a valid Phred value: " + d);
                }
            }
        } else if (!Double.isFinite(hetToHomRatio) || hetToHomRatio <= 0) {
            throw new CommandLineException.BadArgumentValue(HET_TO_HOM_RATIO_FULL_NAME, "must be finite and greater than 0 but found " + hetToHomRatio);
        } else if (minLociCount < 1) {
            throw new CommandLineException.BadArgumentValue(MIN_LOCI_COUNT_FULL_NAME, "must be greater than 0 but found " + minLociCount);
        } else if (minGOP > maxGOP) {
            throw new CommandLineException.BadArgumentValue(MIN_GOP_FULL_NAME, "the min GOP " + minGOP + " cannot be larger than the max " + maxGOP);
        }
    }
}
