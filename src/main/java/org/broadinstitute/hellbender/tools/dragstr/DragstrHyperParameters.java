package org.broadinstitute.hellbender.tools.dragstr;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.utils.dragstr.DoubleSequence;

public final class DragstrHyperParameters {

    public static final int DEFAULT_MAX_PERIOD = 8;
    public static final int DEFAULT_MAX_REPEAT_LENGTH = 20;

    public static final DoubleSequence DEFAULT_PHRED_GP_VALUES = new DoubleSequence("10:1.0:50");

    public static final String GP_VALUES_ARGUMENT_FULL_NAME = "gp-values";
    public static final DoubleSequence DEFAULT_PHRED_API_VALUES = new DoubleSequence("0:1.0:40");
    public static final String API_VALUES_ARGUMENT_FULL_NAME = "api-values";
    public static final DoubleSequence DEFAULT_PHRED_GOP_VALUES = new DoubleSequence("10:.25:50");
    public static final String GOP_VALUES_ARGUMENT_FULL_NAME = "gop-values";
    public static final double DEFAULT_HET_TO_HOM_RATIO = 2.0;
    public static final String HET_TO_HOM_RATIO_FULL_NAME = "het-to-hom-ratio";
    public static final int DEFAULT_MIN_LOCI_COUNT = 50;
    public static final String MIN_LOCI_COUNT_FULL_NAME = "min-loci-count";
    public static final int DEFAULT_API_MONO_THRESHOLD = 3;
    public static final String API_MONO_THRESHOLD_FULL_NAME = "api-mono-threshold";
    public static final String MAX_PERIOD_ARGUMENT_FULL_NAME = "max-period";
    public static final String MAX_REPEATS_ARGUMENT_FULL_NAME = "max-repeats";
    public static final int DEFAULT_MIN_DEPTH = 10;
    public static final String MIN_DEPTH_ARGUMENT_SHORT_NAME = "md";
    public static final String MIN_DEPTH_ARGUMENT_FULL_NAME = "minimum-depth";
    public static final int DEFAULT_SAMPLING_PADDING = 5;
    public static final String SAMPLING_PADDING_ARGUMENT_FULL_NAME = "pileup-padding";
    public static final int DEFAULT_SAMPLING_MIN_MQ = 20;
    public static final String SAMPLING_MIN_MQ_ARGUMENT_FULL_NAME = "sampling-min-mq";



    @Argument(doc = "Possible Gap-Penalty values for the DRAGstr model parameter esimation. " +
            "These are expressed in Phred scaled values with the following format: start:step:end. For example the default '10:1.0:50' indicate the sequence starting at 10 finishing at 50 sampled at 1.0 intervals.",
            optional = true,
            fullName = GP_VALUES_ARGUMENT_FULL_NAME)
    DoubleSequence phredGpValues = DEFAULT_PHRED_GP_VALUES;

    @Argument(doc = "Possible prior probabilities for the heterozygous indel call for the DRAGstr model parameter esimation. " +
            "These are expressed in Phred scaled values with the following format: start:step:end. For example the default '10:1.0:50' indicate the sequence starting at 10 finishing at 50 sampled at 1.0 intervals.",
            optional = true,
            fullName = API_VALUES_ARGUMENT_FULL_NAME)
    DoubleSequence phredApiValues = DEFAULT_PHRED_API_VALUES;

    @Argument(doc = "Possible values for the gop parmeter" +
            "These are expressed in Phred scaled values with the following format: start:step:end. For example the default '10:1.0:50' indicate the sequence starting at 10 finishing at 50 sampled at 1.0 intervals.",
            optional = true,
            fullName = GOP_VALUES_ARGUMENT_FULL_NAME)
    DoubleSequence phredGopValues = DEFAULT_PHRED_GOP_VALUES;

    @Argument(doc = "Possible prior probabilities for the heterozygous indel call for the DRAGstr model parameter esimation. " +
            "These are expressed in Phred scaled values with the following format: start:step:end. For example the default '10:1.0:50' indicate the sequence starting at 10 finishing at 50 sampled at 1.0 intervals.",
            optional = true,
            fullName = HET_TO_HOM_RATIO_FULL_NAME,
            minValue = 0.0,
            maxValue = Double.MAX_VALUE)
    double hetToHomRatio = DEFAULT_HET_TO_HOM_RATIO;

    @Argument(doc = "Minimum number of sites for a repeat count and period length pair. The estimation will combine pairs that have a smaller number of such cases (same period but +/- 1 repeat count)",
            optional = true,
            fullName = MIN_LOCI_COUNT_FULL_NAME,
            minValue = 1.0)
    int minLociCount = DEFAULT_MIN_LOCI_COUNT;

    @Argument(doc = "Maximum drop allowed in the API parameter within a period between consecutive repeat length values in Phred scale",
            optional = true,
            fullName = API_MONO_THRESHOLD_FULL_NAME,
            minValue = 0.0)
    int apiMonothresh = DEFAULT_API_MONO_THRESHOLD;

    @Argument(fullName = MAX_PERIOD_ARGUMENT_FULL_NAME,
            doc = "maximum period considered in STR analysis",
            optional = true, minValue = 2, maxValue = 10)
    public int maxPeriod = DEFAULT_MAX_PERIOD;

    @Argument(fullName = MAX_REPEATS_ARGUMENT_FULL_NAME,
            doc = "maximum repeat length (in repeated units) considered in STR analyses",
            optional = true, minValue = 2, maxValue = 100)
    public int maxRepeatLength = DEFAULT_MAX_REPEAT_LENGTH;

    @Argument(fullName = MIN_DEPTH_ARGUMENT_FULL_NAME,
            shortName = MIN_DEPTH_ARGUMENT_SHORT_NAME,
            doc = "Minimum coverage to consider a locus for sampling",
            optional = true)
    public int minDepth = DEFAULT_MIN_DEPTH;

    @Argument(fullName = SAMPLING_PADDING_ARGUMENT_FULL_NAME,
            doc = "bases on either side of the repeat that are included in the STR pileup",
            optional = true)
    public int strPadding = DEFAULT_SAMPLING_PADDING;

    @Argument(fullName = SAMPLING_MIN_MQ_ARGUMENT_FULL_NAME,
            doc = "the minimum read mapping quality allowed in sampled loci. Any read with a lower MQ will result in discarding that locus",
            optional = true)
    public int minMQ = DEFAULT_SAMPLING_MIN_MQ;

    /**
     * Performs some validation constraint changes on the values provided to this argument collection.
     * <p>
     *     Will throw the appropriate user exception if it detects any inconsistency.
     * </p>
     */
    public void validate() {
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
        }
    }
}