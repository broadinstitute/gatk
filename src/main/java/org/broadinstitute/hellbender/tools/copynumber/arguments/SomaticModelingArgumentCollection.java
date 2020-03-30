package org.broadinstitute.hellbender.tools.copynumber.arguments;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;

public class SomaticModelingArgumentCollection implements Serializable {
    public static final long serialVersionUID = 1L;

    //MCMC argument names
    public static final String MINOR_ALLELE_FRACTION_PRIOR_ALPHA_LONG_NAME = "minor-allele-fraction-prior-alpha";
    public static final String NUMBER_OF_SAMPLES_COPY_RATIO_LONG_NAME = "number-of-samples-copy-ratio";
    public static final String NUMBER_OF_BURN_IN_SAMPLES_COPY_RATIO_LONG_NAME = "number-of-burn-in-samples-copy-ratio";
    public static final String NUMBER_OF_SAMPLES_ALLELE_FRACTION_LONG_NAME = "number-of-samples-allele-fraction";
    public static final String NUMBER_OF_BURN_IN_SAMPLES_ALLELE_FRACTION_LONG_NAME = "number-of-burn-in-samples-allele-fraction";

    //smoothing argument names
    public static final String SMOOTHING_CREDIBLE_INTERVAL_THRESHOLD_COPY_RATIO_LONG_NAME = "smoothing-credible-interval-threshold-copy-ratio";
    public static final String SMOOTHING_CREDIBLE_INTERVAL_THRESHOLD_ALLELE_FRACTION_LONG_NAME = "smoothing-credible-interval-threshold-allele-fraction";
    public static final String MAXIMUM_NUMBER_OF_SMOOTHING_ITERATIONS_LONG_NAME = "maximum-number-of-smoothing-iterations";
    public static final String NUMBER_OF_SMOOTHING_ITERATIONS_PER_FIT_LONG_NAME = "number-of-smoothing-iterations-per-fit";

    @Argument(
            doc = "Alpha hyperparameter for the 4-parameter beta-distribution prior on segment minor-allele fraction. " +
                    "The prior for the minor-allele fraction f in each segment is assumed to be Beta(alpha, 1, 0, 1/2). " +
                    "Increasing this hyperparameter will reduce the effect of reference bias at the expense of sensitivity.",
            fullName = MINOR_ALLELE_FRACTION_PRIOR_ALPHA_LONG_NAME,
            optional = true,
            minValue = 1
    )
    public double minorAlleleFractionPriorAlpha = 25.;

    @Argument(
            doc = "Total number of MCMC samples for copy-ratio model.",
            fullName = NUMBER_OF_SAMPLES_COPY_RATIO_LONG_NAME,
            optional = true,
            minValue = 1
    )
    public int numSamplesCopyRatio = 100;

    @Argument(
            doc = "Number of burn-in samples to discard for copy-ratio model.",
            fullName = NUMBER_OF_BURN_IN_SAMPLES_COPY_RATIO_LONG_NAME,
            optional = true,
            minValue = 0
    )
    public int numBurnInCopyRatio = 50;

    @Argument(
            doc = "Total number of MCMC samples for allele-fraction model.",
            fullName = NUMBER_OF_SAMPLES_ALLELE_FRACTION_LONG_NAME,
            optional = true,
            minValue = 1
    )
    public int numSamplesAlleleFraction = 100;

    @Argument(
            doc = "Number of burn-in samples to discard for allele-fraction model.",
            fullName = NUMBER_OF_BURN_IN_SAMPLES_ALLELE_FRACTION_LONG_NAME,
            optional = true,
            minValue = 0
    )
    public int numBurnInAlleleFraction = 50;

    @Argument(
            doc = "Number of 10% equal-tailed credible-interval widths to use for copy-ratio segmentation smoothing.",
            fullName = SMOOTHING_CREDIBLE_INTERVAL_THRESHOLD_COPY_RATIO_LONG_NAME,
            optional = true,
            minValue = 0.
    )
    public double smoothingCredibleIntervalThresholdCopyRatio = 2.;

    @Argument(
            doc = "Number of 10% equal-tailed credible-interval widths to use for allele-fraction segmentation smoothing.",
            fullName = SMOOTHING_CREDIBLE_INTERVAL_THRESHOLD_ALLELE_FRACTION_LONG_NAME,
            optional = true,
            minValue = 0.
    )
    public double smoothingCredibleIntervalThresholdAlleleFraction = 2.;

    @Argument(
            doc = "Maximum number of iterations allowed for segmentation smoothing.",
            fullName = MAXIMUM_NUMBER_OF_SMOOTHING_ITERATIONS_LONG_NAME,
            optional = true,
            minValue = 0
    )
    public int maxNumSmoothingIterations = 25;

    @Argument(
            doc = "Number of segmentation-smoothing iterations per MCMC model refit. " +
                    "(Increasing this will decrease runtime, but the final number of segments may be higher. " +
                    "Setting this to 0 will completely disable model refitting between iterations.)",
            fullName = NUMBER_OF_SMOOTHING_ITERATIONS_PER_FIT_LONG_NAME,
            optional = true,
            minValue = 0
    )
    public int numSmoothingIterationsPerFit = 0;

    public void validateArguments() {
        Utils.validateArg(numSamplesCopyRatio > numBurnInCopyRatio,
                "Number of copy-ratio samples must be greater than number of copy-ratio burn-in samples.");
        Utils.validateArg(numSamplesAlleleFraction > numBurnInAlleleFraction,
                "Number of allele-fraction samples must be greater than number of allele-fraction burn-in samples.");
    }
}
