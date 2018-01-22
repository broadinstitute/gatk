package org.broadinstitute.hellbender.tools.copynumber.arguments;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.tools.copynumber.GermlineCNVCaller;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class GermlineDenoisingModelArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final String MAX_BIAS_FACTORS_LONG_NAME = "max-bias-factors";
    public static final String MAPPING_ERROR_RATE_LONG_NAME = "mapping-error-rate";
    public static final String INTERVAL_PSI_SCALE_LONG_NAME = "interval-psi-scale";
    public static final String SAMPLE_PSI_SCALE_LONG_NAME = "sample-psi-scale";
    public static final String DEPTH_CORRECTION_TAU_LONG_NAME = "depth-correction-tau";
    public static final String LOG_MEAN_BIAS_STANDARD_DEVIATION_LONG_NAME = "log-mean-bias-standard-deviation";
    public static final String INIT_ARD_REL_UNEXPLAINED_VARIANCE_LONG_NAME = "init-ard-rel-unexplained-variance";
    public static final String NUM_GC_BINS_LONG_NAME = "num-gc-bins";
    public static final String GC_CURVE_STANDARD_DEVIATION_LONG_NAME = "gc-curve-standard-deviation";
    public static final String COPY_NUMBER_POSTERIOR_EXPECTATION_MODE_LONG_NAME = "copy-number-posterior-expectation-mode";
    public static final String ENABLE_BIAS_FACTORS_LONG_NAME = "enable-bias-factors";
    public static final String ACTIVE_CLASS_PADDING_HYBRID_MODE_LONG_NAME = "active-class-padding-hybrid-mode";

    public enum CopyNumberPosteriorExpectationMode {
        MAP("map"),
        EXACT("exact"),
        HYBRID("hybrid");

        final String pythonArgumentString;

        CopyNumberPosteriorExpectationMode(final String pythonArgumentString) {
            this.pythonArgumentString = pythonArgumentString;
        }
    }

    @Argument(
            doc = "Maximum number of bias factors.",
            fullName = MAX_BIAS_FACTORS_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private int maxBiasFactors = 5;

    @Argument(
            doc = "Typical mapping error rate.",
            fullName = MAPPING_ERROR_RATE_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    private double mappingErrorRate = 0.01;

    @Argument(
            doc = "Typical scale of interval-specific unexplained variance.",
            fullName = INTERVAL_PSI_SCALE_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    private double intervalPsiScale = 0.001;

    @Argument(
            doc = "Typical scale of sample-specific correction to the unexplained variance.",
            fullName = SAMPLE_PSI_SCALE_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    private double samplePsiScale = 0.0001;

    @Argument(
            doc = "Precision of read depth pinning to its global value.",
            fullName = DEPTH_CORRECTION_TAU_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    private double depthCorrectionTau = 10000.0;

    @Argument(
            doc = "Standard deviation of log mean bias.",
            fullName = LOG_MEAN_BIAS_STANDARD_DEVIATION_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    private double logMeanBiasStandardDeviation = 0.1;

    @Argument(
            doc = "Initial value of ARD prior precisions relative to the scale of interval-specific " +
                    "unexplained variance.",
            fullName = INIT_ARD_REL_UNEXPLAINED_VARIANCE_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    private double initARDRelUnexplainedVariance = 0.1;

    @Argument(
            doc = "Number of bins used to represent the GC-bias curves.",
            fullName = NUM_GC_BINS_LONG_NAME,
            minValue = 1,
            optional = true
    )
    private int numGCBins = 20;

    @Argument(
            doc = "Prior standard deviation of the GC curve from flat.",
            fullName = GC_CURVE_STANDARD_DEVIATION_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    private double gcCurveStandardDeviation = 1.;

    @Argument(
            doc = "The strategy for calculating copy number posterior expectations in the coverage denoising model.",
            fullName = COPY_NUMBER_POSTERIOR_EXPECTATION_MODE_LONG_NAME,
            optional = true
    )
    private CopyNumberPosteriorExpectationMode copyNumberPosteriorExpectationMode =
            CopyNumberPosteriorExpectationMode.HYBRID;

    @Argument(
            doc = "Enable discovery of bias factors.",
            fullName = ENABLE_BIAS_FACTORS_LONG_NAME,
            optional = true
    )
    private boolean enableBiasFactors = true;

    @Argument(
            doc = "If copy-number-posterior-expectation-mode is set to HYBRID, CNV-active intervals determined " +
                    "at any time will be padded by this value (in the units of bp) in order to obtain the set of " +
                    "intervals on which copy number posterior expectation is performed exactly.",
            fullName = ACTIVE_CLASS_PADDING_HYBRID_MODE_LONG_NAME,
            optional = true
    )
    private int activeClassPaddingHybridMode = 50000;

    /**
     * Generates arguments for the python CLI tool. Note that 'enable_explicit_gc_bias_modeling' is added
     * by {@link GermlineCNVCaller}.
     */
    public List<String> generatePythonArguments(final GermlineCNVCaller.RunMode runMode) {
        final List<String> arguments = new ArrayList<>(Arrays.asList(
                String.format("--psi_s_scale=%e", samplePsiScale),
                String.format("--mapping_error_rate=%e", mappingErrorRate),
                String.format("--depth_correction_tau=%e", depthCorrectionTau),
                String.format("--q_c_expectation_mode=%s", copyNumberPosteriorExpectationMode.pythonArgumentString)));
        if (runMode == GermlineCNVCaller.RunMode.COHORT) {
            arguments.addAll(Arrays.asList(
                    String.format("--max_bias_factors=%d", maxBiasFactors),
                    String.format("--psi_t_scale=%e", intervalPsiScale),
                    String.format("--log_mean_bias_std=%e", logMeanBiasStandardDeviation),
                    String.format("--init_ard_rel_unexplained_variance=%e", initARDRelUnexplainedVariance),
                    String.format("--num_gc_bins=%d", numGCBins),
                    String.format("--gc_curve_sd=%e", gcCurveStandardDeviation),
                    String.format("--active_class_padding_hybrid_mode=%d", activeClassPaddingHybridMode)));
            if (enableBiasFactors) {
                arguments.add("--enable_bias_factors=True");
            } else {
                arguments.add("--enable_bias_factors=False");
            }
            //this python argument is not exposed but we add it for completeness and logging purposes
            arguments.add("--disable_bias_factors_in_active_class=False");
        }
        return arguments;
    }

    public void validate() {
        ParamUtils.isPositive(maxBiasFactors,
                "Maximum number of bias factors must be positive.");
        ParamUtils.isPositive(mappingErrorRate,
                "Mapping error rate must be positive.");
        ParamUtils.isPositive(intervalPsiScale,
                "Typical scale of interval-specific unexplained variance must be positive.");
        ParamUtils.isPositive(samplePsiScale,
                "Typical scale of sample-specific correction to the unexplained variance must be positive.");
        ParamUtils.isPositive(depthCorrectionTau,
                "Precision of read depth pinning to its global value must be positive.");
        ParamUtils.isPositive(logMeanBiasStandardDeviation,
                "Standard deviation of log mean bias must be positive.");
        ParamUtils.isPositive(initARDRelUnexplainedVariance,
                "Initial value of ARD prior precision relative to the scale of " +
                        "interval-specific unexplained variance must be positive.");
        Utils.validateArg(numGCBins >= 2,
                "Number of bins used to represent the GC-bias curves must be at least 2.");
        ParamUtils.isPositive(gcCurveStandardDeviation,
                "Prior standard deviation of the GC curve from flat must be positive.");
        ParamUtils.isPositiveOrZero(activeClassPaddingHybridMode,
                "CNV-active interval padding in HYBRID copy-number posterior expectation mode must " +
                        "be non-negative.");
    }
}
