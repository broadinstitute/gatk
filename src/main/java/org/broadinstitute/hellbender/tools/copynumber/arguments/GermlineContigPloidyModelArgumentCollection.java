package org.broadinstitute.hellbender.tools.copynumber.arguments;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.tools.copynumber.DetermineGermlineContigPloidy;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class GermlineContigPloidyModelArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final String PLOIDY_CONCENTRATION_SCALE_LONG_NAME = "ploidy-concentration-scale";
    public static final String DEPTH_UPPER_BOUND_LONG_NAME = "depth-upper-bound";
    public static final String ERROR_RATE_UPPER_BOUND_LONG_NAME = "error-rate-upper-bound";
    public static final String CONTIG_BIAS_LOWER_BOUND_LONG_NAME = "contig-bias-lower-bound";
    public static final String CONTIG_BIAS_UPPER_BOUND_LONG_NAME = "contig-bias-upper-bound";
    public static final String CONTIG_BIAS_SCALE_LONG_NAME = "contig-bias-scale";

    @Argument(
            doc = "Scale factor for the concentration parameters of the per-contig-set Dirichlet prior on ploidy states.  " +
                    "The relative probabilities given by the ploidy-state priors are normalized and multiplied by this factor " +
                    "to yield the concentration parameters.",
            fullName = PLOIDY_CONCENTRATION_SCALE_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    private double ploidyConcentrationScale = 0.01;

    @Argument(
            doc = "Upper bound of the uniform prior on the per-sample depth.",
            fullName = DEPTH_UPPER_BOUND_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    private double depthUpperBound = 1000.;

    @Argument(
            doc = "Upper bound of the uniform prior on the error rate.",
            fullName = ERROR_RATE_UPPER_BOUND_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    private double errorRateUpperBound = 0.1;

    @Argument(
            doc = "Lower bound of the Gamma prior on the per-contig bias.",
            fullName = CONTIG_BIAS_LOWER_BOUND_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    private double contigBiasLowerBound = 0.1;

    @Argument(
            doc = "Upper bound of the Gamma prior on the per-contig bias.",
            fullName = CONTIG_BIAS_UPPER_BOUND_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    private double contigBiasUpperBound = 2.;

    @Argument(
            doc = "Scale factor for the Gamma prior on the per-contig bias.  " +
                    "Both alpha and beta hyperparameters for the Gamma prior will be set to this factor.",
            fullName = CONTIG_BIAS_SCALE_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    private double contigBiasScale = 10.;

    public List<String> generatePythonArguments(final DetermineGermlineContigPloidy.RunMode runMode) {
        if (runMode == DetermineGermlineContigPloidy.RunMode.COHORT) {
            return Arrays.asList(
                    String.format("--ploidy_concentration_scale=%e", ploidyConcentrationScale),
                    String.format("--depth_upper_bound=%e", depthUpperBound),
                    String.format("--error_rate_upper_bound=%e", errorRateUpperBound),
                    String.format("--contig_bias_lower_bound=%e", contigBiasLowerBound),
                    String.format("--contig_bias_upper_bound=%e", contigBiasUpperBound),
                    String.format("--contig_bias_scale=%e", contigBiasScale));
        }
        return Collections.emptyList();
    }

    public void validate() {
        ParamUtils.isPositive(ploidyConcentrationScale,
                "Ploidy concentration scale must be positive.");
        ParamUtils.isPositive(depthUpperBound,
                "Upper bound of the per-sample depth must be positive.");
        ParamUtils.isPositive(errorRateUpperBound,
                "Upper bound of the error rate must be positive.");
        ParamUtils.isPositive(contigBiasLowerBound,
                "Lower bound of the per-contig bias must be positive.");
        ParamUtils.isPositive(contigBiasUpperBound,
                "Upper bound of the per-contig bias must be positive.");
        ParamUtils.isPositive(contigBiasScale,
                "Scale of the per-contig bias must be positive.");
        Utils.validateArg(contigBiasLowerBound < contigBiasUpperBound,
                "Lower bound of the per-contig bias must be less than the upper bound.");
    }
}
