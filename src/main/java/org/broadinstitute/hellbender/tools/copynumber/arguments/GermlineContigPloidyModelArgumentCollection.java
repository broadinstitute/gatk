package org.broadinstitute.hellbender.tools.copynumber.arguments;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.tools.copynumber.DetermineGermlineContigPloidy;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class GermlineContigPloidyModelArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final String MEAN_BIAS_STANDARD_DEVIATION_LONG_NAME = "mean-bias-standard-deviation";
    public static final String MAPPING_ERROR_RATE_LONG_NAME = "mapping-error-rate";
    public static final String GLOBAL_PSI_SCALE_LONG_NAME = "global-psi-scale";
    public static final String SAMPLE_PSI_SCALE_LONG_NAME = "sample-psi-scale";

    @Argument(
            doc = "Prior standard deviation of the contig-level mean coverage bias.  If a single sample is provided, " +
                    "this input will be ignored.",
            fullName = MEAN_BIAS_STANDARD_DEVIATION_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    private double meanBiasStandardDeviation = 0.01;

    @Argument(
            doc = "Typical mapping error rate.",
            fullName = MAPPING_ERROR_RATE_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    private double mappingErrorRate = 0.01;

    @Argument(
            doc = "Prior scale of contig coverage unexplained variance.  If a single sample is provided, " +
                    "this input will be ignored.",
            fullName = GLOBAL_PSI_SCALE_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    private double globalPsiScale = 0.001;

    @Argument(
            doc = "Prior scale of the sample-specific correction to the coverage unexplained variance.",
            fullName = SAMPLE_PSI_SCALE_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    private double samplePsiScale = 0.0001;

    public List<String> generatePythonArguments(final DetermineGermlineContigPloidy.RunMode runMode) {
        final List<String> arguments = new ArrayList<>(Arrays.asList(
                String.format("--mapping_error_rate=%e", mappingErrorRate),
                String.format("--psi_s_scale=%e", samplePsiScale)));
        if (runMode == DetermineGermlineContigPloidy.RunMode.COHORT) {
            arguments.addAll(Arrays.asList(
                    String.format("--mean_bias_sd=%e", meanBiasStandardDeviation),
                    String.format("--psi_j_scale=%e", globalPsiScale)));
        }
        return arguments;
    }

    public void validate() {
        ParamUtils.isPositive(meanBiasStandardDeviation,
                "Prior standard deviation of the contig-level mean coverage bias must be positive.");
        ParamUtils.isPositive(mappingErrorRate,
                "Typical mapping error rate must be positive.");
        ParamUtils.isPositive(globalPsiScale,
                "Prior scale of contig coverage unexplained variance must be positive.");
        ParamUtils.isPositive(samplePsiScale,
                "Prior scale of the sample-specific correction to the coverage unexplained variance " +
                        "must be positive.");
    }
}
