package org.broadinstitute.hellbender.tools.copynumber.arguments;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.tools.copynumber.DetermineGermlineContigPloidy;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class GermlineContigPloidyModelArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final String MEAN_BIAS_STANDARD_DEVIATION_LONG_NAME = "mean-bias-standard-deviation";
    public static final String MAPPING_ERROR_RATE_LONG_NAME = "mapping-error-rate";

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

    public List<String> generatePythonArguments(final DetermineGermlineContigPloidy.RunMode runMode) {
        final List<String> arguments = new ArrayList<>(Collections.singletonList(
                String.format("--mapping_error_rate=%e", mappingErrorRate)));
        if (runMode == DetermineGermlineContigPloidy.RunMode.COHORT) {
            arguments.addAll(Collections.singletonList(
                    String.format("--mean_bias_sd=%e", meanBiasStandardDeviation)));
        }
        return arguments;
    }

    public void validate() {
        ParamUtils.isPositive(meanBiasStandardDeviation,
                "Prior standard deviation of the contig-level mean coverage bias must be positive.");
        ParamUtils.isPositive(mappingErrorRate,
                "Typical mapping error rate must be positive.");
    }
}
