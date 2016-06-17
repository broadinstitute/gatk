package org.broadinstitute.hellbender.tools.exome.copyratio;

import org.broadinstitute.hellbender.utils.mcmc.ParameterEnum;

/**
 * Enumerates the parameters for {@link CopyRatioState}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public enum CopyRatioParameter implements ParameterEnum {
    VARIANCE("CR_variance"),
    OUTLIER_PROBABILITY("CR_outlier_probability"),
    SEGMENT_MEANS("CR_segment_means"),
    OUTLIER_INDICATORS("CR_outlier_indicators");

    public final String name;

    CopyRatioParameter(final String name) {
        this.name = name;
    }
}
