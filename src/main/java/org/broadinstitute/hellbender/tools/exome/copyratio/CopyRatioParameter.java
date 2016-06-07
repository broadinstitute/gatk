package org.broadinstitute.hellbender.tools.exome.copyratio;

import org.broadinstitute.hellbender.utils.mcmc.ParameterEnum;

/**
 * Enumerates the parameters for {@link CopyRatioState}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public enum CopyRatioParameter implements ParameterEnum {
    VARIANCE, OUTLIER_PROBABILITY, SEGMENT_MEANS, OUTLIER_INDICATORS
}
