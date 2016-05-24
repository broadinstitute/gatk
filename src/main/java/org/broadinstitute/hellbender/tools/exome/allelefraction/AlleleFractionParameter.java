package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.broadinstitute.hellbender.utils.mcmc.ParameterEnum;

/**
 * Enumerates the parameters for {@link AlleleFractionState}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public enum AlleleFractionParameter implements ParameterEnum {
    MEAN_BIAS, BIAS_VARIANCE, OUTLIER_PROBABILITY, MINOR_ALLELE_FRACTIONS
}
