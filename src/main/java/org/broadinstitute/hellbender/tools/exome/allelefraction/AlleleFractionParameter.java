package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.broadinstitute.hellbender.utils.mcmc.ParameterEnum;

/**
 * Enumerates the parameters for {@link AlleleFractionState}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public enum AlleleFractionParameter implements ParameterEnum {
    MEAN_BIAS("AF_reference_bias_mean"),
    BIAS_VARIANCE("AF_reference_bias_variance"),
    OUTLIER_PROBABILITY("AF_outlier_probability"),
    MINOR_ALLELE_FRACTIONS("AF_minor_allele_fractions");

    public final String name;

    AlleleFractionParameter(final String name) {
        this.name = name;
    }
}
