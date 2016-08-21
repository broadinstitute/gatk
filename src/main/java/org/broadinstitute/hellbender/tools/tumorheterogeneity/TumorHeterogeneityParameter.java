package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.broadinstitute.hellbender.utils.mcmc.ParameterEnum;


/**
 * Enumerates the parameters for {@link TumorHeterogeneityState}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public enum TumorHeterogeneityParameter implements ParameterEnum {
    CONCENTRATION("concentration"),
    COPY_RATIO_NOISE_FLOOR("cr_noise_floor"),
    COPY_RATIO_NOISE_FACTOR("cr_noise_factor"),
    MINOR_ALLELE_FRACTION_NOISE_FACTOR("maf_noise_factor"),
    POPULATION_MIXTURE("population_mixture");

    public final String name;

    TumorHeterogeneityParameter(final String name) {
        this.name = name;
    }
}

