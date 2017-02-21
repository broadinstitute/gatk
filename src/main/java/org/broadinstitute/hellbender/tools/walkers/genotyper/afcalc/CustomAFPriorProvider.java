package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;

/**
 * Custom priors provided as an array.
 */
public final class CustomAFPriorProvider extends AFPriorProvider {

    private final double[] priors;

    /**
     *
     * @param priorValues must exactly equal the number of genomes in the samples (the total ploidy).
     */
    public CustomAFPriorProvider(final List<Double> priorValues) {
        Utils.nonNull(priorValues,"the input prior values cannot be null");
        priors = new double[priorValues.size() + 1];

        int i = 1;
        double sum = 0;
        for (final double value : priorValues) {
            if (value <= 0 || value >= 1) {
                throw new IllegalArgumentException("the AF prior value " + value + " is out of the valid interval (0,1)");
            }
            if (Double.isNaN(value)) {
                throw new IllegalArgumentException("NaN is not a valid prior AF value");
            }
            priors[i++] = Math.log10(value);
            sum += value;
        }
        if (sum >= 1) {
            throw new IllegalArgumentException("the AF prior value sum must be less than 1: " + sum);
        }
        priors[0] = Math.log10(1 - sum);
    }

    @Override
    protected double[] buildPriors(final int totalPloidy) {
        return priors;
    }
}
