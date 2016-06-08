package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AFPriorProvider;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.Arrays;

/**
 * Allele frequency prior provider based on heterozygosity.
 */
public final class HeterozygosityAFPriorProvider extends AFPriorProvider {

    private final double heterozygosity;
    private final double log10Heterozygosity;

    /**
     * Construct a new provider given the heterozygosity value.
     * @param heterozygosity must be a valid heterozygosity between larger than 0 and smaller than 1.
     * @throws IllegalArgumentException if {@code heterozygosity} is not valid one in the interval (0,1).
     */
    public HeterozygosityAFPriorProvider(final double heterozygosity) {
        if (heterozygosity <= 0) {
            throw new IllegalArgumentException("the heterozygosity must be greater than 0");
        }
        if (heterozygosity >= 1) {
            throw new IllegalArgumentException("the heterozygosity must be less than 1");
        }
        if (Double.isNaN(heterozygosity)) {
            throw new IllegalArgumentException("the heterozygosity cannot be a NaN");
        }
        this.heterozygosity = heterozygosity;
        this.log10Heterozygosity = Math.log10(heterozygosity);
    }

    @Override
    protected double[] buildPriors(final int totalPloidy) {
        final double[] result = new double [totalPloidy + 1];
        Arrays.fill(result, log10Heterozygosity);
        result[0] = Double.NEGATIVE_INFINITY;
        for (int i = 1; i <= totalPloidy; i++) {
            result[i] -= MathUtils.log10(i);
        }
        final double log10Sum = MathUtils.approximateLog10SumLog10(result);
        if (log10Sum >= 0) {
            throw new IllegalArgumentException("heterozygosity " + heterozygosity + " is too large of total ploidy " + totalPloidy);
        }
        result[0] = MathUtils.log10OneMinusPow10(log10Sum);
        return result;
    }
}
