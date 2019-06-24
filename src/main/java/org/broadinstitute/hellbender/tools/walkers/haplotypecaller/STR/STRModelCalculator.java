package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.STR;


/**
 * Created by valentin on 11/29/16.
 */
public interface STRModelCalculator {

    double logCoefficient(final int actualReapeatCount, final int observedRepeatCount);

    default double log10Coefficient(final int actualRepeatCount, final int observedRepeatCount) {
        final double INV_LN10 = 1.0 / Math.log(10);
        return logCoefficient(actualRepeatCount, observedRepeatCount) * INV_LN10;
    }
}
