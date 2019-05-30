package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.STR;

/**
 * Created by valentin on 11/29/16.
 */
final class NullSTRModelCalculator implements STRModelCalculator {

    static final NullSTRModelCalculator INSTANCE = new NullSTRModelCalculator();

    @Override
    public double logCoefficient(int actualReapeatCount, int observedRepeatCount) {
        return actualReapeatCount == observedRepeatCount ? 0.0 : Double.NEGATIVE_INFINITY;
    }

    private NullSTRModelCalculator() {

    }
}
