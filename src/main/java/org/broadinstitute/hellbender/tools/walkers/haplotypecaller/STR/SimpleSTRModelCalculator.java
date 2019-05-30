package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.STR;

/**
 * Created by valentin on 11/29/16.
 */
public  class SimpleSTRModelCalculator implements  STRModelCalculator {

    private final STRModelParameter pi;
    private final STRModelParameter tau;
    private final STRModelParameter del;
    private final STRModelParameter ins;

    SimpleSTRModelCalculator(
            final STRModelParameter pi,
            final STRModelParameter tau,
            final STRModelParameter del,
            final STRModelParameter ins) {
        this.pi = pi;
        this.tau = tau;
        this.del = del;
        this.ins = ins;
    }

    public double logCoefficient(final int actualReapeatCount, final int observedRepeatCount) {
        if (actualReapeatCount == observedRepeatCount) {
            return Math.log1p(- pi.valueFor(actualReapeatCount));
        } else {
            final double piValue = pi.valueFor(actualReapeatCount);
            if (actualReapeatCount > observedRepeatCount) {
                final double delValue = del.valueFor(actualReapeatCount);
                return Math.log(piValue * tau.valueFor(observedRepeatCount))
                        + (actualReapeatCount - observedRepeatCount - 1) * Math.log(delValue) + Math.log1p(-delValue);
            } else {
                final double insValue = ins.valueFor(observedRepeatCount);
                return Math.log(piValue * (1 - tau.valueFor(observedRepeatCount))) + (observedRepeatCount - actualReapeatCount - 1) * Math.log(insValue)
                        +  Math.log1p(-insValue);
            }
        }
    }
}
