package org.broadinstitute.hellbender.tools.coveragemodel.math;

/**
 * This class stores the description of a solver job
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */

public final class UnivariateSolverJobDescription {

    private final int index, maxEvaluations;
    private final double min, max, x0;

    public UnivariateSolverJobDescription(final int index, final double min, final double max, final double x0,
                                          final int maxEvaluations) {
        this.index = index;
        this.min = min;
        this.max = max;
        this.x0 = x0;
        this.maxEvaluations = maxEvaluations;
    }

    public int getIndex() {
        return index;
    }

    public int getMaxEvaluations() {
        return maxEvaluations;
    }

    public double getMin() {
        return min;
    }

    public double getMax() {
        return max;
    }

    public double getInitialGuess() {
        return x0;
    }

}