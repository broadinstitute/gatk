package org.broadinstitute.hellbender.tools.coveragemodel.math;

import org.apache.commons.math3.exception.TooManyEvaluationsException;

/**
 * This class stores the description of a solver job.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */

public final class UnivariateSolverJobDescription {

    /**
     * Job index (an arbitrary task-specific integer index)
     */
    private final int index;

    /**
     * Maximum function evaluations; solvers will throw a {@link TooManyEvaluationsException} exception if function
     * evaluations exceed this number
     */
    private final int maxEvaluations;

    /**
     * Left endpoint for an interval that brackets the root
     */
    private final double min;

    /**
     * Right endpoint for an interval that brackets the root
     */
    private final double max;

    /**
     * Initial guess (must lie inside the interval)
     */
    private final double x0;

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