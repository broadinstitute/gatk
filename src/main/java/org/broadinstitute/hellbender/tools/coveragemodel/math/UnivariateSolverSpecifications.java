package org.broadinstitute.hellbender.tools.coveragemodel.math;

/**
 * Accuracy specifications of a univariate solver
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class UnivariateSolverSpecifications {

    private final double absoluteAccuracy, relativeAccuracy, functionValueAccuracy;

    public UnivariateSolverSpecifications(final double absoluteAccuracy, final double relativeAccuracy,
                                          final double functionValueAccuracy) {
        this.absoluteAccuracy = absoluteAccuracy;
        this.relativeAccuracy = relativeAccuracy;
        this.functionValueAccuracy = functionValueAccuracy;
    }

    public double getAbsoluteAccuracy() {
        return absoluteAccuracy;
    }

    public double getRelativeAccuracy() {
        return relativeAccuracy;
    }

    public double getFunctionValueAccuracy() {
        return functionValueAccuracy;
    }
}
