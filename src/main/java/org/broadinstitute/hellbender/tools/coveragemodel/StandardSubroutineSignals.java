package org.broadinstitute.hellbender.tools.coveragemodel;

/**
 * This class provides a number of standard keys for setting and getting values into and out of a
 * {@link SubroutineSignal}. The value types are context-dependent.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public enum StandardSubroutineSignals {
    /**
     * Number of iterations until termination
     */
    ITERATIONS,

    /**
     * Residual error norm (e.g. in iterative solvers)
     */
    RESIDUAL_ERROR_NORM,

    /**
     * The final solution (in case a routine communicates its result via a {@link SubroutineSignal} instance)
     */
    SOLUTION,

    /**
     * Exit status or termination criterion of a routine
     */
    EXIT_STATUS
}
