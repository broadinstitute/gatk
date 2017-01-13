package org.broadinstitute.hellbender.tools.coveragemodel.math;

import org.apache.commons.math3.analysis.solvers.*;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.function.Function;

/**
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public enum UnivariateSolverType {
    BRENT,
    BISECTION,
    MULLER,
    MULLER_2,
    ILLINOIS;

    private final Function<UnivariateSolverSpecifications, AbstractUnivariateSolver> solverFactory;

    UnivariateSolverType() {
        switch (name()) {
            case "BRENT":
                solverFactory = job -> new BrentSolver(job.getRelativeAccuracy(), job.getAbsoluteAccuracy(),
                        job.getFunctionValueAccuracy());
                break;
            case "BISECTION":
                solverFactory = job -> new BisectionSolver(job.getRelativeAccuracy(), job.getAbsoluteAccuracy());
                break;
            case "MULLER":
                solverFactory = job -> new MullerSolver(job.getRelativeAccuracy(), job.getAbsoluteAccuracy());
                break;
            case "MULLER_2":
                solverFactory = job -> new MullerSolver2(job.getRelativeAccuracy(), job.getAbsoluteAccuracy());
                break;
            case "ILLINOIS":
                solverFactory = job -> new IllinoisSolver(job.getRelativeAccuracy(), job.getAbsoluteAccuracy(),
                        job.getFunctionValueAccuracy());
                break;
            default:
                throw new GATKException.ShouldNeverReachHereException("Can not create the solver factory");
        }
    }

    public Function<UnivariateSolverSpecifications, AbstractUnivariateSolver> getSolverFactory() {
        return solverFactory;
    }
}