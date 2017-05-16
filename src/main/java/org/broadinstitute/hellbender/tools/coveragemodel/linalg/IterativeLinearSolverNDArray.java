package org.broadinstitute.hellbender.tools.coveragemodel.linalg;

import org.apache.commons.math3.linear.NonPositiveDefiniteOperatorException;
import org.broadinstitute.hellbender.tools.coveragemodel.StandardSubroutineSignals;
import org.broadinstitute.hellbender.tools.coveragemodel.SubroutineSignal;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.util.function.BiFunction;
import java.util.function.Function;

/**
 * A basic implementation for iterative linear system solvers using Nd4j as the liner algebra backend.
 * The class currently provides the preconditioned conjugate gradient solver; see {@link #solveUsingPreconditionedConjugateGradient(INDArray)}.
 *
 * This class is useful for solving linear system of equations $Ax = b$ where $A$ is an instance of
 * {@link GeneralLinearOperator<INDArray>} and represents a symmetric positive-definite linear operator.
 *
 * Common uses cases are:
 *
 *   - when $A$ is _not_ specified by a dense block matrix representation (e.g. for efficiency or memory
 *     constraints), though, the linear operator can be applied to an arbitrary vector _algorithmically_.
 *
 *   - when a decent preconditioner exists (e.g. $A$ is roughly diagonal so that the inverse of the diagonal
 *     part can be used as a preconditioner)
 *
 *   - when a decent guess for the solution vector $x$ is already available, and one requires an approximate
 *     refinement within specified relative/absolute tolerance of the exact solution.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class IterativeLinearSolverNDArray {

    private final GeneralLinearOperator<INDArray> linop;
    private final GeneralLinearOperator<INDArray> precond;
    private final INDArray b;
    private final Function<INDArray, Double> normFunc;
    private final BiFunction<INDArray, INDArray, Double> innerProductFunc;
    private final boolean performDetailedChecks;
    private final double absTol;
    private final double relTol;
    private final int maxIters;

    public static final String RESIDUAL_NORM_TERMINATION_THRESHOLD_SUBROUTINE_SIGNAL_KEY =
            "RESIDUAL_NORM_TERMINATION_THRESHOLD";

    public enum ExitStatus {
        SUCCESS_ABS_TOL,
        SUCCESS_REL_TOL,
        FAIL_MAX_ITERS
    }

    /**
     * Constructs the iterative solver.
     *
     * @param linop a linear operator (a rank-2 matrix or an abstract operator)
     * @param b rhs vector (a rank-1 object)
     * @param precond the preconditioner
     * @param absTol absolute tolerance
     * @param relTol relative tolerance
     * @param maxIters maximum iterations
     * @param normFunc a norm function
     * @param innerProductFunc an inner product function
     * @param performDetailedChecks perform detailed checks (has different meaning in different solvers)
     *                              refer to the documentation of a specific solver
     */
    public IterativeLinearSolverNDArray(@Nonnull final GeneralLinearOperator<INDArray> linop,
                                        @Nonnull final INDArray b,
                                        @Nullable final GeneralLinearOperator<INDArray> precond,
                                        final double absTol, final double relTol,
                                        final int maxIters,
                                        @Nonnull final Function<INDArray, Double> normFunc,
                                        @Nonnull final BiFunction<INDArray, INDArray, Double> innerProductFunc,
                                        final boolean performDetailedChecks) {
        Utils.nonNull(b, "The right hand side vector must be non-null");
        Utils.nonNull(linop, "The linear operator must be non-null");
        Utils.nonNull(normFunc, "The norm function must be non-null");
        Utils.nonNull(innerProductFunc, "The inner product function must be non-null");
        Utils.validateArg(b.length() == linop.getColumnDimension(), "The provided rhs vector has wrong dimensions");
        Utils.validateArg(precond == null || (linop.getColumnDimension() == precond.getColumnDimension() &&
                linop.getRowDimension() == precond.getRowDimension()), "The preconditioner and linear operators" +
                " must have the same row and column dimensions");
        this.performDetailedChecks = performDetailedChecks;
        this.linop = linop;
        this.b = b;
        this.precond = precond;
        this.normFunc = normFunc;
        this.innerProductFunc = innerProductFunc;
        this.absTol = ParamUtils.isPositive(absTol, "Absolute error tolerance must be positive");
        this.relTol = ParamUtils.isPositive(relTol, "Relative error tolerance must be positive");
        this.maxIters = ParamUtils.isPositive(maxIters, "Maximum iterations must be positive");
    }

    /**
     * Solve using the preconditioned conjugate gradient method. The implementation is based on:
     *
     * <dt><a id="BARR1994">Barret et al. (1994)</a></dt>
     * <dd>R. Barrett, M. Berry, T. F. Chan, J. Demmel, J. M. Donato, J. Dongarra,
     * V. Eijkhout, R. Pozo, C. Romine and H. Van der Vorst,
     * <a href="http://www.netlib.org/linalg/html_templates/Templates.html"><em>
     * Templates for the Solution of Linear Systems: Building Blocks for Iterative
     * Methods</em></a>, SIAM</dd>
     *
     * If {@link #performDetailedChecks} is true, the (required) positive-definiteness condition
     * on the linear operator is checked during each iteration (with no computation cost) and an
     * exception is thrown if the condition is violated.
     *
     * The final solution is contained in the output {@link SubroutineSignal}, along with the norm
     * of the residual vector, exit status, etc; see {@link #generateConjugateGradientsSubroutineSignal}
     *
     * @param x0 initial guess
     * @return an instance of {@link SubroutineSignal} containing the solution and residual error
     * @throws IllegalArgumentException for bad initial guess or illegal linear operators
     * @throws NonPositiveDefiniteOperatorException if the linear operator is deemed to be non-positive-definite
     */
    public SubroutineSignal solveUsingPreconditionedConjugateGradient(@Nonnull final INDArray x0) throws IllegalArgumentException {
        Utils.nonNull(x0, "The initial guess must be non-null");
        Utils.validateArg(x0.length() == linop.getRowDimension(), "The initial guess has wrong dimensions");
        Utils.validateArg(linop.getRowDimension() == linop.getColumnDimension(), "The CG routine only works for square" +
                " linear operators");

        /* initialize */
        int iter = 0;
        final double rmax = relTol * normFunc.apply(b);
        final INDArray x = x0.dup();
        final INDArray p = x.dup();
        final INDArray q = linop.operate(p);
        final INDArray r = b.sub(q);
        final INDArray z;
        z = precond == null ? r : Nd4j.create(r.shape());
        double rnorm = normFunc.apply(r);
        if (rnorm <= rmax) {
            return generateConjugateGradientsSubroutineSignal(x, iter, rnorm, rmax, ExitStatus.SUCCESS_REL_TOL);
        }

        double rhoPrev = 0;
        while (true) {
            /* z_{i-1} = M r_{i-1} */
            if (precond != null) {
                z.assign(precond.operate(r));
            }
            /* rho_{i-1} = r_{i-1}^T z_{i-1} */
            final double rhoNext = innerProductFunc.apply(r, z);
            if (performDetailedChecks && rhoNext <= 0) {
                throw new NonPositiveDefiniteOperatorException();
            }
            if (iter == 0) {
                /* p_1 = z_0 */
                p.assign(z);
            } else {
                /* p_{i} = z_{i-1}  + (\rho_{i-1}/\rho_{i-2}) p_{i-1} */
                p.muli(rhoNext / rhoPrev).addi(z);
            }
            /* q_{i} = A p_{i} */
            q.assign(linop.operate(p));
            final double pq = innerProductFunc.apply(p, q);
            if (performDetailedChecks && pq <= 0) {
                throw new NonPositiveDefiniteOperatorException();
            }
            /* \alpha_{i} = \rho_{i-1}/(p_{i}^T q_{i}^T) */
            final double alpha = rhoNext / pq;
            /* x_{i} = x_{i-1}  + \alpha_{i} p_{i} */
            x.addi(p.mul(alpha));
            /* r_{i} = r_{i-1}  - \alpha_{i} q_{i} */
            r.addi(q.mul(-alpha));
            rhoPrev = rhoNext;
            rnorm = normFunc.apply(r);
            if (rnorm <= rmax) {
                return generateConjugateGradientsSubroutineSignal(x, iter, rnorm, rmax, ExitStatus.SUCCESS_REL_TOL);
            }
            if (rnorm <= absTol) {
                return generateConjugateGradientsSubroutineSignal(x, iter, rnorm, rmax, ExitStatus.SUCCESS_ABS_TOL);
            }
            if (++iter >= maxIters) {
                return generateConjugateGradientsSubroutineSignal(x, iter, rnorm, rmax, ExitStatus.FAIL_MAX_ITERS);
            }
        }
    }

    /**
     * Generate an exit signal for {@link #solveUsingPreconditionedConjugateGradient(INDArray)}
     *
     * @param x solution
     * @param iter number of iterations
     * @param resNorm residual error norm
     * @param resNormThreshold latest solver stopping criterion on residual error norm
     * @param status the exit status of the solver
     * @return a instance of {@link SubroutineSignal}
     */
    private static SubroutineSignal generateConjugateGradientsSubroutineSignal(@Nonnull final INDArray x,
                                                                               final int iter,
                                                                               final double resNorm,
                                                                               final double resNormThreshold,
                                                                               @Nonnull final ExitStatus status) {
        return SubroutineSignal.builder()
                .put(StandardSubroutineSignals.SOLUTION, x)
                .put(StandardSubroutineSignals.ITERATIONS, iter)
                .put(StandardSubroutineSignals.RESIDUAL_ERROR_NORM, resNorm)
                .put(RESIDUAL_NORM_TERMINATION_THRESHOLD_SUBROUTINE_SIGNAL_KEY, resNormThreshold)
                .put(StandardSubroutineSignals.EXIT_STATUS, status).build();
    }
}
