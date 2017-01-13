package org.broadinstitute.hellbender.tools.coveragemodel.linalg;

import org.apache.commons.math3.linear.NonPositiveDefiniteOperatorException;
import org.broadinstitute.hellbender.tools.coveragemodel.SubroutineSignal;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.util.function.BiFunction;
import java.util.function.Function;

/**
 * A basic implementation for iterative linear system solvers
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
    private final double absTol, relTol;
    private final int maxIters;

    public enum ExitStatus {
        SUCCESS_ABS_TOL,
        SUCCESS_REL_TOL,
        FAIL_MAX_ITERS
    }

    /**
     * Constructor
     *
     * @param linop a linear operator (a rank-2 matrix or an abstract operator)
     * @param b rhs vector (a rank-1 object)
     * @param precond the preconditioner
     * @param absTol absolute tolerance
     * @param relTol relative tolerance
     * @param maxIters maximum iterations
     * @param normFunc a norm function
     * @param innerProductFunc an inner product function
     * @param performDetailedChecks perform detailed checks
     */
    public IterativeLinearSolverNDArray(@Nonnull final GeneralLinearOperator<INDArray> linop,
                                        @Nonnull final INDArray b,
                                        @Nullable final GeneralLinearOperator<INDArray> precond,
                                        final double absTol, final double relTol,
                                        final int maxIters,
                                        @Nonnull final Function<INDArray, Double> normFunc,
                                        @Nonnull final BiFunction<INDArray, INDArray, Double> innerProductFunc,
                                        final boolean performDetailedChecks) {
        this.performDetailedChecks = performDetailedChecks;
        this.linop = linop;
        if (b.length() != linop.getColumnDimension()) {
            throw new IllegalArgumentException("The provided rhs vector has wrong dimensions");
        } else {
            this.b = b;
        }

        this.precond = precond;
        if (precond != null) {
            if (linop.getColumnDimension() != precond.getColumnDimension() ||
                    linop.getRowDimension() != precond.getRowDimension()) {
                throw new IllegalArgumentException("The preconditioner and linear operators must have the same " +
                        "row and column dimensions");
            }
        }

        this.normFunc = normFunc;
        this.innerProductFunc = innerProductFunc;

        this.absTol = ParamUtils.isPositive(absTol, "Absolute error tolerance must be positive");
        this.relTol = ParamUtils.isPositive(relTol, "Relative error tolerance must be pisitive");
        this.maxIters = ParamUtils.isPositive(maxIters, "Maximum iterations must be positive");
    }

    /**
     * Preconditioned conjugate gradients method
     *
     * @param x0 initial guess
     * @return an instance of {@link SubroutineSignal}
     * @throws IllegalArgumentException for bad intial guess or illegla linear operators
     */
    public SubroutineSignal cg(@Nonnull final INDArray x0)
            throws IllegalArgumentException {
        if (x0.length() != linop.getRowDimension()) {
            throw new IllegalArgumentException("The initial guess has wrong dimensions");
        }
        if (linop.getRowDimension() != linop.getColumnDimension()) {
            throw new UnsupportedOperationException("The CG routine only works for square linear operators");
        }

        /* initialize */
        int iter = 0;
        final double rmax = relTol * normFunc.apply(b);
        final INDArray x = x0.dup();
        final INDArray p = x.dup();
        final INDArray q = linop.operate(p);
        final INDArray r = b.sub(q);
        final INDArray z;
        if (precond == null) {
            z = r; /* by reference */
        } else {
            z = Nd4j.create(r.shape());
        }
        double rnorm = normFunc.apply(r);
        if (rnorm <= rmax) {
            return generateCGSignal(x, iter, rnorm, rmax, ExitStatus.SUCCESS_REL_TOL);
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
                return generateCGSignal(x, iter, rnorm, rmax, ExitStatus.SUCCESS_REL_TOL);
            }
            if (rnorm <= absTol) {
                return generateCGSignal(x, iter, rnorm, rmax, ExitStatus.SUCCESS_ABS_TOL);
            }
            iter += 1;
            if (iter >= maxIters) {
                return generateCGSignal(x, iter, rnorm, rmax, ExitStatus.FAIL_MAX_ITERS);
            }
        }
    }

    /**
     * Generate an exit signal for {@link #cg(INDArray)}
     *
     * @param x solution
     * @param iter number of iterations
     * @param resNorm residual error norm
     * @param resNormThreshold latest solver stopping criterion on residual error norm
     * @param status the exist status of the solver
     * @return a instance of {@link SubroutineSignal}
     */
    private static SubroutineSignal generateCGSignal(@Nonnull final INDArray x,
                                                    final int iter,
                                                    final double resNorm,
                                                    final double resNormThreshold,
                                                    @Nonnull final ExitStatus status) {
        return SubroutineSignal.builder()
                .put("x", x)
                .put("iterations", iter)
                .put("residual_norm", resNorm)
                .put("residual_norm_termination_threshold", resNormThreshold)
                .put("status", status).build();
    }
}
