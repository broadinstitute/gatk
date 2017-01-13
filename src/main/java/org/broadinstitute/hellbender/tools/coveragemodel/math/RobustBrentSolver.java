package org.broadinstitute.hellbender.tools.coveragemodel.math;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.AbstractUnivariateSolver;
import org.apache.commons.math3.analysis.solvers.BaseAbstractUnivariateSolver;
import org.apache.commons.math3.exception.NoBracketingException;
import org.apache.commons.math3.exception.NumberIsTooLargeException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Precision;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A robust version of the Brent solver that tries to avoid spurious non-bracketing conditions.
 *
 * The root is bracketed by combing the solution search interval. If multiple roots are found,
 * the one that maximizes a provided "merit" function is chosen. If a merit function is not provided,
 * then root is chosen based on a merit "policy" (i.e. largest root, smallest root), in which case,
 * it is mandatory.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class RobustBrentSolver {

    /**
     * If a merit function is not provided, these are the possible policies
     */
    public enum MeritPolicy {
        SMALLEST_ROOT, LARGEST_ROOT
    }

    private final double relativeAccuracy, absoluteAccuracy, functionValueAccuracy;
    private int evaluations = 0;

    public RobustBrentSolver(final double relativeAccuracy, final double absoluteAccuracy,
                             final double functionValueAccuracy) {
        this.relativeAccuracy = relativeAccuracy;
        this.absoluteAccuracy = absoluteAccuracy;
        this.functionValueAccuracy = functionValueAccuracy;
    }

    /**
     * The main solver routine
     *
     * @param maxEval maximum function evaluations
     * @param objFunc the objective function to perform root finding on
     * @param meritFunc a merit function to choose the best root (if multiple roots are found)
     * @param min lower bound for the root
     * @param max upper bound for the root
     * @param numBisections number of bisections for "combing" the search interval
     * @param depth the fractal depth of the search grid
     * @return the root
     * @throws TooManyEvaluationsException if too many calls are made to the objective function
     * @throws NoBracketingException if the root could not be ultimately bracketed
     */
    public double solve(final int maxEval,
                        @Nonnull final UnivariateFunction objFunc,
                        @Nullable final UnivariateFunction meritFunc,
                        @Nullable final MeritPolicy meritPolicy,
                        final double min, final double max,
                        final int numBisections, final int depth)
            throws TooManyEvaluationsException, NoBracketingException {
        Utils.validateArg(meritFunc != null || meritPolicy != null, "Either a merit function or a merit policy" +
                        " must be provided");
        Utils.nonNull(objFunc, "The objective function must be non-null");
        ParamUtils.isPositive(maxEval, "Maximum function evaluations must be positive");
        ParamUtils.isPositiveOrZero(max - min, "The search interval is ill-posed");
        ParamUtils.isPositiveOrZero(numBisections, "Number of search interval bisections must be non-negative");
        ParamUtils.isPositiveOrZero(depth, "The search depth must be non-negative");

        evaluations = 0;

        /* create a search grid */
        final double[] xSearchGrid = createFractalSearchGrid(min, max, numBisections, depth);

        /* evaluate the function on the search grid */
        final double[] fSearchGrid = new double[xSearchGrid.length];
        for (int k = 0; k < xSearchGrid.length; k++) {
            fSearchGrid[k] = objFunc.value(xSearchGrid[k]);
            evaluations++;
        }

        final List<Bracket> bracketsList = detectBrackets(xSearchGrid, fSearchGrid);
        if (bracketsList.isEmpty()) {
            throw new NoBracketingException(min, max, fSearchGrid[0], fSearchGrid[fSearchGrid.length-1]);
        } else {
            final List<Double> roots = new ArrayList<>(bracketsList.size());
            for (final Bracket bracket : bracketsList) {
                final BrentSolverWithInitialFuncEvals solver = new BrentSolverWithInitialFuncEvals(relativeAccuracy,
                        absoluteAccuracy, functionValueAccuracy);
                try {
                    final double x0 = 0.5 * (bracket.xMin + bracket.xMax);
                    final double x = solver.solve(maxEval, objFunc, bracket.xMin, bracket.xMax, x0,
                            bracket.fMin, bracket.fMax, objFunc.value(x0));
                    roots.add(x);
                } catch (final TooManyEvaluationsException | NoBracketingException ex) {
                    /* nothing to do */
                }
                evaluations += solver.getEvaluations();
            }
            if (roots.isEmpty()) {
                throw new NoBracketingException(min, max, fSearchGrid[0], fSearchGrid[fSearchGrid.length-1]);
            } else if (roots.size() == 1) {
                return roots.get(0);
            } else {
                /* calculate merits */
                if (meritFunc != null) {
                    final List<Double> merits = roots.stream().mapToDouble(meritFunc::value).boxed().collect(Collectors.toList());
                    evaluations += merits.size();
                    double maxMerit = merits.get(0);
                    int maxMeritIndex = 0;
                    for (int rootIndex = 1; rootIndex < roots.size(); rootIndex++) {
                        if (merits.get(rootIndex) > maxMerit) {
                            maxMerit = merits.get(rootIndex);
                            maxMeritIndex = rootIndex;
                        }
                    }
                    return roots.get(maxMeritIndex);
                } else {
                    switch (meritPolicy) {
                        case SMALLEST_ROOT:
                            return Collections.min(roots);
                        case LARGEST_ROOT:
                            return Collections.max(roots);
                        default:
                            return Double.NaN;
                    }
                }
            }
        }
    }

    public int getEvaluations() {
        return evaluations;
    }

    @VisibleForTesting
    double[] createFractalSearchGrid(final double min, final double max, final int numBisections,
                                     final int refinementOrder) {
        final double[] baseGrid = createPowerOfTwoGrid(min, max, numBisections);
        final int nBaseGrid = numBisections + 2;
        if (refinementOrder > 1) {
            final double[] refinedGrid = new double[(nBaseGrid - 1) * refinementOrder + 1];
            for (int j = 0; j < nBaseGrid - 1; j++) {
                final double len = (baseGrid[j + 1] - baseGrid[j]) / refinementOrder;
                for (int k = 0; k < refinementOrder; k++) {
                    refinedGrid[j * refinementOrder + k] = baseGrid[j] + k * len;
                }
            }
            refinedGrid[(nBaseGrid - 1) * refinementOrder] = max;
            return refinedGrid;
        } else {
            return baseGrid;
        }
    }

    private double[] createPowerOfTwoGrid(final double min, final double max, final int numBisections) {
        final double[] grid = new double[numBisections + 2];
        grid[0] = 0;
        grid[numBisections + 1] = max - min;
        for (int j = numBisections; j > 0; j--) {
            grid[j] = grid[j+1]/2;
        }
        for (int j = 0; j < numBisections + 2; j++) {
            grid[j] += min;
        }
        return grid;
    }

    @VisibleForTesting
    static List<Bracket> detectBrackets(final double[] x, final double[] f) {
        final List<Bracket> brackets = new ArrayList<>();
        final double[] signs = new double[f.length];
        for (int i = 0; i < f.length; i++) {
            signs[i] = FastMath.signum(f[i]);
        }
        double prevSignum = signs[0];
        int prevIdx = 0;
        int idx = 1;
        while (idx < f.length) {
            if (signs[idx]*prevSignum <= 0) {
                brackets.add(new Bracket(x[prevIdx], x[idx], f[prevIdx], f[idx]));
                prevIdx = idx;
                prevSignum = signs[idx];
            }
            idx++;
        }
        return brackets;
    }

    @VisibleForTesting
    static final class Bracket {
        final double xMin, xMax, fMin, fMax;

        Bracket(final double xMin, final double xMax, final double fMin, final double fMax) {
            this.xMin = xMin;
            this.xMax = xMax;
            this.fMin = fMin;
            this.fMax = fMax;
        }

        @Override
        public String toString() {
            return "Bracket{" +
                    "xMin=" + xMin +
                    ", xMax=" + xMax +
                    ", fMin=" + fMin +
                    ", fMax=" + fMax +
                    '}';
        }
    }

    static private class BrentSolverWithInitialFuncEvals extends AbstractUnivariateSolver {
        /**
         * Construct a solver.
         *
         * @param relativeAccuracy Relative accuracy.
         * @param absoluteAccuracy Absolute accuracy.
         * @param functionValueAccuracy Function value accuracy.
         *
         * @see BaseAbstractUnivariateSolver#BaseAbstractUnivariateSolver(double,double,double)
         */
        public BrentSolverWithInitialFuncEvals(double relativeAccuracy,
                                               double absoluteAccuracy,
                                               double functionValueAccuracy) {
            super(relativeAccuracy, absoluteAccuracy, functionValueAccuracy);
        }

        public double solve(int maxEval, UnivariateFunction f, double min, double max, double startValue,
                            double fMin, double fMax, double fStartValue)
                throws TooManyEvaluationsException, NoBracketingException {
            setup(maxEval, f, min, max, startValue);
            return doSolve(fMin, fMax, fStartValue);
        }

        /**
         * {@inheritDoc}
         */
        @Override
        protected double doSolve()
                throws NoBracketingException,
                TooManyEvaluationsException,
                NumberIsTooLargeException {
            double min = getMin();
            double max = getMax();
            final double initial = getStartValue();
            final double functionValueAccuracy = getFunctionValueAccuracy();

            verifySequence(min, initial, max);

            // Return the initial guess if it is good enough.
            double yInitial = computeObjectiveValue(initial);
            if (FastMath.abs(yInitial) <= functionValueAccuracy) {
                return initial;
            }

            // Return the first endpoint if it is good enough.
            double yMin = computeObjectiveValue(min);
            if (FastMath.abs(yMin) <= functionValueAccuracy) {
                return min;
            }

            // Reduce interval if min and initial bracket the root.
            if (yInitial * yMin < 0) {
                return brent(min, initial, yMin, yInitial);
            }

            // Return the second endpoint if it is good enough.
            double yMax = computeObjectiveValue(max);
            if (FastMath.abs(yMax) <= functionValueAccuracy) {
                return max;
            }

            // Reduce interval if initial and max bracket the root.
            if (yInitial * yMax < 0) {
                return brent(initial, max, yInitial, yMax);
            }

            throw new NoBracketingException(min, max, yMin, yMax);
        }

        private double doSolve(final double fMin, final double fMax, final double fStartValue)
                throws NoBracketingException,
                TooManyEvaluationsException,
                NumberIsTooLargeException {
            double min = getMin();
            double max = getMax();
            final double initial = getStartValue();
            final double functionValueAccuracy = getFunctionValueAccuracy();

            verifySequence(min, initial, max);

            // Return the initial guess if it is good enough.
            if (FastMath.abs(fStartValue) <= functionValueAccuracy) {
                return initial;
            }

            // Return the first endpoint if it is good enough.
            if (FastMath.abs(fMin) <= functionValueAccuracy) {
                return min;
            }

            // Reduce interval if min and initial bracket the root.
            if (fStartValue * fMin < 0) {
                return brent(min, initial, fMin, fStartValue);
            }

            // Return the second endpoint if it is good enough.
            if (FastMath.abs(fMax) <= functionValueAccuracy) {
                return max;
            }

            // Reduce interval if initial and max bracket the root.
            if (fStartValue * fMax < 0) {
                return brent(initial, max, fStartValue, fMax);
            }

            throw new NoBracketingException(min, max, fMin, fMax);
        }

        /**
         * Search for a zero inside the provided interval.
         * This implementation is based on the algorithm described at page 58 of
         * the book
         * <blockquote>
         *  <b>Algorithms for Minimization Without Derivatives</b>
         *  <it>Richard P. Brent</it>
         *  Dover 0-486-41998-3
         * </blockquote>
         *
         * @param lo Lower bound of the search interval.
         * @param hi Higher bound of the search interval.
         * @param fLo Function value at the lower bound of the search interval.
         * @param fHi Function value at the higher bound of the search interval.
         * @return the value where the function is zero.
         */
        private double brent(double lo, double hi,
                             double fLo, double fHi) {
            double a = lo;
            double fa = fLo;
            double b = hi;
            double fb = fHi;
            double c = a;
            double fc = fa;
            double d = b - a;
            double e = d;

            final double t = getAbsoluteAccuracy();
            final double eps = getRelativeAccuracy();

            while (true) {
                if (FastMath.abs(fc) < FastMath.abs(fb)) {
                    a = b;
                    b = c;
                    c = a;
                    fa = fb;
                    fb = fc;
                    fc = fa;
                }

                final double tol = 2 * eps * FastMath.abs(b) + t;
                final double m = 0.5 * (c - b);

                if (FastMath.abs(m) <= tol ||
                        Precision.equals(fb, 0))  {
                    return b;
                }
                if (FastMath.abs(e) < tol ||
                        FastMath.abs(fa) <= FastMath.abs(fb)) {
                    // Force bisection.
                    d = m;
                    e = d;
                } else {
                    double s = fb / fa;
                    double p;
                    double q;
                    // The equality test (a == c) is intentional,
                    // it is part of the original Brent's method and
                    // it should NOT be replaced by proximity test.
                    if (a == c) {
                        // Linear interpolation.
                        p = 2 * m * s;
                        q = 1 - s;
                    } else {
                        // Inverse quadratic interpolation.
                        q = fa / fc;
                        final double r = fb / fc;
                        p = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
                        q = (q - 1) * (r - 1) * (s - 1);
                    }
                    if (p > 0) {
                        q = -q;
                    } else {
                        p = -p;
                    }
                    s = e;
                    e = d;
                    if (p >= 1.5 * m * q - FastMath.abs(tol * q) ||
                            p >= FastMath.abs(0.5 * s * q)) {
                        // Inverse quadratic interpolation gives a value
                        // in the wrong direction, or progress is slow.
                        // Fall back to bisection.
                        d = m;
                        e = d;
                    } else {
                        d = p / q;
                    }
                }
                a = b;
                fa = fb;

                if (FastMath.abs(d) > tol) {
                    b += d;
                } else if (m > 0) {
                    b += tol;
                } else {
                    b -= tol;
                }
                fb = computeObjectiveValue(b);
                if ((fb > 0 && fc > 0) ||
                        (fb <= 0 && fc <= 0)) {
                    c = a;
                    fc = fa;
                    d = b - a;
                    e = d;
                }
            }
        }
    }
}
