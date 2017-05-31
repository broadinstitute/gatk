package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;

import java.util.function.Function;

/**
 * Created by davidben on 4/27/16.
 */
public final class OptimizationUtils {
    protected static final MaxEval DEFAULT_MAX_EVAL = new MaxEval(1000);
    private static final double DEFAULT_RELATIVE_TOLERANCE = 0.001;
    private static final double DEFAULT_ABSOLUTE_TOLERANCE = 0.001;
    protected static final BrentOptimizer DEFAULT_OPTIMIZER = new BrentOptimizer(DEFAULT_RELATIVE_TOLERANCE, DEFAULT_ABSOLUTE_TOLERANCE);;

    private static final double DEFAULT_EPSILON_FOR_NUMERIC_DERIVATIVES = 1.0e-5;
    private static final double DEFAULT_DERIVATIVE_THRESHOLD = 1.0e-8;

    private OptimizationUtils() { }

    public static double argmax(final Function<Double, Double> function, final double min, final double max, final double guess) {
        final SearchInterval interval = new SearchInterval(min, max, guess);
        return DEFAULT_OPTIMIZER.optimize(new UnivariateObjectiveFunction(function::apply), GoalType.MAXIMIZE, interval, DEFAULT_MAX_EVAL).getPoint();
    }

    public static double argmax(final Function<Double, Double> function, final double min, final double max, final double guess,
                                final double relativeTolerance, final double absoluteTolerance, final int maxEvaluations) {
        final BrentOptimizer optimizer = new BrentOptimizer(relativeTolerance, absoluteTolerance);
        final SearchInterval interval = new SearchInterval(min, max, guess);
        return optimizer.optimize(new UnivariateObjectiveFunction(function::apply), GoalType.MAXIMIZE, interval, new MaxEval(maxEvaluations)).getPoint();
    }


    /**
     * One iteration of Newton's method for univariate optimization (i.e. not repeated until convergence)
     * @param function  the function to maximize
     * @param min       minimum allowable value
     * @param max       maximum allowable value
     * @param guess     current guess
     * @param epsilon   infinitesimal for numerical first and second derivatives
     * @param derivativeThreshold   first derivative below which current guess is considered stationary
     */
    public static double singleNewtonArgmaxUpdate(final Function<Double, Double> function, final double min, final double max,
                                                  final double guess, final double epsilon, final double derivativeThreshold) {
        final double fx = function.apply(guess);
        final double fxPlusEps = function.apply(guess + epsilon);
        final double fxMinusEps = function.apply(guess - epsilon);
        final double secondDerivative = (fxPlusEps + fxMinusEps - 2 * fx) / (epsilon * epsilon);
        final double firstDerivative = (fxPlusEps - fxMinusEps)/(2*epsilon);
        final double argmax = guess - firstDerivative / secondDerivative;
        if (Math.abs(firstDerivative) < derivativeThreshold) {
            return guess;
        }
        return argmax < min ? min : (argmax > max ? max : argmax);
    }

    public static double singleNewtonArgmaxUpdate(final Function<Double, Double> function, final double min, final double max, final double guess) {
        return singleNewtonArgmaxUpdate(function, min, max, guess, DEFAULT_EPSILON_FOR_NUMERIC_DERIVATIVES, DEFAULT_DERIVATIVE_THRESHOLD);
    }
}
