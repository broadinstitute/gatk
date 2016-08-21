package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.distribution.EnumeratedDistribution;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.MathArrays;
import org.apache.commons.math3.util.Pair;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.coveragemodel.linalg.FourierLinearOperator;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import javax.annotation.Nonnull;
import java.io.Serializable;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Extra MathUtils that should be moved to gatk-public
 * Created by davidben on 1/22/16.
 */
public class GATKProtectedMathUtils implements Serializable {

    private static final long serialVersionUID = -7587940242569731513L;

    private GATKProtectedMathUtils() {
    }

    public static final double INV_LOG_2 = 1.0 / Math.log(2.0);

    public static final double INV_LOG_10 = 1.0 / Math.log(10);

    /**
     * Threshold used to determine best way to calculate log(1-exp(a))
     * based on https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
     */
    private static final double LN_1_M_EXP_THRESHOLD = - Math.log(2);

    /**
     * Computes $\log(\sum_i e^{a_i})$ trying to avoid underflow issues by using the log-sum-exp trick.
     *
     * <p>
     * This trick consists on shift all the log values by the maximum so that exponent values are
     * much larger (close to 1) before they are summed. Then the result is shifted back down by
     * the same amount in order to obtain the correct value.
     * </p>
     * @return any double value.
     */
    public static double logSumExp(final double ... values) {
        double max = MathUtils.arrayMax(Utils.nonNull(values));
        double sum = 0.0;
        for (int i = 0; i < values.length; ++i) {
            if (values[i] != Double.NEGATIVE_INFINITY) {
                sum += Math.exp(values[i] - max);
            }
        }
        return max + Math.log(sum);
    }

    public static double logSumExp(final Collection<Double> values) {
        double max = Collections.max(values);
        double sum = 0.0;
        for (final double val : values) {
            if (val != Double.NEGATIVE_INFINITY) {
                sum += Math.exp(val - max);
            }
        }
        return max + Math.log(sum);
    }

    public static double mean(final double ... values) {
        Utils.nonNull(values);
        return MathUtils.mean(values, 0, values.length);
    }

    public static double[] rowMeans(final RealMatrix matrix) {
        Utils.nonNull(matrix);
        return IntStream.range(0, matrix.getRowDimension())
                .mapToDouble(r -> mean(matrix.getRow(r))).toArray();
    }

    public static double[] rowVariances(final RealMatrix matrix) {
        Utils.nonNull(matrix);
        final Variance varianceEvaluator = new Variance();
        return IntStream.range(0, matrix.getRowDimension())
                .mapToDouble(r -> varianceEvaluator.evaluate(matrix.getRow(r))).toArray();
    }

    /**
     * Calculates the standard deviation per row from a matrix.
     * @param matrix the input matrix.
     * @return never {@code null}, an array with as many positions as rows in {@code matrix}.
     * @throws IllegalArgumentException if {@code matrix} is {@code null}.
     */
    public static double[] rowStdDevs(final RealMatrix matrix) {
        Utils.nonNull(matrix);
        final Variance varianceEvaluator = new Variance();
        return IntStream.range(0, matrix.getRowDimension())
                .mapToDouble(r -> Math.sqrt(varianceEvaluator.evaluate(matrix.getRow(r)))).toArray();
    }

    /**
     * Calculates the mean per column from a matrix.
     * @param matrix the input matrix.
     * @return never {@code null}, an array with as many positions as columns in {@code matrix}.
     * @throws IllegalArgumentException if {@code matrix} is {@code null}.
     */
    public static double[] columnMeans(final RealMatrix matrix) {
        Utils.nonNull(matrix);
        return IntStream.range(0, matrix.getColumnDimension())
                .mapToDouble(c -> mean(matrix.getColumn(c))).toArray();
    }

    /**
     * Calculates the standard deviation per column from a matrix.
     * @param matrix the input matrix.
     * @return never {@code null}, an array with as many positions as columns in {@code matrix}.
     * @throws IllegalArgumentException if {@code matrix} is {@code null}.
     */
    public static double[] columnStdDevs(final RealMatrix matrix) {
        Utils.nonNull(matrix);
        final Variance varianceEvaluator = new Variance();
        return IntStream.range(0, matrix.getColumnDimension())
                .mapToDouble(c -> Math.sqrt(varianceEvaluator.evaluate(matrix.getColumn(c)))).toArray();
    }

    /**
     * Calculate the standard deviation of a collection of {@link Number} instances.
     * @param values the input values.
     * @return the standard deviation.
     * @throws IllegalArgumentException if {@code values} is {@code null} or it contains {@code null}.
     */
    public static double stdDev(final Collection<? extends Number> values) {
        Utils.nonNull(values);
        Utils.containsNoNull(values, "input values must not contain a null");
        final double[] doubleValues = values.stream()
                .mapToDouble(Number::doubleValue).toArray();
        return stdDev(doubleValues);
    }

    /**
     * Calculate the standard deviation of a double array.
     * @param values the input values.
     * @return the standard deviation.
     * @throws IllegalArgumentException if {@code values} is {@code null}.
     */
    public static double stdDev(final double ... values) {
        Utils.nonNull(values);
        return Math.sqrt(new Variance().evaluate(values));
    }

    /**
     * Create a set number of points, linearly spaced, between a minimum and maximum.
     *
     * Inspired by http://stackoverflow.com/questions/6208878/java-version-of-matlabs-linspace
     *
     * NaN are allowed, but will likely give useless results.
     *
     * @param min starting value
     * @param max ending value
     * @param points number of points, must be greater than -1
     * @return Never {@code null}
     */
    public static double[] createEvenlySpacedPoints(final double min, final double max, int points) {
        ParamUtils.isPositiveOrZero(points, "Number of points must be >= 0");
        if (points == 1) {
            return new double[] {max};
        }
        if (points == 0) {
            return new double[] {};
        }
        return IntStream.range(0, points).mapToDouble(n -> min + n * (max - min) / (points - 1)).toArray();
    }

    // given a list of options and a function for calculating their probabilities (these must sum to 1 over the whole list)
    // randomly choose one option from the implied categorical distribution
    public static <E> E randomSelect(final List<E> choices, final Function<E, Double> probabilityFunction, final RandomGenerator rng) {
        final List<Pair<E, Double>> pmf = choices.stream()
                .map(e -> new Pair<>(e, probabilityFunction.apply(e))).collect(Collectors.toList());
        return new EnumeratedDistribution<>(rng, pmf).sample();
    }

    /**
     *  Return an array with column sums in each entry
     *
     * @param matrix Never {@code null}
     * @return Never {@code null}
     */
    public static double[] columnSums(final RealMatrix matrix) {
        Utils.nonNull(matrix);

        return IntStream.range(0, matrix.getColumnDimension())
                .mapToDouble(c -> MathUtils.sum(matrix.getColumn(c))).toArray();
    }

    /**
     *  Return an array with row sums in each entry
     *
     * @param matrix Never {@code null}
     * @return Never {@code null}
     */
    public static double[] rowSums(final RealMatrix matrix) {
        Utils.nonNull(matrix);

        return IntStream.range(0, matrix.getRowDimension())
                .mapToDouble(r -> MathUtils.sum(matrix.getRow(r))).toArray();
    }

    /**
     *  Return sum of 3d array
     *
     * @param array Never {@code null}
     * @return sum of array
     */
    public static double sum(final double[][][] array) {
        Utils.nonNull(array);
        double result = 0;
        for (double[][] d: array) {
            for (double[] da: d){
                for (double daa: da) {
                    result += daa;
                }
            }

        }
        return result;
    }

    public static int minIndex(final int ... values) {
        Utils.nonNull(values);
        if (values.length == 0) {
            return -1;
        }
        int minValue = values[0];
        int minIndex = 0;
        for (int i = 0; i < values.length; i++) {
            final int nextValue = values[i];
            if (nextValue < minValue) {
                minValue = nextValue;
                minIndex = i;
            }
        }
        return minIndex;
    }

    /**
     * Given an int array returns the difference between the second smallest element
     * and the smallest element.
     * @param values the input values.
     * @param defaultValue value to return in case the input array has zero or only one element.
     * @return 0 if the input array has less than 2 elements, otherwise the second smallest
     * - the smallest.
     * @throws IllegalArgumentException if {@code values} is {@code null}.
     */
    public static int secondSmallestMinusSmallest(final int[] values, final int defaultValue) {
        Utils.nonNull(values);
        if (values.length <= 1) {
            return defaultValue;
        } else {
            int smallest = values[0];
            int secondSmallest = Integer.MAX_VALUE;
            for (int i = 1; i < values.length; i++) {
                if (values[i] < smallest) {
                    secondSmallest = smallest;
                    smallest = values[i];
                } else if (values[i] < secondSmallest) {
                    secondSmallest = values[i];
                }
            }
            return secondSmallest - smallest;
        }
    }

    /**
     * Returns the smallest power of 2 that exceeds or equals {@code val}
     * @param val input value
     * @return power of 2 integer
     * @throws IllegalArgumentException if an int overlow is encountered
     */
    public static int smallestPowerOfTwoGreaterThan(final int val)
            throws IllegalArgumentException {
        if (val > Integer.MAX_VALUE/2 || val < 0)
            throw new IllegalArgumentException("The mallest power of 2 greater than " + val + " is greater than Integer.MAX_VALUE" +
                    " or negative input.");
        return val > 1 ? Integer.highestOneBit(val - 1) << 1 : 1;
    }

    /**
     * Shrinks or expands a double array uniformly to a given length >= 2 ({@code newLength}) using nearest
     * neighbor interpolation. The routine works as follows: it puts the array on a uniform partition of the unit
     * interval on the real axis and and overlaps it with the desired grid. The value at a new grid point is
     * set to its nearest neighbor from the original grid. For example, if data = {1, 2, 3} and newLength = 4,
     * the output will be {1, 2, 2, 3}:
     *
     *      data     = 1     2     3
     *      old grid = x     x     x
     *      new grid = x   x   x   x
     *      new data = 1   2   2   3
     *
     * @param data original array
     * @param newLength length of the new array
     * @return interpolated array
     * @throws IllegalArgumentException if the input array is empty
     */
    public static double[] nearestNeighborUniform1DInterpolate(@Nonnull final double[] data, final int newLength)
            throws IllegalArgumentException {
        Utils.validateArg(data.length > 0, "The input array is empty.");
        Utils.validateArg(newLength >= 2, "The new length of the array must be >= 2");
        final double fact = (double)(data.length - 1)/(newLength - 1);
        return IntStream.range(0, newLength).mapToDouble(i -> data[(int) FastMath.floor(i*fact + 0.5)]).toArray();
    }

    /**
     * Find the maximum difference between entries of two arrays.  This is useful for testing convergence, for example
     */
    public static double maxDifference(final List<Double> array1, final List<Double> array2) {
        Utils.validateArg(array1.size() == array2.size(), "arrays must have same length.");
        Utils.validateArg(array1.size() > 0, "arrays must be non-empty");
        return IntStream.range(0, array1.size()).mapToDouble(n -> Math.abs(array1.get(n) - array2.get(n))).max().getAsDouble();
    }

    public static double[] posteriors(double[] log10Priors, double[] log10Likelihoods) {
        return MathUtils.normalizeFromLog10ToLinearSpace(MathArrays.ebeAdd(log10Priors, log10Likelihoods));
    }

    @FunctionalInterface
    public interface IntToDoubleArrayFunction {
        double[] apply(int value);
    }

    // sum of int -> double[] function mapped to an index range
    public static double[] sumArrayFunction(final int min, final int max, final IntToDoubleArrayFunction function) {
        Utils.validateArg(max >= min, "max must be at least as great as min");
        final double[] result = function.apply(min);
        for (int n = min + 1; n < max; n++) {
            final double[] newValues = function.apply(n);
            Utils.validateArg(newValues.length == result.length, "array function returns different sizes for different inputs!");
            for (int i = 0; i < result.length; i++) {
                result[i] += newValues[i];
            }
        }
        return result;
    }

    public static void addToArrayInPlace(final double[] array, final double[] summand) {
        Utils.validateArg(array.length == summand.length, "Arrays must have same length");
        for (int n = 0; n < array.length; n++) {
            array[n] += summand[n];
        }
    }

	public static int median(final int[] values) {
        return (int) FastMath.round(new Median().evaluate(Arrays.stream(values).mapToDouble(n -> n).toArray()));
    }

    public static double dotProduct(double[] a, double[] b){
        return MathUtils.sum(MathArrays.ebeMultiply(a, b));
    }

    /**
     * Calculates the complement of a log probability.
     *
     * <p>
     *     With complement of {@code x} we mean: {@code log(1-log(x))}.
     * </p>
     * @param x the input log probability.
     * @return {@code log(1-log(x))}
     */
    public static double logProbComplement(final double x) {
        return x >= LN_1_M_EXP_THRESHOLD
                ? Math.log(-Math.expm1(x))
                : Math.log1p(-Math.exp(x));
    }

    /**
     * Transform a log scaled probability (x) into the Phred scaled
     * equivalent or its complement (1-x) Phred scaled equivalent.
     * <p>
     *     This method tolerates probabilities slightly larger than 1.0
     *     (> 0.0 in log scale) which may occur occasionally due to
     *     float point calculation rounding.
     * </p>
     * <p>
     *     The value returned is a phred score capped by {@code maxPhredScore}.
     * </p>
     *
     * @param rawLogProb the probability.
     * @param complement whether to return the direct Phred transformation ({@code false})
     *                    or its complement ({@code true)}.
     * @param phredScorePrecision Phred score precision (i.e. quantization unit)
     * @param maxLogProb maximum tolerated log probability
     * @param maxPhredScore maximum reported Phred score
     * @return a values between 0 and {@code maxPhredScore}.
     * @throws GATKException if {@code rawLogProb} is larger than {@code maxLogProb}.
     */
    public static double logProbToPhredScore(final double rawLogProb, final boolean complement,
                                             final double maxPhredScore, final double maxLogProb,
                                             final double phredScorePrecision) {
        if (rawLogProb > maxLogProb) {
            throw new GATKException(String.format("possible numerical instability problem: the log-probability is too" +
                    " large: %g > 0.0 (with maximum tolerance %g)", rawLogProb, maxLogProb));
        }
        /* make sure that the probability is less than 1 in linear scale. there are cases where
         * log probability exceeds 0.0 due to floating point errors. */
        final double logProbEqOrLessThan0 = Math.min(0.0, rawLogProb);

        // Accurate way to calculate log(1-exp(a))
        // based on https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
        final double finalLogProb = complement
                ? logProbComplement(logProbEqOrLessThan0)
                : logProbEqOrLessThan0;

        final double absoluteQualScore = QualityUtils.phredScaleLog10ErrorRate(finalLogProb * INV_LOG_10);
        final double exactValue = Math.min(maxPhredScore, absoluteQualScore);
        // We round the value to the required precession.
        return roundPhred(exactValue, phredScorePrecision);
    }

    /**
     * Round a Phred scaled score to precision {@code phredScorePrecision}
     *
     * @param value raw Phred score
     * @param phredScorePrecision Phred score precision
     * @return rounded Phred score
     */
    public static double roundPhred(final double value, final double phredScorePrecision) {
        return Math.round(value / phredScorePrecision) * phredScorePrecision;
    }

    /**
     * Generate Fourier factors for a midpass filter; refer to {@link FourierLinearOperator}
     *
     * @param dimension dimension of the signal
     * @param freqLowerCutoff lower frequency cutoff (excluded)
     * @param freqUpperCutoff upper frequency cutoff (included)
     * @return Fourier factors
     */
    public static double[] getMidpassFilterFourierFactors(final int dimension, final int freqLowerCutoff,
                                                          final int freqUpperCutoff) {
        ParamUtils.isPositiveOrZero(freqUpperCutoff - freqLowerCutoff, "Upper cutoff must be >= lower cutoff.");
        ParamUtils.isPositive(dimension/2 + 1 - freqUpperCutoff, "For dimension " + dimension +
                ", the upper frequency cutoff must be <= " + (dimension/2));
        return IntStream.range(0, dimension/2 + 1)
                .mapToDouble(i -> i > freqLowerCutoff && i <= freqUpperCutoff ? 1.0 : 0).toArray();
    }
}
