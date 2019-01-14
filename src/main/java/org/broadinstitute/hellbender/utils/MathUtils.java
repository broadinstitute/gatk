package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.distribution.EnumeratedDistribution;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.exception.NumberIsTooLargeException;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.MathArrays;
import org.apache.commons.math3.util.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import javax.annotation.Nonnull;
import java.util.*;
import java.util.function.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * MathUtils is a static class (no instantiation allowed!) with some useful math methods.
 */
public final class MathUtils {

    /**
     * The smallest log10 value we'll emit from normalizeFromLog10 and other functions
     * where the real-space value is 0.0.
     */
    public static final double LOG10_P_OF_ZERO = -1000000.0;

    public static final double LOG10_ONE_HALF = Math.log10(0.5);

    public static final double LOG_ONE_HALF = FastMath.log(0.5);

    public static final double LOG10_ONE_THIRD = -Math.log10(3.0);
    public static final double INV_LOG_2 = 1.0 / Math.log(2.0);
    public static final double INV_LOG_10 = 1.0 / Math.log(10);

    private static final double LN_10 = Math.log(10);

    private static final double LOG1MEXP_THRESHOLD = Math.log(0.5);

    public static final double INV_SQRT_2_PI = 1.0 / Math.sqrt(2.0 * Math.PI);

    private static final double NATURAL_LOG_OF_TEN = Math.log(10.0);

    private static final double SQUARE_ROOT_OF_TWO_TIMES_PI = Math.sqrt(2.0 * Math.PI);

    /**
     * Log10 of the e constant.
     */
    public static final double LOG10_OF_E = Math.log10(Math.E);
    public static final double FAIR_BINOMIAL_PROB_LOG10_0_5 = Math.log10(0.5);
    /**
     * Threshold used to determine best way to calculate log(1-exp(a))
     * based on https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
     */
    private static final double LN_1_M_EXP_THRESHOLD = - Math.log(2);

    /**
     * Private constructor.  No instantiating this class!
     */
    private MathUtils() { }

    /**
     * Computes $\log(\sum_i e^{a_i})$ trying to avoid underflow issues by using the log-sum-exp trick.
     *
     * <p>
     * This trick consists of shifting all the log values by the maximum so that exponent values are
     * much larger (close to 1) before they are summed. Then the result is shifted back down by
     * the same amount in order to obtain the correct value.
     * </p>
     * @return any double value.
     */
    public static double logSumExp(final double ... values) {
        double max = arrayMax(Utils.nonNull(values));
        double sum = 0.0;
        for (int i = 0; i < values.length; ++i) {
            if (values[i] != Double.NEGATIVE_INFINITY) {
                sum += Math.exp(values[i] - max);
            }
        }
        return max + Math.log(sum);
    }

    public static double logSumExp(final Collection<Double> values) {
        double max = Collections.max(Utils.nonNull(values));
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
        return mean(values, 0, values.length);
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
        Utils.nonNull(choices);
        Utils.nonNull(probabilityFunction);
        Utils.nonNull(rng);
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
                .mapToDouble(c -> sum(matrix.getColumn(c))).toArray();
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
                .mapToDouble(r -> sum(matrix.getRow(r))).toArray();
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
     * @throws IllegalArgumentException if an int overflow is encountered
     */
    public static int smallestPowerOfTwoGreaterThan(final int val) {
        ParamUtils.inRange(val, 0, Integer.MAX_VALUE/2, "The smallest power of 2 greater than this value is greater than Integer.MAX_VALUE or negative input.");
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
    public static double[] nearestNeighborUniform1DInterpolate(@Nonnull final double[] data, final int newLength) {
        Utils.nonNull(data);
        Utils.validateArg(data.length > 0, "The input array is empty.");
        Utils.validateArg(newLength >= 2, "The new length of the array must be >= 2");
        final double fact = (double)(data.length - 1)/(newLength - 1);
        return IntStream.range(0, newLength).mapToDouble(i -> data[(int) FastMath.floor(i*fact + 0.5)]).toArray();
    }

    /**
     * Find the maximum difference between entries of two arrays.  This is useful for testing convergence, for example
     */
    public static double maxDifference(final List<Double> array1, final List<Double> array2) {
        Utils.nonNull(array1);
        Utils.nonNull(array2);
        Utils.validateArg(array1.size() == array2.size(), "arrays must have same length.");
        Utils.validateArg(array1.size() > 0, "arrays must be non-empty");
        return IntStream.range(0, array1.size()).mapToDouble(n -> Math.abs(array1.get(n) - array2.get(n))).max().getAsDouble();
    }

    public static double[] posteriors(double[] log10Priors, double[] log10Likelihoods) {
        return normalizeFromLog10ToLinearSpace(MathArrays.ebeAdd(log10Priors, log10Likelihoods));
    }

    public static int[] normalizePLs(int[] PLs) {
        final int[] newPLs = new int[PLs.length];
        final int smallest = arrayMin(PLs);
        for(int i=0; i<PLs.length; i++) {
            newPLs[i] = PLs[i]-smallest;
        }
        return newPLs;
    }

    public static int[] ebeAdd(final int[] a, final int[] b) throws DimensionMismatchException {
            if (a.length != b.length) {
                throw new DimensionMismatchException(a.length, b.length);
            }

            final int[] result = a.clone();
            for (int i = 0; i < a.length; i++) {
                result[i] += b[i];
            }
            return result;
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
        Utils.nonNull(values);
        return (int) FastMath.round(new Median().evaluate(Arrays.stream(values).mapToDouble(n -> n).toArray()));
    }

    public static double dotProduct(double[] a, double[] b){
        return sum(MathArrays.ebeMultiply(Utils.nonNull(a), Utils.nonNull(b)));
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

    @FunctionalInterface
    public interface IntToDoubleArrayFunction {
        double[] apply(int value);
    }


    /**
     * A utility class that computes on the fly average and standard deviation for a stream of numbers.
     * The number of observations does not have to be known in advance, and can be also very big (so that
     * it could overflow any naive summation-based scheme or cause loss of precision).
     * Instead, adding a new number <code>observed</code>
     * to a sample with <code>add(observed)</code> immediately updates the instance of this object so that
     * it contains correct mean and standard deviation for all the numbers seen so far. Source: Knuth, vol.2
     * (see also e.g. http://www.johndcook.com/standard_deviation.html for online reference).
     */
    public static class RunningAverage {
        private double mean = 0.0;
        private double s = 0.0;
        private long obs_count = 0;

        public void add(double obs) {
            obs_count++;
            double oldMean = mean;
            mean += (obs - mean) / obs_count; // update mean
            s += (obs - oldMean) * (obs - mean);
        }

        public void addAll(Collection<Number> col) {
            for (Number o : col) {
                add(o.doubleValue());
            }
        }

        public double mean() {
            return mean;
        }

        public double stddev() {
            return Math.sqrt(s / (obs_count - 1));
        }

        public double var() {
            return s / (obs_count - 1);
        }

        public long observationCount() {
            return obs_count;
        }

        @Override
        public RunningAverage clone() {
            RunningAverage ra = new RunningAverage();
            ra.mean = this.mean;
            ra.s = this.s;
            ra.obs_count = this.obs_count;
            return ra;
        }

        public void merge(RunningAverage other) {
            if (this.obs_count > 0 || other.obs_count > 0) { // if we have any observations at all
                this.mean = (this.mean * this.obs_count + other.mean * other.obs_count) / (this.obs_count + other.obs_count);
                this.s += other.s;
            }
            this.obs_count += other.obs_count;
        }
    }


    /**
     * Get a random int between min and max (inclusive) using the global GATK random number generator.
     * It is required that max > min - 1.
     *
     * @param min lower bound of the range
     * @param max upper bound of the range
     * @return a random int >= min and <= max
     */
    public static int randomIntegerInRange( final int min, final int max ) {
        Utils.validateArg(max > min - 1, () -> "invalid arguments min:" + min + " max:" + max);
        return Utils.getRandomGenerator().nextInt(max - min + 1) + min;
    }

    /**
     * Computes the root mean square of the given array or 0.0 if the array is empty.
     */
    public static double rms(final Collection<Integer> list) {
        Utils.nonNull(list);
        return Math.sqrt(list.stream().mapToInt(i -> i * i).average().orElse(0.0));
    }

    /**
     * Computes the sum of squares of the given collection of ints or 0.0 if the collection is empty.
     */
    public static double sumOfSquares(final Collection<Integer> collection) {
        Utils.nonNull(collection);
        return collection.stream().mapToInt(i -> i * i).sum();
    }

    /**
     * Creates a new sample of k ints from [0..n-1], without duplicates.
     * @throws NumberIsTooLargeException if {@code k > n}.
     * @throws NotStrictlyPositiveException if {@code k <= 0}.
     */
    public static int[] sampleIndicesWithoutReplacement(final int n, final int k) {
        //No error checking : RandomDataGenetator.nextPermutation does it
        return Utils.getRandomDataGenerator().nextPermutation(n, k);
    }

    /**
     * Calculates {@code log10(1-10^a)} without losing precision.
     *
     * <p>
     *     This is based on the approach described in:
     *
     * </p>
     * <p>
     *     Maechler M, Accurately Computing log(1-exp(-|a|)) Assessed by the Rmpfr package, 2012 <br/>
     *     <a ref="http://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf">Online document</a>.
     * </p>
     *
     * @param a the input exponent.
     * @return {@link Double#NaN NaN} if {@code a > 0}, otherwise the corresponding value.
     */
    public static double log10OneMinusPow10(final double a) {
        if (a > 0) return Double.NaN;
        if (a == 0) return Double.NEGATIVE_INFINITY;
        final double b = a * LN_10;
        return log1mexp(b) / LN_10;
    }

    /**
     * Calculates {@code log(1-exp(a))} without losing precision.
     *
     * <p>
     *     This is based on the approach described in:
     *
     * </p>
     * <p>
     *     Maechler M, Accurately Computing log(1-exp(-|a|)) Assessed by the Rmpfr package, 2012 <br/>
     *     <a ref="http://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf">Online document</a>.
     *
     * </p>
     *
     * @param a the input exponent.
     * @return {@link Double#NaN NaN} if {@code a > 0}, otherwise the corresponding value.
     */
    public static double log1mexp(final double a) {
        if (a > 0) return Double.NaN;
        if (a == 0) return Double.NEGATIVE_INFINITY;
        return (a < LOG1MEXP_THRESHOLD) ? Math.log1p(-Math.exp(a)) : Math.log(-Math.expm1(a));
    }

    /**
     * Check that the log10 prob vector vector is well formed
     *
     * @param vector
     * @param expectedSize
     * @param shouldSumToOne
     *
     * @return true if vector is well-formed, false otherwise
     */
    public static boolean goodLog10ProbVector(final double[] vector, final int expectedSize, final boolean shouldSumToOne) {
        Utils.nonNull(vector);
        return vector.length == expectedSize &&
                allMatch(vector, MathUtils::goodLog10Probability) &&
                !( shouldSumToOne && compareDoubles(sumLog10(vector), 1.0, 1e-4) != 0 );
    }

    /**
     *  Returns the sum of values whose log10s we have. That is, returns sum(10^x_i).
     */
    public static double sumLog10(final double[] log10values) {
        return Math.pow(10.0, log10SumLog10(Utils.nonNull(log10values)));
    }

    /** Compute Z=X-Y for two numeric vectors X and Y
     *
     * @param x                 First vector
     * @param y                 Second vector
     * @return Vector of same length as x and y so that z[k] = x[k]-y[k]
     */
    public static int[] vectorDiff(final int[] x, final int[] y) {
        Utils.nonNull(x, "x is null");
        Utils.nonNull(y, "y is null");
        Utils.validateArg(x.length == y.length, "Lengths of x and y must be the same");
        return new IndexRange(0, x.length).mapToInteger(k -> x[k] - y[k]);
    }


    /**
     * Calculates the log10 of the multinomial coefficient. Designed to prevent
     * overflows even with very large numbers.
     *
     * @param n total number of trials
     * @param k array of any size with the number of successes for each grouping (k1, k2, k3, ..., km)
     * @return {@link Double#NaN NaN} if {@code a > 0}, otherwise the corresponding value.
     */
    public static double log10MultinomialCoefficient(final int n, final int[] k) {
        Utils.validateArg(n >= 0, "n: Must have non-negative number of trials");
        Utils.validateArg(allMatch(k, x -> x >= 0), "Elements of k must be non-negative");
        Utils.validateArg(sum(k) == n, "Sum of observations k must sum to total number of trials n");
        return log10Factorial(n) -  new IndexRange(0, k.length).sum(j -> log10Factorial(k[j]));
    }

    public static double log10(int i) {
        return Log10Cache.get(i);
    }

    public static double log10sumLog10(final double[] log10values) {
        return log10sumLog10(Utils.nonNull(log10values), 0);
    }

    public static double log10sumLog10(final double[] log10p, final int start) {
        return log10sumLog10(Utils.nonNull(log10p), start, log10p.length);
    }

    public static double log10sumLog10(final double[] log10p, final int start, final int finish) {
        Utils.nonNull(log10p);
        if (start >= finish) {
            return Double.NEGATIVE_INFINITY;
        }
        final int maxElementIndex = maxElementIndex(log10p, start, finish);
        final double maxValue = log10p[maxElementIndex];
        if(maxValue == Double.NEGATIVE_INFINITY) {
            return maxValue;
        }
        final double sum = 1.0 + new IndexRange(start, finish).sum(i -> i == maxElementIndex ? 0 : Math.pow(10.0, log10p[i] - maxValue));
        Utils.validateArg(!Double.isNaN(sum) && sum != Double.POSITIVE_INFINITY, "log10p values must be non-infinite and non-NAN");
        return maxValue + Math.log10(sum);
    }


    /**
     * A helper class to maintain a cache of log10 values.
     * The cache expands when a number is not available.
     * NOTE: this cache is thread safe and it may be accessed from multiple threads.
     */
    private static final class Log10Cache {
        private static final Logger logger = LogManager.getLogger(Log10Cache.class);

        //initialize with the special case: log(0) = NEGATIVE_INFINITY
        private static double[] cache = new double[] { Double.NEGATIVE_INFINITY };

        /**
         * Get the value of log10(n), expanding the cache as necessary
         * @param i operand
         * @return log10(n)
         */
        public static double get(final int i) {
            Utils.validateArg(i >= 0, () -> String.format("Can't take the log of a negative number: %d", i));
            if (i >= cache.length) {
                final int newCapacity = Math.max(i + 10, 2 * cache.length);
                logger.debug("cache miss " + i + " > " + (cache.length-1) + " expanding to " + newCapacity);
                expandCache(newCapacity);
            }
            /*
               Array lookups are not atomic.  It's possible that the reference to cache could be
               changed between the time the reference is loaded and the data is fetched from the correct
               offset.  However, the value retrieved can't change, and it's guaranteed to be present in the
               old reference by the conditional above.
             */
            return cache[i];
        }

        /**
         * Ensures that the cache contains a value for n.  After completion of expandCache(n),
         * #get(n) is guaranteed to return without causing a cache expansion
         * @param newCapacity desired value to be precomputed
         */
        private static synchronized void expandCache(final int newCapacity) {
            if (newCapacity < cache.length) {
                //prevents a race condition when multiple threads want to expand the cache at the same time.
                //in that case, one of them will be first to enter the synchronized method expandCache and
                //so the others may end up in this method even if n < cache.length
                return;
            }
            final double[] newCache = new double[newCapacity + 1];
            System.arraycopy(cache, 0, newCache, 0, cache.length);
            for (int i = cache.length; i < newCache.length; i++) {
                newCache[i] = Math.log10(i);
            }
            cache = newCache;
        }

        public static int size() {
            return cache.length;
        }

    }

    /**
     * Encapsulates the second term of Jacobian log identity for differences up to MAX_TOLERANCE
     */
    private static final class JacobianLogTable {
        // if log(a) - log(b) > MAX_TOLERANCE, b is effectively treated as zero in approximateLogSumLog
        // MAX_TOLERANCE = 8.0 introduces an error of at most one part in 10^8 in sums
        public static final double MAX_TOLERANCE = 8.0;

        //  Phred scores Q and Q+1 differ by 0.1 in their corresponding log-10 probabilities, and by
        // 0.1 * log(10) in natural log probabilities.  Setting TABLE_STEP to an exact divisor of this
        // quantity ensures that approximateSumLog in fact caches exact values for integer phred scores
        private static final double TABLE_STEP = 0.0001;
        private static final double INV_STEP = 1.0 / TABLE_STEP;
        private static final double[] cache = new IndexRange(0, (int) (MAX_TOLERANCE / TABLE_STEP) + 1)
                .mapToDouble(k -> Math.log10(1.0 + Math.pow(10.0, -k * TABLE_STEP)));

        public static double get(final double difference) {
            final int index = fastRound(difference * INV_STEP);
            return cache[index];
        }
    }

    // A fast implementation of the Math.round() method.  This method does not perform
    // under/overflow checking, so this shouldn't be used in the general case (but is fine
    // if one is already make those checks before calling in to the rounding).
    public static int fastRound(final double d) {
        return (d > 0.0) ? (int) (d + 0.5d) : (int) (d - 0.5d);
    }

    /**
     * Converts LN to LOG10
     * @param ln log(x)
     * @return log10(x)
     */
    public static double logToLog10(final double ln) {
        return ln * LOG10_OF_E;
    }

    public static double approximateLog10SumLog10(final double[] vals) {
        return approximateLog10SumLog10(Utils.nonNull(vals), vals.length);
    }

    public static double approximateLog10SumLog10(final double[] vals, final int endIndex) {
        Utils.nonNull(vals);
        final int maxElementIndex = MathUtils.maxElementIndex(vals, endIndex);
        double approxSum = vals[maxElementIndex];

        for (int i = 0; i < endIndex; i++) {
            if (i == maxElementIndex || vals[i] == Double.NEGATIVE_INFINITY) {
                continue;
            }

            // if vals[i] isn't too tiny relative to the sum so far, add it; otherwise ignore it
            final double diff = approxSum - vals[i];
            approxSum += diff < JacobianLogTable.MAX_TOLERANCE ? JacobianLogTable.get(diff) : 0.0;
        }

        return approxSum;
    }

    public static double approximateLog10SumLog10(final double a, final double b, final double c) {
        return approximateLog10SumLog10(a, approximateLog10SumLog10(b, c));
    }

    public static double approximateLog10SumLog10(final double a, final double b) {
        // this code works only when a <= b so we flip them if the order is opposite
        if (a > b) {
            return approximateLog10SumLog10(b, a);
        } else if (a == Double.NEGATIVE_INFINITY) {
            return b;
        }

        // if |b-a| < tol we need to compute log(e^a + e^b) = log(e^b(1 + e^(a-b))) = b + log(1 + e^(-(b-a)))
        // we compute the second term as a table lookup with integer quantization
        // we have pre-stored correction for 0,0.1,0.2,... 10.0
        final double diff = b - a;
        return b + (diff < JacobianLogTable.MAX_TOLERANCE ? JacobianLogTable.get(diff) : 0.0);
    }

    /**
     * Calculate the approximate log10 sum of an array range.
     * @param vals the input values.
     * @param fromIndex the first inclusive index in the input array.
     * @param toIndex index following the last element to sum in the input array (exclusive).
     * @return the approximate sum.
     * @throws IllegalArgumentException if {@code vals} is {@code null} or  {@code fromIndex} is out of bounds
     * or if {@code toIndex} is larger than
     * the length of the input array or {@code fromIndex} is larger than {@code toIndex}.
     */
    public static double approximateLog10SumLog10(final double[] vals, final int fromIndex, final int toIndex) {
        Utils.nonNull(vals);
        if (fromIndex == toIndex) return Double.NEGATIVE_INFINITY;
        final int maxElementIndex = MathUtils.maxElementIndex(vals,fromIndex,toIndex);
        double approxSum = vals[maxElementIndex];

        for (int i = fromIndex; i < toIndex; i++) {
            final double val;
            if (i == maxElementIndex || (val = vals[i]) == Double.NEGATIVE_INFINITY)
                continue;
            final double diff = approxSum - val;
            if (diff < JacobianLogTable.MAX_TOLERANCE)
                approxSum += JacobianLogTable.get(diff);
        }
        return approxSum;
    }

    public static double sum(final double[] values) {
        Utils.nonNull(values);
        double s = 0.0;
        for (double v : values)
            s += v;
        return s;
    }

    public static long sum(final int[] x) {
        Utils.nonNull(x);
        long total = 0;
        for (int v : x)
            total += v;
        return total;
    }

    public static int sum(final byte[] x) {
        Utils.nonNull(x);
        int total = 0;
        for (byte v : x)
            total += (int)v;
        return total;
    }

    public static long sum(final long[] x) {
        Utils.nonNull(x);
        int total = 0;
        for (long v : x)
            total += v;
        return total;
    }

    /** Returns the sum of the elements in the array starting with start and ending before stop. */
    public static long sum(final long[] arr, final int start, final int stop) {
        long result = 0;
        for (int n = start; n < stop; n++) {
            result += arr[n];
        }
        return result;
    }

    /** Returns the sum of the elements in the array starting with start and ending before stop. */
    public static double sum(final double[] arr, final int start, final int stop) {
        Utils.nonNull(arr);
        Utils.validateArg(start <= stop, () -> start + " > " + stop);
        Utils.validateArg(start >= 0, () -> start + " < " + 0);
        Utils.validateArg(stop <= arr.length, () -> stop + " >  " + arr.length);
        double result = 0.0;
        for (int n = start; n < stop; n++) {
            result += arr[n];
        }
        return result;
    }

    public static <E> double sumDoubleFunction(final Collection<E> collection, final ToDoubleFunction<E> function) {
        double result = 0;
        for (final E e: collection) {
            result += function.applyAsDouble(e);
        }
        return result;
    }

    public static <E> int sumIntFunction(final Collection<E> collection, final ToIntFunction<E> function) {
        int result = 0;
        for (final E e: collection) {
            result += function.applyAsInt(e);
        }
        return result;
    }

    /**
     * Compares double values for equality (within 1e-6), or inequality.
     *
     * @param a the first double value
     * @param b the second double value
     * @return -1 if a is greater than b, 0 if a is equal to be within 1e-6, 1 if b is greater than a.
     */
    public static byte compareDoubles(final double a, final double b) {
        return compareDoubles(a, b, 1e-6);
    }

    /**
     * Compares double values for equality (within epsilon), or inequality.
     *
     * @param a       the first double value
     * @param b       the second double value
     * @param epsilon the precision within which two double values will be considered equal
     * @return -1 if a is greater than b, 0 if a is equal to be within epsilon, 1 if b is greater than a.
     */
    public static byte compareDoubles(final double a, final double b, final double epsilon) {
        if (Math.abs(a - b) < epsilon) {
            return 0;
        }
        if (a > b) {
            return -1;
        }
        return 1;
    }

    /**
     * Checks that value is between min and max, inclusive, with a specified tolerance.
     *
     * Does not check for NaNs or infinities.
     *
     * @param value value to check
     * @param min start of allowable range
     * @param max end of allowable range
     * @param tolerance perform double comparisons to within this tolerance
     * @return true if value is within the specified range (+/- tolerance), otherwise false
     */
    public static boolean doubleWithinRangeWithTolerance( final double value, final double min, final double max, final double tolerance ) {
        return value >= min - tolerance && value <= max + tolerance;
    }

    /**
     * Calculates the binomial coefficient. Designed to prevent
     * overflows even with very large numbers.
     *
     * @param n total number of trials
     * @param k number of successes
     * @return the binomial coefficient
     */
    public static double binomialCoefficient(final int n, final int k) {
        return Math.pow(10, log10BinomialCoefficient(n, k));
    }

    /**
     * @see #binomialCoefficient(int, int) with log10 applied to result
     */
    public static double log10BinomialCoefficient(final int n, final int k) {
        Utils.validateArg(n >= 0, "Must have non-negative number of trials");
        Utils.validateArg( k <= n && k >= 0, "k: Must have non-negative number of successes, and no more successes than number of trials");
        return log10Factorial(n) - log10Factorial(k) - log10Factorial(n - k);
    }

    /**
     * Computes a binomial probability.  This is computed using the formula
     * <p/>
     * B(k; n; p) = [ n! / ( k! (n - k)! ) ] (p^k)( (1-p)^k )
     * <p/>
     * where n is the number of trials, k is the number of successes, and p is the probability of success
     *
     * @param n number of Bernoulli trials
     * @param k number of successes
     * @param p probability of success
     * @return the binomial probability of the specified configuration.  Computes values down to about 1e-237.
     */
    public static double binomialProbability(final int n, final int k, final double p) {
        return Math.pow(10.0, log10BinomialProbability(n, k, Math.log10(p)));
    }

    /**
     * binomial Probability(int, int, double) with log applied to result
     */
    public static double log10BinomialProbability(final int n, final int k, final double log10p) {
        Utils.validateArg(log10p < 1.0e-18, "log10p: Log10-probability must be 0 or less");
        if (log10p == Double.NEGATIVE_INFINITY){
            return k == 0 ? 0 : Double.NEGATIVE_INFINITY;
        }
        double log10OneMinusP = Math.log10(1 - Math.pow(10.0, log10p));
        return log10BinomialCoefficient(n, k) + log10p * k + log10OneMinusP * (n - k);
    }

    /**
     * @see #binomialProbability(int, int, double) with p=0.5 and log10 applied to result
     */
    public static double log10BinomialProbability(final int n, final int k) {
        return log10BinomialCoefficient(n, k) + (n * FAIR_BINOMIAL_PROB_LOG10_0_5);
    }

    public static double log10SumLog10(final double[] log10Values, final int start) {
        return log10SumLog10(Utils.nonNull(log10Values), start, log10Values.length);
    }

    public static double log10SumLog10(final double[] log10Values) {
        return log10SumLog10(Utils.nonNull(log10Values), 0);
    }

    public static double log10SumLog10(final double[] log10Values, final int start, final int finish) {
        Utils.nonNull(log10Values);
        if (start >= finish) {
            return Double.NEGATIVE_INFINITY;
        }
        final int maxElementIndex = maxElementIndex(log10Values, start, finish);
        final double maxValue = log10Values[maxElementIndex];
        if(maxValue == Double.NEGATIVE_INFINITY) {
            return maxValue;
        }
        double sum = 1.0;
        for (int i = start; i < finish; i++) {
            final double curVal = log10Values[i];
            if (i == maxElementIndex || curVal == Double.NEGATIVE_INFINITY) {
                continue;
            } else {
                final double scaled_val = curVal - maxValue;
                sum += Math.pow(10.0, scaled_val);
            }
        }
        if ( Double.isNaN(sum) || sum == Double.POSITIVE_INFINITY ) {
            throw new IllegalArgumentException("log10 p: Values must be non-infinite and non-NAN");
        }
        return maxValue + (sum != 1.0 ? Math.log10(sum) : 0.0);
    }

    public static double log10SumLog10(final double a, final double b) {
        return a > b ? a + Math.log10(1 + Math.pow(10.0, b - a)) : b + Math.log10(1 + Math.pow(10.0, a - b));
    }

    public static double logSumLog(final double a, final double b) {
        return a > b ? a + FastMath.log(1 + FastMath.exp(b - a)) : b + FastMath.log(1 + FastMath.exp(a - b));
    }

    /**
     * Calculate f(x) = log10 ( Normal(x | mu = mean, sigma = sd) )
     * @param mean the desired mean of the Normal distribution
     * @param sd the desired standard deviation of the Normal distribution
     * @param x the value to evaluate
     * @return a well-formed double
     */
    public static double normalDistributionLog10(final double mean, final double sd, final double x) {
        Utils.validateArg( sd >= 0, "sd: Standard deviation of normal must be > 0");
        if ( ! wellFormedDouble(mean) || ! wellFormedDouble(sd) || ! wellFormedDouble(x) )
            throw new IllegalArgumentException("mean, sd, or, x : Normal parameters must be well formatted (non-INF, non-NAN)");
        final double a = -1.0 * Math.log10(sd * SQUARE_ROOT_OF_TWO_TIMES_PI);
        final double b = -1.0 * (square(x - mean) / (2.0 * square(sd))) / NATURAL_LOG_OF_TEN;
        return a + b;
    }

    /**
     * Calculate f(x) = x^2
     * @param x the value to square
     * @return x * x
     */
    public static double square(final double x) {
        return x * x;
    }

    public static double distanceSquared(final double[] x, final double[] y) {
        Utils.nonNull(x);
        Utils.nonNull(y);
        return new IndexRange(0, x.length).sum(n -> square(x[n] - y[n]));
    }

    /**
     * normalizes the log-probability array in-place.
     *
     * @param array the array to be normalized
     * @return the normalized-in-place array
     */
    public static double[] normalizeFromLog10ToLinearSpace(final double[] array) {
        return normalizeLog10(Utils.nonNull(array), false, true);
    }

    /**
     * normalizes the log-probability array in-place.
     *
     * @param array             the array to be normalized
     * @return the normalized-in-place array, maybe log transformed
     */
    public static double[] normalizeLog10(final double[] array) {
        return normalizeLog10(Utils.nonNull(array), true, true);
    }


    /**
     * See #normalizeFromLog but with the additional option to use an approximation that keeps the calculation always in log-space
     *
     * @param array
     * @param takeLog10OfOutput
     * @param inPlace           if true, modify the input array in-place
     *
     * @return
     */
    public static double[] normalizeLog10(final double[] array, final boolean takeLog10OfOutput, final boolean inPlace) {
        final double log10Sum = log10SumLog10(Utils.nonNull(array));
        final double[] result = inPlace ? applyToArrayInPlace(array, x -> x - log10Sum) : applyToArray(array, x -> x - log10Sum);
        return takeLog10OfOutput ? result : applyToArrayInPlace(result, x -> Math.pow(10.0, x));
    }


    //TODO: delete after we are satisfied with the concordance of VQSR with GATK3
    public static double[] normalizeLog10DeleteMePlease(final double[] array, final boolean takeLog10OfOutput) {
        Utils.nonNull(array);
        final double maxValue = arrayMax(array);
        final double[] normalized = applyToArray(array, x -> Math.pow(10.0, x - maxValue));
        final double sum = sum(normalized);
        if (!takeLog10OfOutput) {
            return applyToArrayInPlace(normalized, x -> x/sum);
        } else {
            final double log10Sum = Math.log10(sum);
            return applyToArrayInPlace(array, x -> x - maxValue - log10Sum);
        }
    }

    /**
     * Given an array of log space (log or log10) values, subtract all values by the array maximum so that the max element in log space
     * is zero.  This is equivalent to dividing by the maximum element in real space and is useful for avoiding underflow/overflow
     * when the array's values matter only up to an arbitrary normalizing factor, for example, an array of likelihoods.
     *
     * @param array
     * @return the scaled-in-place array
     */
    public static double[] scaleLogSpaceArrayForNumericalStability(final double[] array) {
        Utils.nonNull(array);
        final double maxValue = arrayMax(array);
        return applyToArrayInPlace(array, x -> x - maxValue);
    }

    /**
     * normalizes the real-space probability array.
     *
     * Does not assume anything about the values in the array, beyond that no elements are below 0.  It's ok
     * to have values in the array of > 1, or have the sum go above 0.
     *
     * @param array the array to be normalized
     * @return a newly allocated array corresponding the normalized values in array
     */
    public static double[] normalizeFromRealSpace(final double[] array) {
        Utils.nonNull(array);
        if ( array.length == 0 )
            return array;

        final double sum = sum(array);
        Utils.validateArg(sum >= 0.0, () -> "Values in probability array sum to a negative number " + sum);
        return applyToArray(array, x -> x/sum);
    }

    public static int maxElementIndex(final double[] array) {
        return maxElementIndex(Utils.nonNull(array), array.length);
    }

    public static int maxElementIndex(final double[] array, final int start, final int endIndex) {
        Utils.nonNull(array);
        Utils.validateArg(array.length > 0, "array may not be empty");
        Utils.validateArg(start <= endIndex, "Start cannot be after end.");

        int maxI = start;
        for (int i = (start+1); i < endIndex; i++) {
            if (array[i] > array[maxI])
                maxI = i;
        }
        return maxI;
    }

    public static int maxElementIndex(final double[] array, final int endIndex) {
        return maxElementIndex(Utils.nonNull(array), 0, endIndex);
    }

    public static double arrayMax(final double[] array) {
        Utils.nonNull(array);
        return array[maxElementIndex(array)];
    }

    public static int arrayMin(final int[] array) {
        Utils.nonNull(array);
        int min=array[0];
        for(int i=0; i<array.length; i++) {
            if(array[i] < min) {
                min = array[i];
            }
        }
        return min;
    }

    public static double arrayMin(final double[] array) {
        Utils.nonNull(array);
        double min=array[0];
        for(int i=0; i<array.length; i++) {
            if(array[i] < min) {
                min = array[i];
            }
        }
        return min;
    }

    /**
     * Checks that the result is a well-formed log10 probability
     *
     * @param result a supposedly well-formed log10 probability value
     * @return true if result is really well formed
     */
    public static boolean goodLog10Probability(final double result) {
        return result <= 0.0;
    }

    /**
     * Checks that the result is a well-formed probability
     *
     * @param result a supposedly well-formed probability value
     * @return true if result is really well formed
     */
    public static boolean goodProbability(final double result) {
        return result >= 0.0 && result <= 1.0;
    }

    public static double log10ToLog(final double log10){
        return log10 * LN_10;
    }

    /**
      * Calculates the log10 of the gamma function for x.
      *
      * @param x the x parameter
      * @return the log10 of the gamma function at x.
      */
    public static double log10Gamma(final double x) {
       return logToLog10(Gamma.logGamma(x));
    }

    public static double log10Factorial(final int n) {
        return n >=0 && n < Log10FactorialCache.size() ? Log10FactorialCache.get(n) : log10Gamma(n+1);
    }

    /**
     * Converts a real space array of numbers (typically probabilities) into a log10 array
     *
     * @param prRealSpace
     * @return
     */
    public static double[] toLog10(final double[] prRealSpace) {
        Utils.nonNull(prRealSpace);
        return applyToArray(prRealSpace, Math::log10);
    }

    /**
     * Wrapper class so that the log10Factorial array is only calculated if it's used
     */
    private static class Log10FactorialCache {

        /**
         * The size of the precomputed cache.  Must be a positive number!
         */
        private static final int CACHE_SIZE = 10_000;

        private static double[] cache = null;

        public static int size() { return CACHE_SIZE; }

        public static double get(final int n) {
            if (cache == null) {
                initialize();
            }
            return cache[n];
        }

        private static synchronized void initialize() {
            if (cache == null) {//this null check is here to prevent a race condition
                // when multiple threads want to initialize the cache
                Log10Cache.expandCache(CACHE_SIZE);
                cache = new double[CACHE_SIZE];
                cache[0] = 0.0;
                for (int k = 1; k < cache.length; k++) {
                    cache[k] = cache[k - 1] + Log10Cache.get(k);
                }
            }
        }
    }

    /**
     * Compute in a numerically correct way the quantity log10(1-x)
     *
     * Uses the approximation log10(1-x) = log10(1/x - 1) + log10(x) to avoid very quick underflow
     * in 1-x when x is very small
     *
     * @param x a positive double value between 0.0 and 1.0
     * @return an estimate of log10(1-x)
     */
    public static double log10OneMinusX(final double x) {
        if ( x == 1.0 )
            return Double.NEGATIVE_INFINITY;
        else if ( x == 0.0 )
            return 0.0;
        else {
            final double d = Math.log10(1 / x - 1) + Math.log10(x);
            return Double.isInfinite(d) || d > 0.0 ? 0.0 : d;
        }
    }

    /** Calculate the mean of an array of doubles. */
    public static double mean(final double[] in, final int start, final int stop) {
        Utils.nonNull(in);
        return stop <= start ? Double.NaN : Arrays.stream(in, start, stop).average().getAsDouble();
    }

    /** "Promotes" an int[] into a double array with the same values (or as close as precision allows). */
    public static double[] promote(final int[] is) {
        Utils.nonNull(is);
        return new IndexRange(0, is.length).mapToDouble(n -> (double) is[n]);
    }


    /**
     * Compute the median of a list of numbers
     *
     * If values.length is even, this will be the middle value when the elements are sorted
     * If values.length is odd then it will be the mean of the two values closest to the middle.
     *
     * @param values a list of numbers
     * @return the median element of values
     */
    public static <T extends Number & Comparable<T>> double median(final Collection<T> values) {
        Utils.nonEmpty(values, "cannot take the median of a collection with no values.");
        return new Median().evaluate(values.stream().mapToDouble(Number::doubleValue).toArray());
    }

    /**
     * Rounds the double to the given number of decimal places.
     * For example, rounding 3.1415926 to 3 places would give 3.142.
     * The requirement is that it works exactly as writing a number down with string.format and reading back in.
     */
    public static double roundToNDecimalPlaces(final double in, final int n) {
         Utils.validateArg(n > 0, "must round to at least one decimal place");
         final double mult = Math.pow(10,n);
         return Math.round( (in+Math.ulp(in))*mult )/mult;
    }

    public static boolean wellFormedDouble(final double val) {
        return !Double.isInfinite(val) && !Double.isNaN(val);
    }

    /**
     * Calculate f(x) = Normal(x | mu = mean, sigma = sd)
     * @param mean the desired mean of the Normal distribution
     * @param sd the desired standard deviation of the Normal distribution
     * @param x the value to evaluate
     * @return a well-formed double
     */
    public static double normalDistribution(final double mean, final double sd, final double x) {
        Utils.validateArg(sd >= 0, "sd: Standard deviation of normal must be >= 0");
        Utils.validateArg(wellFormedDouble(mean) && wellFormedDouble(sd) && wellFormedDouble(x),
                          "mean, sd, or, x : Normal parameters must be well formatted (non-INF, non-NAN)");

        return (INV_SQRT_2_PI / sd) *  Math.exp(-(x - mean)*(x-mean)/ (2.0 * sd * sd));
    }

    /**
     * Return the likelihood of observing the counts of categories having sampled a population
     * whose categorical frequencies are distributed according to a Dirichlet distribution
     * @param params - params of the prior dirichlet distribution
     * @param counts - the counts of observation in each category
     * @return - associated likelihood
     */
    public static double dirichletMultinomial(double[] params, int[] counts) {
        Utils.nonNull(params);
        Utils.nonNull(counts);
        Utils.validateArg(params.length == counts.length, "The number of dirichlet parameters must match the number of categories");
        final double dirichletSum = sum(params);
        final int countSum = (int) sum(counts);
        double prefactor = log10MultinomialCoefficient(countSum,counts) + log10Gamma(dirichletSum) - log10Gamma(dirichletSum+countSum);
        return prefactor + new IndexRange(0, counts.length).sum(n -> log10Gamma(counts[n] + params[n]) - log10Gamma(params[n]));
    }

    /**
     * The following method implements Arrays.stream(array).map(func).toArray(), which is concise but performs poorly due
     * to the overhead of creating a stream, especially with small arrays.  Thus we wrap the wordy but fast array code
     * in the following method which permits concise Java 8 code.
     *
     * Returns a new array -- the original array in not modified.
     *
     * This method has been benchmarked and performs as well as array-only code.
     */
    public static double[] applyToArray(final double[] array, final DoubleUnaryOperator func) {
        Utils.nonNull(func);
        Utils.nonNull(array);
        final double[] result = new double[array.length];
        for (int m = 0; m < result.length; m++) {
            result[m] = func.applyAsDouble(array[m]);
        }
        return result;
    }

    /**
     * The following method implements Arrays.stream(array).map(func).toArray(), which is concise but performs poorly due
     * to the overhead of creating a stream, especially with small arrays.  Thus we wrap the wordy but fast array code
     * in the following method which permits concise Java 8 code.
     *
     * Returns a new array -- the original array in not modified.
     *
     * This method has been benchmarked and performs as well as array-only code.
     */
    public static double[] applyToArray(final int[] array, final IntToDoubleFunction func) {
        Utils.nonNull(func);
        Utils.nonNull(array);
        final double[] result = new double[array.length];
        for (int m = 0; m < result.length; m++) {
            result[m] = func.applyAsDouble(array[m]);
        }
        return result;
    }

    /**
     * The following method implements Arrays.stream(array).map(func).toArray(), which is concise but performs poorly due
     * to the overhead of creating a stream, especially with small arrays.  Thus we wrap the wordy but fast array code
     * in the following method which permits concise Java 8 code.
     *
     * The original array is modified in place.
     *
     * This method has been benchmarked and performs as well as array-only code.
     */
    public static double[] applyToArrayInPlace(final double[] array, final DoubleUnaryOperator func) {
        Utils.nonNull(array);
        Utils.nonNull(func);
        for (int m = 0; m < array.length; m++) {
            array[m] = func.applyAsDouble(array[m]);
        }
        return array;
    }

    /**
     * Test whether all elements of a double[] array satisfy a double -> boolean predicate
     */
    public static boolean allMatch(final double[] array, final DoublePredicate pred) {
        Utils.nonNull(array);
        Utils.nonNull(pred);
        for (final double x : array) {
            if (!pred.test(x)) {
                return false;
            }
        }
        return true;
    }

    /**
     * Test whether all elements of an int[] array satisfy an int -> boolean predicate
     */
    public static boolean allMatch(final int[] array, final IntPredicate pred) {
        Utils.nonNull(array, "array may not be null");
        Utils.nonNull(pred, "predicate may not be null");
        for (final int x : array) {
            if (!pred.test(x)) {
                return false;
            }
        }
        return true;
    }

    /**
     *
     * @param array array of integers
     * @return index of the max. In case of a tie, return the smallest index
     */
    public static int maxElementIndex(final int[] array){
        int maxIndex = 0;
        int currentMax = Integer.MIN_VALUE;
        for (int i = 0; i < array.length; i++){
            if (array[i] > currentMax){
                maxIndex = i;
                currentMax = array[i];
            }
        }
        return maxIndex;
    }

    /**
     *
     * Computes the log10 probability density of BetaBinomial(k|n, alpha, beta)
     *
     * @param alpha pseudocount of number of heads
     * @param beta pseudocount of number of tails
     * @param k value to evaluate
     * @param n number of coin flips
     * @return probability density function evaluated at k
     */
    public static double log10BetaBinomialProbability(final int k, final int n, final double alpha, final double beta){
        Utils.validateArg(k <= n, String.format("k must be less than or equal to n but got k = %d, n = %d", k, n));
        return log10BinomialCoefficient(n, k) + Beta.logBeta(k + alpha, n - k + beta) * LOG10_OF_E -
                Beta.logBeta(alpha, beta) * LOG10_OF_E;
    }

    public static boolean isAProbability(final double p){
        return p >= 0 && p <= 1;
    }
}
