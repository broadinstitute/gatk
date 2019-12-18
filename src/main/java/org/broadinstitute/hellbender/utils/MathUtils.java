package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.distribution.EnumeratedDistribution;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.exception.NumberIsTooLargeException;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.MathArrays;
import org.apache.commons.math3.util.Pair;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.function.*;
import java.util.stream.Collectors;

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
    public static final double LOG10_ONE_THIRD = -Math.log10(3.0);
    public static final double LOG_ONE_THIRD = -Math.log(3.0);
    public static final double INV_LOG_2 = 1.0 / Math.log(2.0);
    private static final double LOG_10 = Math.log(10);
    private static final double INV_LOG_10 = 1.0 / LOG_10;
    public static final double LOG10_E = Math.log10(Math.E);

    private static final double ROOT_TWO_PI = Math.sqrt(2.0 * Math.PI);

    private static final Log10Cache LOG_10_CACHE = new Log10Cache();
    private static final Log10FactorialCache LOG_10_FACTORIAL_CACHE = new Log10FactorialCache();
    private static final DigammaCache DIGAMMA_CACHE = new DigammaCache();

    /**
     * Private constructor.  No instantiating this class!
     */
    private MathUtils() { }

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
     * Given an double array returns the difference between the second smallest element
     * and the smallest element.
     * @param values the input values.
     * @param defaultValue value to return in case the input array has zero or only one element.
     * @return 0 if the input array has less than 2 elements, otherwise the second smallest
     * - the smallest.
     * @throws IllegalArgumentException if {@code values} is {@code null}.
     */
    public static double secondSmallestMinusSmallest(final double[] values, final double defaultValue) {
        Utils.nonNull(values);
        if (values.length <= 1) {
            return defaultValue;
        } else {
            double smallest = values[0];
            double secondSmallest = Double.MAX_VALUE;
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

    public static double[] ebeAdd(final double[] a, final double[] b) throws DimensionMismatchException {
        if (a.length != b.length) {
            throw new DimensionMismatchException(a.length, b.length);
        }

        final double[] result = a.clone();
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

    /**
     * Add two double arrays, modifying the first
     * @param array the array to store the result
     * @param summand   the array to be added
     */
    public static void addToArrayInPlace(final double[] array, final double[] summand) {
        Utils.validateArg(array.length == summand.length, "Arrays must have same length");
        for (int n = 0; n < array.length; n++) {
            array[n] += summand[n];
        }
    }

    /**
     * Add two int arrays, modifying the first
     * @param array the array to store the result
     * @param summand   the array to be added
     */
    public static void addToArrayInPlace(final int[] array, final int[] summand) {
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
     * Composes and array of doubles given the start, end (enclosed*) and the step or interval between consecutive
     * values in the sequence.
     * <p>
     *     For example {@code doubles(1, 10, 1)} results in the sequence {@code  1.0, 2.0, 3.0, ..., 9.0, 10.0}.
     * </p>
     * <p>
     *     It also works with a negative step as long as end < start. For example {@code doubles(3, 1, -.5)} results in
     *     {@code 3.0, 2.5, 2.0, 1.5, 1.0}.
     * </p>
     * <p>
     *     The 'end' values might be included if there is a N >= 1 such that N * step + start == end. There is a
     *     difference "tolerances" three orders of magnitude below thet step absolut value's so if the step is 1. and
     *     the last value is under 0.001 absolute difference with 'end', the 'end' value is used instead.
     * </p>
     * <p>
     *     A sequence request where the difference between the start and end is less than the 'tolerance' described above
     *     results in a single value sequence equal to the start.
     * </p>
     *
     * @param start the first values of the sequence.
     * @param limit the limit of the sequence
     * @param step the step increment (or decrement) between elements of the sequence.
     * @return never {@code null}.
     */
    public static double[] doubles(final double start, final double limit, final double step) {
        ParamUtils.isFinite(start, "the start must be finite");
        ParamUtils.isFinite(limit, "the limit must be finite");
        ParamUtils.isFinite(step, "the step must be finite");
        final double tolerance = Math.pow(10, Math.floor(Math.min(0, Math.log10(Math.abs(step)) - 3)));
        final double diff = limit - start;
        if (Math.abs(diff) < tolerance) {
            return new double[] {start};
        }
        Utils.validateArg(diff * step > 0, "the difference between start and end must have the same sign as the step");
        if (diff == 0) {
            return new double[] {start};
        } else if ((diff > 0) == (step > 0)) {
            final long lengthAsLong = Math.round(Math.floor(1.0 + tolerance + diff / step));
            if (lengthAsLong > Integer.MAX_VALUE) {
                throw new IllegalArgumentException("cannot produce such a large sequence with " + lengthAsLong + " elements");
            }
            final int length = (int) lengthAsLong;
            final double[] result = new double[length];
            for (int i = 0; i < length; i++) {
                result[i] = start + step * i;
            }
            if (Math.abs(result[result.length - 1] - limit) <= tolerance) {
                result[result.length - 1] = limit;
            }
            return result;
        } else {
            throw new IllegalArgumentException("the max - min difference and increment must have the same sign");
        }
    }

    /**
     * Returns an array of doubles filled with a particular value.
     * @param repeats number of elements in the returned array.
     * @param val the fill value.
     * @return never {@code null}.
     */
    public static double[] doubles(final int repeats, final double val) {
        ParamUtils.isPositiveOrZero(repeats, "repeats must be 0 or greater");
        final double[] result = new double[repeats];
        Arrays.fill(result, val);
        return result;
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
     * @param a the input exponent.
     * @return {@link Double#NaN NaN} if {@code a > 0}, otherwise the corresponding value.
     */
    public static double log10OneMinusPow10(final double a) {
        if (a > 0) return Double.NaN;
        if (a == 0) return Double.NEGATIVE_INFINITY;
        final double b = a * LOG_10;
        return NaturalLogUtils.log1mexp(b) * INV_LOG_10;
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
    public static boolean isValidLog10ProbabilityVector(final double[] vector, final int expectedSize, final boolean shouldSumToOne) {
        Utils.nonNull(vector);
        return vector.length == expectedSize && allMatch(vector, MathUtils::isValidLog10Probability) &&
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
        return LOG_10_CACHE.get(i);
    }

    public static double digamma(int i) {
        return DIGAMMA_CACHE.get(i);
    }

    public static double log10sumLog10(final double[] log10values) {
        return log10sumLog10(Utils.nonNull(log10values), 0);
    }

    public static double log10sumLog10(final double[] log10p, final int start) {
        return log10sumLog10(Utils.nonNull(log10p), start, log10p.length);
    }

    public static double log10sumLog10(final double[] log10p, final int start, final int finish) {
        Utils.nonNull(log10p);
        if (finish - start < 2) {
            return finish == start ? Double.NEGATIVE_INFINITY : log10p[start];
        } else {
            final int maxElementIndex = maxElementIndex(log10p, start, finish);
            final double maxValue = log10p[maxElementIndex];
            if (maxValue == Double.NEGATIVE_INFINITY) {
                return maxValue;
            }
            final double sum = 1.0 + new IndexRange(start, finish).sum(i -> i == maxElementIndex ? 0 : Math.pow(10.0, log10p[i] - maxValue));
            Utils.validateArg(!Double.isNaN(sum) && sum != Double.POSITIVE_INFINITY, "log10p values must be non-infinite and non-NAN");
            return maxValue + Math.log10(sum);
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
        return ln * LOG10_E;
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

    public static long sum(final long[] x) {
        Utils.nonNull(x);
        int total = 0;
        for (long v : x)
            total += v;
        return total;
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
        } else if (log10p == 0) {
            return k == n ? 0 : Double.NEGATIVE_INFINITY;
        }
        double log10OneMinusP = Math.log10(1 - Math.pow(10.0, log10p));
        return log10BinomialCoefficient(n, k) + log10p * k + log10OneMinusP * (n - k);
    }

    /**
     * @see #binomialProbability(int, int, double) with p=0.5 and log10 applied to result
     */
    public static double log10BinomialProbability(final int n, final int k) {
        return log10BinomialCoefficient(n, k) + (n * LOG10_ONE_HALF);
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

    /**
     * Do the log-sum trick for three double values.
     * @param a
     * @param b
     * @param c
     * @return the sum... perhaps NaN or infinity if it applies.
     */
    public static double log10SumLog10(final double a, final double b, final double c) {
        if (a >= b && a >= c)  {
            return a + Math.log10(1 + Math.pow(10.0, b - a) + Math.pow(10.0, c - a));
        } else if (b >= c) {
            return b + Math.log10(1 + Math.pow(10.0, a - b) + Math.pow(10.0, c - b));
        } else {
            return c + Math.log10(1 + Math.pow(10.0, a - c) + Math.pow(10.0, b - c));
        }
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
        final double a = -1.0 * Math.log10(sd * ROOT_TWO_PI);
        final double b = -1.0 * (square(x - mean) / (2.0 * square(sd))) / LOG_10;
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
    public static double[] normalizeSumToOne(final double[] array) {
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

    public static int arrayMax(final int[] array) {
        Utils.nonNull(array);
        return array[maxElementIndex(array)];
    }

    /**
     * Returns the maximum value within and int array interval.
     * <p>
     *     A default value must be provided in case the requested interval is empty.
     * </p>
     *
     * @param array the source array.
     * @param from first position to be considered.
     * @param to position after the last to be considered (i.e. exclusive).
     * @param defaultValue the default value to return in case that the interval provided is empty.
     * @throws IndexOutOfBoundsException if the from-to interval indexes point to a invalid index range.
     * @return any integer value is a valid return for this method.
     */
    public static int arrayMax(final int[] array, final int from, final int to, final int defaultValue) {
        if (to > from) {
            int value = array[from];
            for (int i = from + 1; i < to; i++) {
                final int candidate = array[i];
                if (candidate > value) {
                    value = candidate;
                }
            }
            return value;
        } else if (from >= 0) {
            return defaultValue;
        } else {
            throw new ArrayIndexOutOfBoundsException(from);
        }
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

    public static boolean isValidLog10Probability(final double result) { return result <= 0.0; }

    public static boolean isValidProbability(final double result) {
        return result >= 0.0 && result <= 1.0;
    }

    public static double log10ToLog(final double log10){ return log10 * LOG_10; }

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
        return LOG_10_FACTORIAL_CACHE.get(n);
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
     * Compute in a numerically correct way the quantity log10(1-x)
     *
     * Uses the approximation log10(1-x) = log10(1/x - 1) + log10(x) to avoid very quick underflow
     * in 1-x when x is very small
     *
     * log10(1-x) = log10( x * (1-x) / x )
     *            = log10(x) + log10( (1-x) / x)
     *            = log10(x) + log10(1/x - 1)
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

    /**
     * Compute the median of a list of numbers
     *
     * If values.length is odd, this will be the middle value when the elements are sorted
     * If values.length is even then it will be the mean of the two values closest to the middle.
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

        return Math.exp(-(x - mean)*(x-mean)/ (2.0 * sd * sd)) / (sd * ROOT_TWO_PI);
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
     * @param array array of integers
     * @return index of the min. In case of a tie, return the smallest index
     */
    public static int minElementIndex(final int[] array){
        int minIndex = 0;
        int currentMin = Integer.MAX_VALUE;
        for (int i = 0; i < array.length; i++){
            if (array[i] < currentMin){
                minIndex = i;
                currentMin = array[i];
            }
        }
        return minIndex;
    }

    /**
     *
     * @param array array of doubles
     * @return index of the min. In case of a tie, return the smallest index
     */
    public static int minElementIndex(final double[] array){
        int minIndex = 0;
        double currentMin = Double.MAX_VALUE;
        for (int i = 0; i < array.length; i++){
            if (array[i] < currentMin){
                minIndex = i;
                currentMin = array[i];
            }
        }
        return minIndex;
    }

    /**
     * Convert a long to the exact same int value. If it can't be represented exactly as an int, then throw an exception
     * produced by the given supplier.
     */
    public static int toIntExactOrThrow(long toConvert, Supplier<RuntimeException> exceptionSupplier){
        if( toConvert == (int)toConvert ){
            return (int)toConvert;
        } else {
            throw exceptionSupplier.get();
        }
    }

    /**
     * Computes the entropy -p*ln(p) - (1-p)*ln(1-p) of a Bernoulli distribution with success probability p
     * using an extremely fast Pade approximation that is very accurate for all values of 0 <= p <= 1.
     *
     * See http://www.nezumi.demon.co.uk/consult/logx.htm
     */
    public static final double fastBernoulliEntropy(final double p) {
        final double product = p * (1 - p);
        return product * (11 + 33 * product) / (2 + 20 * product);
    }
}
