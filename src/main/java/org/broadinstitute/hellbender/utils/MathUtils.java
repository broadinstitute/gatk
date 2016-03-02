package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.exception.NumberIsTooLargeException;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.Arrays;
import java.util.Collection;

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

    private static final double LN_10 = Math.log(10);

    private static final double LOG1MEXP_THRESHOLD = Math.log(0.5);

    /**
     * Log10 of the e constant.
     */
    public static final double LOG10_OF_E = Math.log10(Math.E);
    public static final double FAIR_BINOMIAL_PROB_LOG10_0_5 = Math.log10(0.5);

    /**
     * Private constructor.  No instantiating this class!
     */
    private MathUtils() {
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
        Utils.validateArg(max > min - 1, "invalid arguments min:" + min + " max:" + max);
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
     * Creates a new sample of k ints from [0..n-1], without duplicates.
     * @throws NumberIsTooLargeException if {@code k > n}.
     * @throws NotStrictlyPositiveException if {@code k <= 0}.
     */
    public static int[] sampleIndicesWithoutReplacement(final int n, final int k) {
        //No error checking : RandomDataGenetator.nextPermutation does it
        return Utils.getRandomDataGenerator().nextPermutation(n, k);
    }

    /**
     * Calculates {@code log10(1-10^a)} without loosing precision.
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
     * Calculates {@code log(1-exp(a))} without loosing precision.
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
        if ( vector.length != expectedSize ) return false;

        for ( final double pr : vector ) {
            if ( ! goodLog10Probability(pr) )
                return false;
        }

        if ( shouldSumToOne && compareDoubles(sumLog10(vector), 1.0, 1e-4) != 0 )
            return false;

        return true; // everything is good
    }

    /**
     *  Returns the sum of values whose log10s we have. That is, returns sum(10^x_i).
     */
    public static double sumLog10(final double[] log10values) {
        return Math.pow(10.0, log10SumLog10(log10values));
    }

    /** Compute Z=X-Y for two numeric vectors X and Y
     *
     * @param x                 First vector
     * @param y                 Second vector
     * @return Vector of same length as x and y so that z[k] = x[k]-y[k]
     */
    public static int[] vectorDiff(final int[]x, final int[] y) {
        Utils.nonNull(x, "x is null");
        Utils.nonNull(y, "y is null");
        if (x.length != y.length)
            throw new IllegalArgumentException("BUG: Lengths of x and y must be the same");

        final int[] result = new int[x.length];
        for (int k=0; k <x.length; k++) {
            result[k] = x[k] - y[k];
        }

        return result;
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
        if ( n < 0 )
            throw new IllegalArgumentException("n: Must have non-negative number of trials");
        double denominator = 0.0;
        int sum = 0;
        for (int x : k) {
            if ( x < 0 )
                throw new IllegalArgumentException("x element of k: Must have non-negative observations of group");
            if ( x > n )
                throw new IllegalArgumentException("x element of k, n: Group observations must be bounded by k");
            denominator += log10Factorial(x);
            sum += x;
        }
        if ( sum != n )
            throw new IllegalArgumentException("k and n: Sum of observations in multinomial must sum to total number of trials");
        return log10Factorial(n) - denominator;
    }

    public static double log10(int i) {
        return log10Cache.get(i);
    }

    /**
     * A helper class to maintain a cache of log values.
     * The cache is immutable after creation.
     */
    private static final class Log10Cache {

        private final double[] cache;

        public Log10Cache(final int capacity) {
            cache = new double[capacity + 1];
            cache[0] = Double.NEGATIVE_INFINITY;    //initialize with the special case: log(0) = NEGATIVE_INFINITY
            for (int i = 1; i < cache.length; i++) {
                cache[i] = Math.log10(i);
            }
        }

        /**
         * Get the value of log10(i), fetching it from the cache or computing it afresh
         * @param i operand
         * @return log10(i)
         */
        public double get(final int i) {
            if (i < 0) {
                throw new IllegalArgumentException(String.format("Can't take the log of a negative number: %d", i));
            }
            if (i >= cache.length) {
                return Math.log10(i);
            }
            return cache[i];
        }

        public int size() {
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

        public static double get(final double difference) {
            if (cache == null) {
                initialize();
            }
            final int index = fastRound(difference * INV_STEP);
            return cache[index];
        }

        private static void initialize() {
            if (cache == null) {
                final int tableSize = (int) (MAX_TOLERANCE / TABLE_STEP) + 1;
                cache = new double[tableSize];
                for (int k = 0; k < cache.length; k++) {
                    cache[k] = Math.log10(1.0 + Math.pow(10.0, -((double) k) * TABLE_STEP));
                }
            }
        }

        //  Phred scores Q and Q+1 differ by 0.1 in their corresponding log-10 probabilities, and by
        // 0.1 * log(10) in natural log probabilities.  Setting TABLE_STEP to an exact divisor of this
        // quantity ensures that approximateSumLog in fact caches exact values for integer phred scores
        private static final double TABLE_STEP = 0.0001;
        private static final double INV_STEP = 1.0 / TABLE_STEP;
        private static double[] cache = null;
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
        return approximateLog10SumLog10(vals, vals.length);
    }

    public static double approximateLog10SumLog10(final double[] vals, final int endIndex) {

        final int maxElementIndex = MathUtils.maxElementIndex(vals, endIndex);
        double approxSum = vals[maxElementIndex];

        for (int i = 0; i < endIndex; i++) {
            if (i == maxElementIndex || vals[i] == Double.NEGATIVE_INFINITY) {
                continue;
            }

            final double diff = approxSum - vals[i];
            if (diff < JacobianLogTable.MAX_TOLERANCE) {
                // See notes from the 2-inout implementation below
                approxSum += JacobianLogTable.get(diff);
            }
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
        }

        if (a == Double.NEGATIVE_INFINITY) {
            return b;
        }

        final double diff = b - a;
        if (diff >= JacobianLogTable.MAX_TOLERANCE) {
            return b;
        }

        // OK, so |b-a| < tol
        // we need to compute log(e^a + e^b) = log(e^b(1 + e^(a-b))) = b + log(1 + e^(-(b-a)))
        // we compute the second term as a table lookup with integer quantization
        // we have pre-stored correction for 0,0.1,0.2,... 10.0
        return b + JacobianLogTable.get(diff);
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
        double s = 0.0;
        for (double v : values)
            s += v;
        return s;
    }

    public static long sum(final int[] x) {
        long total = 0;
        for (int v : x)
            total += v;
        return total;
    }

    public static int sum(final byte[] x) {
        int total = 0;
        for (byte v : x)
            total += (int)v;
        return total;
    }

    public static long sum(final long[] x) {
        int total = 0;
        for (long v : x)
            total += v;
        return total;
    }

    /** Returns the sum of the elements in the array starting with start and ending before stop. */
    public static long sum(final long[] arr, final int start, final int stop) {
        return sum(Arrays.copyOfRange(arr, start, stop));
    }

    /** Returns the sum of the elements in the array starting with start and ending before stop. */
    public static double sum(final double[] arr, final int start, final int stop) {
        Utils.nonNull(arr);
        Utils.validateArg(start <= stop, start + " > " + stop);
        Utils.validateArg(start >= 0, start + " < " + 0);
        Utils.validateArg(stop <= arr.length, stop + " >  " + arr.length);
        return sum(Arrays.copyOfRange(arr, start, stop));
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
     */
    public static double log10BinomialCoefficient(final int n, final int k) {
        if ( n < 0 ) {
            throw new IllegalArgumentException("n: Must have non-negative number of trials");
        }
        if ( k > n || k < 0 ) {
            throw new IllegalArgumentException("k: Must have non-negative number of successes, and no more successes than number of trials");
        }

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
        if ( log10p > 1e-18 )
            throw new IllegalArgumentException("log10p: Log10-probability must be 0 or less");
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
        return log10SumLog10(log10Values, start, log10Values.length);
    }

    public static double log10SumLog10(final double[] log10Values) {
        return log10SumLog10(log10Values, 0);
    }

    public static double log10SumLog10(final double[] log10Values, final int start, final int finish) {
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
            double curVal = log10Values[i];
            double scaled_val = curVal - maxValue;
            if (i == maxElementIndex || curVal == Double.NEGATIVE_INFINITY) {
                continue;
            } else {
                sum += Math.pow(10.0, scaled_val);
            }
        }
        if ( Double.isNaN(sum) || sum == Double.POSITIVE_INFINITY ) {
            throw new IllegalArgumentException("log10 p: Values must be non-infinite and non-NAN");
        }
        return maxValue + (sum != 1.0 ? Math.log10(sum) : 0.0);
    }

    /**
     * normalizes the log-probability array.  ASSUMES THAT ALL ARRAY ENTRIES ARE <= 0 (<= 1 IN REAL-SPACE).
     *
     * @param array the array to be normalized
     * @return a newly allocated array corresponding the normalized values in array
     */
    public static double[] normalizeFromLog10(final double[] array) {
        return normalizeFromLog10(array, false);
    }

    /**
     * normalizes the log-probability array.  ASSUMES THAT ALL ARRAY ENTRIES ARE <= 0 (<= 1 IN REAL-SPACE).
     *
     * @param array             the array to be normalized
     * @param takeLog10OfOutput if true, the output will be transformed back into log units
     * @return a newly allocated array corresponding the normalized values in array, maybe log transformed
     */
    public static double[] normalizeFromLog10(final double[] array, final boolean takeLog10OfOutput) {
        return normalizeFromLog10(array, takeLog10OfOutput, false);
    }


    /**
     * See #normalizeFromLog but with the additional option to use an approximation that keeps the calculation always in log-space
     *
     * @param array
     * @param takeLog10OfOutput
     * @param keepInLogSpace
     *
     * @return
     */
    public static double[] normalizeFromLog10(final double[] array, final boolean takeLog10OfOutput, final boolean keepInLogSpace) {
        // for precision purposes, we need to add (or really subtract, since they're
        // all negative) the largest value; also, we need to convert to normal-space.
        double maxValue = arrayMax(array);

        // we may decide to just normalize in log space without converting to linear space
        if (keepInLogSpace) {
            for (int i = 0; i < array.length; i++) {
                array[i] -= maxValue;
            }
            return array;
        }

        // default case: go to linear space
        double[] normalized = new double[array.length];

        for (int i = 0; i < array.length; i++)
            normalized[i] = Math.pow(10.0, array[i] - maxValue);

        // normalize
        double sum = 0.0;
        for (int i = 0; i < array.length; i++)
            sum += normalized[i];
        for (int i = 0; i < array.length; i++) {
            double x = normalized[i] / sum;
            if (takeLog10OfOutput) {
                x = Math.log10(x);
                if ( x < LOG10_P_OF_ZERO || Double.isInfinite(x) )
                    x = array[i] - maxValue;
            }

            normalized[i] = x;
        }

        return normalized;
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
        if ( array.length == 0 )
            return array;

        final double sum = sum(array);
        final double[] normalized = new double[array.length];
        if ( sum < 0.0 ) throw new IllegalArgumentException("Values in probability array sum to a negative number " + sum);
        for ( int i = 0; i < array.length; i++ ) {
            normalized[i] = array[i] / sum;
        }
        return normalized;
    }

    public static int minElementIndex(final int[] array) {
        Utils.nonNull(array);
        Utils.validateArg(array.length != 0, "Array cannot be empty");

        int minI = 0;
        for (int i = 1; i < array.length; i++) {
            if (array[i] < array[minI])
                minI = i;
        }

        return minI;
    }

    public static int arrayMin(final int[] array) {
        return array[minElementIndex(array)];
    }

    public static int maxElementIndex(final double[] array) {
        return maxElementIndex(array, array.length);
    }

    public static int maxElementIndex(final double[] array, final int start, final int endIndex) {
        if (array == null || array.length == 0)
            throw new IllegalArgumentException("Array cannot be null!");

        if (start > endIndex) {
            throw new IllegalArgumentException("Start cannot be after end.");
        }

        int maxI = start;
        for (int i = (start+1); i < endIndex; i++) {
            if (array[i] > array[maxI])
                maxI = i;
        }
        return maxI;
    }

    public static int maxElementIndex(final double[] array, final int endIndex) {
        return maxElementIndex(array, 0, endIndex);
    }

    public static double arrayMax(final double[] array) {
        return array[maxElementIndex(array)];
    }

    /**
     * Checks that the result is a well-formed log10 probability
     *
     * @param result a supposedly well-formed log10 probability value.  By default allows
     *               -Infinity values, as log(0.0) == -Infinity.
     * @return true if result is really well formed
     */
    public static boolean goodLog10Probability(final double result) {
        return goodLog10Probability(result, true);
    }

    /**
     * Checks that the result is a well-formed log10 probability
     *
     * @param result a supposedly well-formed log10 probability value
     * @param allowNegativeInfinity should we consider a -Infinity value ok?
     * @return true if result is really well formed
     */
    public static boolean goodLog10Probability(final double result, final boolean allowNegativeInfinity) {
        return result <= 0.0 && result != Double.POSITIVE_INFINITY && (allowNegativeInfinity || result != Double.NEGATIVE_INFINITY) && ! Double.isNaN(result);
    }

    /**
     * Checks that the result is a well-formed probability
     *
     * @param result a supposedly well-formed probability value
     * @return true if result is really well formed
     */
    public static boolean goodProbability(final double result) {
        return result >= 0.0 && result <= 1.0 && ! Double.isInfinite(result) && ! Double.isNaN(result);
    }


    /**
     * The size of the precomputed cache of logs.
     * The caches are immutable after creation and so it's no big deal that they are static.
     */
    private static final int PRECOMPUTED_LOGS = 10_000;
    private static final Log10Cache log10Cache = new Log10Cache(PRECOMPUTED_LOGS);
    private static final Log10FactorialCache log10FactorialCache = new Log10FactorialCache(log10Cache);

    public static double log10ToLog(final double log10){
        return log10 * LN_10;
    }

    /**
     * Converts LN to LOG10
     * @param ln log(x)
     * @return log10(x)
     */
     public static double lnToLog10(final double ln) {
        return ln * LOG10_OF_E;
     }

    /**
      * Calculates the log10 of the gamma function for x.
      *
      * @param x the x parameter
      * @return the log10 of the gamma function at x.
      */
    public static double log10Gamma(final double x) {
       return lnToLog10(Gamma.logGamma(x));
    }

    public static double log10Factorial(final int x) {
       if (x >= log10FactorialCache.size() || x < 0)
          return log10Gamma(x + 1);
       else
          return log10FactorialCache.get(x);
    }

    /**
     * Converts a real space array of numbers (typically probabilities) into a log10 array
     *
     * @param prRealSpace
     * @return
     */
    public static double[] toLog10(final double[] prRealSpace) {
        final double[] log10s = new double[prRealSpace.length];
        for (int i = 0; i < prRealSpace.length; i++) {
            log10s[i] = Math.log10(prRealSpace[i]);
        }
        return log10s;
    }

    /**
     * Wrapper class so that the log10Factorial array is only calculated if it's used
     */
    private static class Log10FactorialCache {

        private final double[] cache;

        public Log10FactorialCache(final Log10Cache logCache) {
            cache = new double[logCache.size()];
            cache[0] = 0.0;
            for (int k = 1; k < cache.length; k++) {
                cache[k] = cache[k - 1] + logCache.get(k);
            }
        }

        public int size() { return cache.length; }

        /**
         * Retrieves the precomputed result or computes it afresh.
         * @return log of factorial.
         */
        public double get(final int n) {
            if (n >= size() || n < 0) {
                return log10Gamma(n + 1);
            } else {
                return cache[n];
            }
        }
    }

    /**
     * Compute in a numerical correct way the quantity log10(1-x)
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

    /**
     * Now for some matrix methods
     */

    /**
     *
     * @param m a real-valued matrix
     * @return whether m is symmetric
     */
    public static boolean isSymmetric(RealMatrix m) {
        return m.equals(m.transpose());
    }

    /**
     *
     * @param m a real-valued matrix
     * @return whether m is positive semi-definite i.e. has no negative eigenvalues
     */
    public static boolean isPositiveSemiDefinite(RealMatrix m) {
        EigenDecomposition ed = new EigenDecomposition(m);
        for (final double eigval : ed.getRealEigenvalues()) {
            if (eigval < 0) return false;
        }
        return true;
    }

    /**
     * Compute the logarithm of a square matrix.  Unfortunately, Aoache Commons does not have this method.
     *
     * We compute the matrix logarithm by diagonalizing, taking logarithms of the diagonal entries, and
     * reversing the diagonalizing change of basis
     *
     * @param M
     * @return the matrix logarithm of M
     */
    public static RealMatrix matrixLog(RealMatrix M) {
        EigenDecomposition ed = new EigenDecomposition(M);
        RealMatrix D = ed.getD();   //D is diagonal
        RealMatrix V = ed.getV();   //M = V*D*V^T; V is the diagonalizing change of basis

        //replace D (in-place) by its logarithm
        for (int i = 0; i < M.getColumnDimension(); i++) {
            D.setEntry(i, i, Math.log(D.getEntry(i, i)));
        }

        return V.multiply(D).multiply(V.transpose());   //reverse the change of basis
    }

    /**
     * Measure the difference between two covariance matrices in terms of the Kullback-Leibler
     * divergence between associated Gaussians.
     *
     * If d is the dimension of these matrices, the KL divergence between zero-centered Gaussians
     * with covariances A and B is (1/2){tr[A^(-1)B] + ln(det(A) - ln(det(B)) - d}.  Note: the KL
     * divergence is not symmetric.  Switching A <--> B and averaging gives (1/2){tr[A^(-1)B] + tr[B^(-1)A] - d}
     *
     * @param cov1 a matrix covariance
     * @param cov2 a matrix covariance
     * @return the average of KL divergences, (KL(p|q) + KL(q|p))/2, where p and q are probability densities
     * of zero-centered Gaussians with the give covariance
     */
    public static double covarianceKLDivergence(RealMatrix cov1, RealMatrix cov2) {
        if (!isSymmetric(cov1) || !isSymmetric(cov2)) {
            throw new GATKException("Covariance matrices must be symmetric.");
        }

        if (!isPositiveSemiDefinite(cov1) || !isPositiveSemiDefinite(cov2)) {
            throw new GATKException("Covariance matrices must be positive semidefinite.");
        }

        int d = cov1.getRowDimension();

        if (cov1.getRowDimension() != cov2.getRowDimension()) {
            throw new GATKException("Can only compare covariance matrices of equal dimension.");
        }

        LUDecomposition LU1 = new LUDecomposition(cov1);
        LUDecomposition LU2 = new LUDecomposition(cov2);

        return (LU1.getSolver().solve(cov2).getTrace() + LU2.getSolver().solve(cov1).getTrace() - d)/2;
    }

    /**
     * Measure the geodesic distance between the two covariances within the manifold of symmetric,
     * positive-definite matrices.  This is also called the affine-invariant metric.
     *
     * The formula is ||log(A^(-1/2)*B*A^(-1/2)||_F, where ||    ||_F is the Frobenius norm.  This formula
     * is symmetric despite its appearance.
     *
     * For positive semidefinite matrices with eigendecomposition M = V*D*V^(-1), where D is diagonal
     * the matrix inverse square root is M^(-1/2) = V*D^(-1/2)*V^(-1)
     *
     * @param cov1 a covariance matrix
     * @param cov2 a covariance matrix
     * @return the geodesic distance between cov1 and cov2 in the manifold of positive semi-definite
     * symmetric matrices, which is more natural than the Euclidean distance inherited from the embedding
     * in R^(d^2)
     */
    public static double covarianceGeodesicDistance(RealMatrix cov1, RealMatrix cov2) {
        if (!isSymmetric(cov1) || !isSymmetric(cov2)) {
            throw new GATKException("Covariance matrices must be symmetric.");
        }

        if (!isPositiveSemiDefinite(cov1) || !isPositiveSemiDefinite(cov2)) {
            throw new GATKException("Covariance matrices must be positive semidefinite.");
        }

        if (cov1.getRowDimension() != cov2.getRowDimension()) {
            throw new GATKException("Can only compare covariance matrices of equal dimension.");
        }

        RealMatrix sqrt = (new EigenDecomposition(cov1)).getSquareRoot();
        RealMatrix inverseSqrt = (new LUDecomposition(sqrt)).getSolver().getInverse();

        //the thing inside the matrix logarithm
        RealMatrix mat = inverseSqrt.multiply(cov2).multiply(inverseSqrt);
        return matrixLog(mat).getFrobeniusNorm();

    }

    /** Calculate the mean of an array of doubles. */
    public static double mean(final double[] in, final int start, final int stop) {
        if ((stop - start) <= 0 ) return Double.NaN;

        double total = 0;
        for (int i = start; i < stop; ++i) {
            total += in[i];
        }

        return total / (stop - start);
    }

    /** Calculate the (population) standard deviation of an array of doubles. */
    public static double stddev(final double[] in, final int start, final int length) {
        return stddev(in, start, length, mean(in, start, length));
    }

    /** Calculate the (population) standard deviation of an array of doubles. */
    public static double stddev(final double[] in, final int start, final int stop, final double mean) {
        if ((stop - start) <= 0) return Double.NaN;

        double total = 0;
        for (int i = start; i < stop; ++i) {
            total += (in[i] * in[i]);
        }

        return Math.sqrt((total / (stop - start)) - (mean * mean));
    }

    /** "Promotes" an int[] into a double array with the same values (or as close as precision allows). */
    public static double[] promote(final int[] is) {
        final double[] ds = new double[is.length];
        for (int i = 0; i < is.length; ++i) ds[i] = is[i];
        return ds;
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
        DescriptiveStatistics stats = new DescriptiveStatistics();
        values.stream().mapToDouble(Number::doubleValue).forEach(stats::addValue);
        return stats.apply(new Median());
    }

    /**
     * Rounds the double to the given number of decimal places.
     * For example, rounding 3.1415926 to 3 places would give 3.142.
     * The requirement is that it works exactly as writing a number down with string.format and reading back in.
     */
    public static double roundToNDecimalPlaces(final double in, final int n) {
        if (n < 1) {
            throw new IllegalArgumentException("cannot round to " + n + " decimal places");
        }

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

        double a = 1.0 / (sd * Math.sqrt(2.0 * Math.PI));
        double b = Math.exp(-1.0 * (Math.pow(x - mean, 2.0) / (2.0 * sd * sd)));
        return a * b;
    }
}
