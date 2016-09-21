package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.exception.NumberIsTooLargeException;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;

import java.util.Arrays;
import java.util.Collection;
import java.util.function.*;

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

    public static final double INV_SQRT_2_PI = 1.0 / Math.sqrt(2.0 * Math.PI);

    private static final double NATURAL_LOG_OF_TEN = Math.log(10.0);

    private static final double SQUARE_ROOT_OF_TWO_TIMES_PI = Math.sqrt(2.0 * Math.PI);

    /**
     * Log10 of the e constant.
     */
    public static final double LOG10_OF_E = Math.log10(Math.E);
    public static final double FAIR_BINOMIAL_PROB_LOG10_0_5 = Math.log10(0.5);

    /**
     * Private constructor.  No instantiating this class!
     */
    private MathUtils() { }


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
        return vector.length == expectedSize &&
                allMatch(vector, MathUtils::goodLog10Probability) &&
                !( shouldSumToOne && compareDoubles(sumLog10(vector), 1.0, 1e-4) != 0 );
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
        return log10sumLog10(log10values, 0);
    }

    public static double log10sumLog10(final double[] log10p, final int start) {
        return log10sumLog10(log10p, start, log10p.length);
    }

    public static double log10sumLog10(final double[] log10p, final int start, final int finish) {
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
        return approximateLog10SumLog10(vals, vals.length);
    }

    public static double approximateLog10SumLog10(final double[] vals, final int endIndex) {
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
        Utils.validateArg(start <= stop, () -> start + " > " + stop);
        Utils.validateArg(start >= 0, () -> start + " < " + 0);
        Utils.validateArg(stop <= arr.length, () -> stop + " >  " + arr.length);
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

    /**
     * Calculate f(x) = log10 ( Normal(x | mu = mean, sigma = sd) )
     * @param mean the desired mean of the Normal distribution
     * @param sd the desired standard deviation of the Normal distribution
     * @param x the value to evaluate
     * @return a well-formed double
     */
    public static double normalDistributionLog10(final double mean, final double sd, final double x) {
        if( sd < 0 )
            throw new IllegalArgumentException("sd: Standard deviation of normal must be >0");
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
        double dist = 0.0;
        for (int iii = 0; iii < x.length; iii++) {
            dist += (x[iii] - y[iii]) * (x[iii] - y[iii]);
        }
        return dist;
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
     * @param keepInLogSpace    if true, we don't normalize in the sense of probabilities summing to one but rather in the
     *                          sense of numericla stability where we scale entries so that the largest is 0 in log space
     *
     * @return
     */
    public static double[] normalizeFromLog10(final double[] array, final boolean takeLog10OfOutput, final boolean keepInLogSpace) {
        if (keepInLogSpace) {
            double maxValue = arrayMax(array);
            return applyToArray(array, x -> x - maxValue);
        } else {
            final double log10Sum = log10SumLog10(array);
            final double[] result = applyToArray(array, x -> x - log10Sum);
            return takeLog10OfOutput ? result : applyToArrayInPlace(result, x -> Math.pow(10.0, x));
        }
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
        Utils.validateArg(sum >= 0.0, () -> "Values in probability array sum to a negative number " + sum);
        return applyToArray(array, x -> x/sum);
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
        Utils.nonNull(array, "array may not be null");
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
        return stop <= start ? Double.NaN : Arrays.stream(in, start, stop).average().getAsDouble();
    }

    /** "Promotes" an int[] into a double array with the same values (or as close as precision allows). */
    public static double[] promote(final int[] is) {
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
        Utils.nonNull(func, "function may not be null");
        Utils.nonNull(array, "array may not be null");
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
        Utils.nonNull(array, "array may not be null");
        Utils.nonNull(func, "function may not be null");
        for (int m = 0; m < array.length; m++) {
            array[m] = func.applyAsDouble(array[m]);
        }
        return array;
    }

    /**
     * Test whether all elements of a double[] array satisfy a double -> boolean predicate
     */
    public static boolean allMatch(final double[] array, final DoublePredicate pred) {
        Utils.nonNull(array, "array may not be null");
        Utils.nonNull(pred, "predicate may not be null");
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
}
