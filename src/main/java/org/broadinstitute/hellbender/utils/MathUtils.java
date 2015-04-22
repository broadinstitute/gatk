package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.LUDecomposition;
import java.util.*;

/**
 * MathUtils is a static class (no instantiation allowed!) with some useful math methods.
 */
public final class MathUtils {

    /**
     * Private constructor.  No instantiating this class!
     */
    private MathUtils() {
    }

    /**
     * The smallest log10 value we'll emit from normalizeFromLog10 and other functions
     * where the real-space value is 0.0.
     */
    public static final double LOG10_P_OF_ZERO = -1000000.0;
    public static final double NATURAL_LOG_OF_TEN = Math.log(10.0);


    /**
     * A helper class to maintain a cache of log10 values
     */
    public static class Log10Cache {
        /**
         * Get the value of log10(n), expanding the cache as necessary
         * @param n operand
         * @return log10(n)
         */
        public static double get(final int n) {
            if (n < 0)
                throw new GATKException(String.format("Can't take the log of a negative number: %d", n));
            if (n >= cache.length)
                ensureCacheContains(Math.max(n+10, 2*cache.length));
            /*
               Array lookups are not atomic.  It's possible that the reference to cache could be
               changed between the time the reference is loaded and the data is fetched from the correct
               offset.  However, the value retrieved can't change, and it's guaranteed to be present in the
               old reference by the conditional above.
             */
            return cache[n];
        }

        /**
         * Ensures that the cache contains a value for n.  After completion of ensureCacheContains(n),
         * #get(n) is guaranteed to return without causing a cache expansion
         * @param n desired value to be precomputed
         */
        public static void ensureCacheContains(final int n) {
            if (n < cache.length)
                return;
            final double[] newCache = new double[n + 1];
            System.arraycopy(cache, 0, newCache, 0, cache.length);
            for (int i=cache.length; i < newCache.length; i++)
                newCache[i] = Math.log10(i);
            cache = newCache;
        }

        //initialize with the special case: log10(0) = NEGATIVE_INFINITY
        private static double[] cache = new double[] { Double.NEGATIVE_INFINITY };
    }

    // A fast implementation of the Math.round() method.  This method does not perform
    // under/overflow checking, so this shouldn't be used in the general case (but is fine
    // if one is already make those checks before calling in to the rounding).
    public static int fastRound(final double d) {
        return (d > 0.0) ? (int) (d + 0.5d) : (int) (d - 0.5d);
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
     * binomial Probability(int, int, double) with log10 applied to result
     */
    public static double log10BinomialProbability(final int n, final int k, final double log10p) {
        if ( log10p > 1e-18 )
            throw new IllegalArgumentException("log10p: Log-probability must be 0 or less");
        double log10OneMinusP = Math.log10(1 - Math.pow(10, log10p));
        return log10BinomialCoefficient(n, k) + log10p * k + log10OneMinusP * (n - k);
    }

    /**
     * A version of Log10SumLog10 that tolerates NaN values in the array
     *
     * In the case where one or more of the values are NaN, this function returns NaN
     *
     * @param values a non-null vector of doubles
     * @return log10 of the sum of the log10 values, or NaN
     */
    public static double nanTolerantLog10SumLog10(final double[] values) {
        for ( final double value : values ) {
            if (Double.isNaN(value)) {
                return Double.NaN;
            }
        }
        return MathUtils.log10sumLog10(values);
    }

    public static double log10sumLog10(final double[] log10p, final int start) {
        return log10sumLog10(log10p, start, log10p.length);
    }

    public static double log10sumLog10(final double[] log10values) {
        return log10sumLog10(log10values, 0);
    }


    public static double log10sumLog10(final double[] log10p, final int start, final int finish) {

        if (start >= finish)
            return Double.NEGATIVE_INFINITY;
        final int maxElementIndex = MathUtils.maxElementIndex(log10p, start, finish);
        final double maxValue = log10p[maxElementIndex];
        if(maxValue == Double.NEGATIVE_INFINITY)
            return maxValue;
        double sum = 1.0;
        for (int i = start; i < finish; i++) {
            double curVal = log10p[i];
            double scaled_val = curVal - maxValue;
            if (i == maxElementIndex || curVal == Double.NEGATIVE_INFINITY) {
                continue;
            }
            else {
                sum += Math.pow(10.0, scaled_val);
            }
        }
        if ( Double.isNaN(sum) || sum == Double.POSITIVE_INFINITY ) {
            throw new IllegalArgumentException("log10p: Values must be non-infinite and non-NAN");
        }
        return maxValue + (sum != 1.0 ? Math.log10(sum) : 0.0);
    }

    private static final double LOG1MEXP_THRESHOLD = Math.log(0.5);

    private static final double LN_10 = Math.log(10);

    /**
     * normalizes the log10-based array.  ASSUMES THAT ALL ARRAY ENTRIES ARE <= 0 (<= 1 IN REAL-SPACE).
     *
     * @param array             the array to be normalized
     * @param takeLog10OfOutput if true, the output will be transformed back into log10 units
     * @return a newly allocated array corresponding the normalized values in array, maybe log10 transformed
     */
    public static double[] normalizeFromLog10(final double[] array, final boolean takeLog10OfOutput) {
        return normalizeFromLog10(array, takeLog10OfOutput, false);
    }

    public static double arraySum(final double[] arr){
        double sum = 0.0;
        for (double anArr : arr) {
            sum += anArr;
        }
        return sum;
    }
    /**
     * See #normalizeFromLog10 but with the additional option to use an approximation that keeps the calculation always in log-space
     */
    public static double[] normalizeFromLog10(final double[] array, final boolean takeLog10OfOutput, final boolean keepInLogSpace) {
        // for precision purposes, we need to add (or really subtract, since they're
        // all negative) the largest value; also, we need to convert to normal-space.
        final double maxValue = arrayMax(array);

        // we may decide to just normalize in log space without converting to linear space
        if (keepInLogSpace) {
            for (int i = 0; i < array.length; i++) {
                array[i] = array[i] - maxValue;
            }
            return array;
        }

        // default case: go to linear space
        final double[] normalized = new double[array.length];
        for (int i = 0; i < array.length; i++) {
            normalized[i] = Math.pow(10, array[i] - maxValue);
        }

        // normalize
        final double sum = arraySum(normalized);
        for (int i = 0; i < array.length; i++) {
            double x = normalized[i] / sum;
            if (takeLog10OfOutput) {
                x = Math.log10(x);
                if ( x < LOG10_P_OF_ZERO || Double.isInfinite(x) ) {
                    x = array[i] - maxValue;
                }
            }

            normalized[i] = x;
        }

        return normalized;
    }

    /**
     * normalizes the log10-based array.  ASSUMES THAT ALL ARRAY ENTRIES ARE <= 0 (<= 1 IN REAL-SPACE).
     *
     * @param array the array to be normalized
     * @return a newly allocated array corresponding the normalized values in array
     */
    public static double[] normalizeFromLog10(final double[] array) {
        return normalizeFromLog10(array, false);
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

    public static int sum(final List<Integer> list ) {
        int sum = 0;
        for ( Integer i : list ) {
            sum += i;
        }
        return sum;
    }

    public static int countOccurrences(final boolean element, final boolean[] array) {
        int count = 0;
        for (final boolean b : array) {
            if (element == b)
                count++;
        }

        return count;
    }

    /**
     * Checks that the result is a well-formed log10 probability
     *
     * @param result a supposedly well-formed log10 probability value.  By default allows
     *               -Infinity values, as log10(0.0) == -Infinity.
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

    //
    // useful common utility routines
    //

    static public double max(double x0, double x1, double x2) {
        double a = Math.max(x0, x1);
        return Math.max(a, x2);
    }

    /**
     * Converts LN to LOG10
     *
     * @param ln log(x)
     * @return log10(x)
     */
    public static double lnToLog10(final double ln) {
        return ln * Math.log10(Math.E);
    }

    /**
     * Constants to simplify the log gamma function calculation.
     */
    private static final double zero = 0.0, one = 1.0, half = .5, a0 = 7.72156649015328655494e-02, a1 = 3.22467033424113591611e-01, a2 = 6.73523010531292681824e-02, a3 = 2.05808084325167332806e-02, a4 = 7.38555086081402883957e-03, a5 = 2.89051383673415629091e-03, a6 = 1.19270763183362067845e-03, a7 = 5.10069792153511336608e-04, a8 = 2.20862790713908385557e-04, a9 = 1.08011567247583939954e-04, a10 = 2.52144565451257326939e-05, a11 = 4.48640949618915160150e-05, tc = 1.46163214496836224576e+00, tf = -1.21486290535849611461e-01, tt = -3.63867699703950536541e-18, t0 = 4.83836122723810047042e-01, t1 = -1.47587722994593911752e-01, t2 = 6.46249402391333854778e-02, t3 = -3.27885410759859649565e-02, t4 = 1.79706750811820387126e-02, t5 = -1.03142241298341437450e-02, t6 = 6.10053870246291332635e-03, t7 = -3.68452016781138256760e-03, t8 = 2.25964780900612472250e-03, t9 = -1.40346469989232843813e-03, t10 = 8.81081882437654011382e-04, t11 = -5.38595305356740546715e-04, t12 = 3.15632070903625950361e-04, t13 = -3.12754168375120860518e-04, t14 = 3.35529192635519073543e-04, u0 = -7.72156649015328655494e-02, u1 = 6.32827064025093366517e-01, u2 = 1.45492250137234768737e+00, u3 = 9.77717527963372745603e-01, u4 = 2.28963728064692451092e-01, u5 = 1.33810918536787660377e-02, v1 = 2.45597793713041134822e+00, v2 = 2.12848976379893395361e+00, v3 = 7.69285150456672783825e-01, v4 = 1.04222645593369134254e-01, v5 = 3.21709242282423911810e-03, s0 = -7.72156649015328655494e-02, s1 = 2.14982415960608852501e-01, s2 = 3.25778796408930981787e-01, s3 = 1.46350472652464452805e-01, s4 = 2.66422703033638609560e-02, s5 = 1.84028451407337715652e-03, s6 = 3.19475326584100867617e-05, r1 = 1.39200533467621045958e+00, r2 = 7.21935547567138069525e-01, r3 = 1.71933865632803078993e-01, r4 = 1.86459191715652901344e-02, r5 = 7.77942496381893596434e-04, r6 = 7.32668430744625636189e-06, w0 = 4.18938533204672725052e-01, w1 = 8.33333333333329678849e-02, w2 = -2.77777777728775536470e-03, w3 = 7.93650558643019558500e-04, w4 = -5.95187557450339963135e-04, w5 = 8.36339918996282139126e-04, w6 = -1.63092934096575273989e-03;

    /**
     * Efficient rounding functions to simplify the log gamma function calculation
     * double to long with 32 bit shift
     */
    private static final int HI(final double x) {
        return (int) (Double.doubleToLongBits(x) >> 32);
    }

    /**
     * Efficient rounding functions to simplify the log gamma function calculation
     * double to long without shift
     */
    private static final int LO(final double x) {
        return (int) Double.doubleToLongBits(x);
    }

    /**
     * Most efficent implementation of the lnGamma (FDLIBM)
     * Use via the log10Gamma wrapper method.
     */
    @SuppressWarnings("fallthrough")
    private static double lnGamma(final double x) {
        double t, y, z, p, p1, p2, p3, q, r, w;
        int i;

        int hx = HI(x);
        int lx = LO(x);

        /* purge off +-inf, NaN, +-0, and negative arguments */
        int ix = hx & 0x7fffffff;
        if (ix >= 0x7ff00000)
            return Double.POSITIVE_INFINITY;
        if ((ix | lx) == 0 || hx < 0)
            return Double.NaN;
        if (ix < 0x3b900000) {    /* |x|<2**-70, return -log(|x|) */
            return -Math.log(x);
        }

        /* purge off 1 and 2 */
        if ((((ix - 0x3ff00000) | lx) == 0) || (((ix - 0x40000000) | lx) == 0))
            r = 0;
            /* for x < 2.0 */
        else if (ix < 0x40000000) {
            if (ix <= 0x3feccccc) {     /* lgamma(x) = lgamma(x+1)-log(x) */
                r = -Math.log(x);
                if (ix >= 0x3FE76944) {
                    y = one - x;
                    i = 0;
                }
                else if (ix >= 0x3FCDA661) {
                    y = x - (tc - one);
                    i = 1;
                }
                else {
                    y = x;
                    i = 2;
                }
            }
            else {
                r = zero;
                if (ix >= 0x3FFBB4C3) {
                    y = 2.0 - x;
                    i = 0;
                } /* [1.7316,2] */
                else if (ix >= 0x3FF3B4C4) {
                    y = x - tc;
                    i = 1;
                } /* [1.23,1.73] */
                else {
                    y = x - one;
                    i = 2;
                }
            }

            switch (i) {
                case 0:
                    z = y * y;
                    p1 = a0 + z * (a2 + z * (a4 + z * (a6 + z * (a8 + z * a10))));
                    p2 = z * (a1 + z * (a3 + z * (a5 + z * (a7 + z * (a9 + z * a11)))));
                    p = y * p1 + p2;
                    r += (p - 0.5 * y);
                    break;
                case 1:
                    z = y * y;
                    w = z * y;
                    p1 = t0 + w * (t3 + w * (t6 + w * (t9 + w * t12)));    /* parallel comp */
                    p2 = t1 + w * (t4 + w * (t7 + w * (t10 + w * t13)));
                    p3 = t2 + w * (t5 + w * (t8 + w * (t11 + w * t14)));
                    p = z * p1 - (tt - w * (p2 + y * p3));
                    r += (tf + p);
                    break;
                case 2:
                    p1 = y * (u0 + y * (u1 + y * (u2 + y * (u3 + y * (u4 + y * u5)))));
                    p2 = one + y * (v1 + y * (v2 + y * (v3 + y * (v4 + y * v5))));
                    r += (-0.5 * y + p1 / p2);
            }
        }
        else if (ix < 0x40200000) {             /* x < 8.0 */
            i = (int) x;
            t = zero;
            y = x - (double) i;
            p = y * (s0 + y * (s1 + y * (s2 + y * (s3 + y * (s4 + y * (s5 + y * s6))))));
            q = one + y * (r1 + y * (r2 + y * (r3 + y * (r4 + y * (r5 + y * r6)))));
            r = half * y + p / q;
            z = one;    /* lgamma(1+s) = log(s) + lgamma(s) */
            switch (i) {
                case 7:
                    z *= (y + 6.0);    /* FALLTHRU */
                case 6:
                    z *= (y + 5.0);    /* FALLTHRU */
                case 5:
                    z *= (y + 4.0);    /* FALLTHRU */
                case 4:
                    z *= (y + 3.0);    /* FALLTHRU */
                case 3:
                    z *= (y + 2.0);    /* FALLTHRU */
                    r += Math.log(z);
                    break;
            }
            /* 8.0 <= x < 2**58 */
        }
        else if (ix < 0x43900000) {
            t = Math.log(x);
            z = one / x;
            y = z * z;
            w = w0 + z * (w1 + y * (w2 + y * (w3 + y * (w4 + y * (w5 + y * w6)))));
            r = (x - half) * (t - one) + w;
        }
        else
            /* 2**58 <= x <= inf */
            r = x * (Math.log(x) - one);
        return r;
    }

    /**
     * Calculates the log10 of the gamma function for x using the efficient FDLIBM
     * implementation to avoid overflows and guarantees high accuracy even for large
     * numbers.
     *
     * @param x the x parameter
     * @return the log10 of the gamma function at x.
     */
    public static double log10Gamma(final double x) {
        return lnToLog10(lnGamma(x));
    }


    public static double log10Factorial(final int x) {
        if (x >= Log10FactorialCache.size() || x < 0)
            return log10Gamma(x + 1);
        else
            return Log10FactorialCache.get(x);
    }

    /**
     * Wrapper class so that the log10Factorial array is only calculated if it's used
     */
    private static class Log10FactorialCache {

        /**
         * The size of the precomputed cache.  Must be a positive number!
         */
        private static final int CACHE_SIZE = 10_000;

        public static int size() { return CACHE_SIZE; }

        public static double get(final int n) {
            if (cache == null)
                initialize();
            return cache[n];
        }

        private static void initialize() {
            if (cache == null) {
                Log10Cache.ensureCacheContains(CACHE_SIZE);
                cache = new double[CACHE_SIZE];
                cache[0] = 0.0;
                for (int k = 1; k < cache.length; k++)
                    cache[k] = cache[k-1] + Log10Cache.get(k);
            }
        }

        private static double[] cache = null;
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
     * Calculate f(x) = log10 ( Normal(x | mu = mean, sigma = sd) )
     * @param mean the mean of the Normal distribution
     * @param sd the standard deviation of the Normal distribution
     * @param x the value to evaluate
     * @return a well-formed double
     */
    public static double normalDistributionLog10(final double mean, final double sd, final double x) {
        return new NormalDistribution(mean, sd).logDensity(x)/NATURAL_LOG_OF_TEN;
    }

    public static RealMatrix identityMatrix(int n){
        double[] ones = new double[n];
        Arrays.fill(ones, 1);
        return new DiagonalMatrix(ones);
    }

}
