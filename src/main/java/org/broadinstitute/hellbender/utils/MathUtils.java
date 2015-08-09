package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.special.Gamma;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.Arrays;

/**
 * MathUtils is a static class (no instantiation allowed!) with some useful math methods.
 */
public final class MathUtils {

    /**
     * The smallest log value we'll emit from normalizeFromLog and other functions
     * where the real-space value is 0.0.
     */
    public static final double LOG_P_OF_ZERO = -1000000.0;

    /**
     * Private constructor.  No instantiating this class!
     */
    private MathUtils() {
    }

    /**
     * A helper class to maintain a cache of log values
     */
    public static final class LogCache {
        /**
         * Get the value of log(n), expanding the cache as necessary
         * @param n operand
         * @return log(n)
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
                newCache[i] = Math.log(i);
            cache = newCache;
        }

        //initialize with the special case: log(0) = NEGATIVE_INFINITY
        private static double[] cache = { Double.NEGATIVE_INFINITY };
    }

    /**
     * Encapsulates the second term of Jacobian log identity for differences up to MAX_TOLERANCE
     */
    private static final class JacobianLogTable {

        // if log(a) - log(b) > MAX_TOLERANCE, b is effectively treated as zero in approximateLogSumLog
        // The following cutoff is the equivalent of 8.0 for base-10 logarithms, which means MAX_TOLERANCE
        // introduces an error of at most one part in 10^8 in sums
        public static final double MAX_TOLERANCE = 8.0 * Math.log(10.0);

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
                    cache[k] = Math.log(1.0 + Math.exp(-((double) k) * TABLE_STEP));
                }
            }
        }

        //  Phred scores Q and Q+1 differ by 0.1 in their corresponding log-10 probabilities, and by
        // 0.1 * log(10) in natural log probabilities.  Setting TABLE_STEP to an exact divisor of this
        // quantity ensures that approximateSumLog in fact caches exact values for integer phred scores
        private static final double TABLE_STEP = (0.1 * Math.log(10.0))/1000;
        private static final double INV_STEP = 1.0 / TABLE_STEP;
        private static double[] cache = null;
    }

    // A fast implementation of the Math.round() method.  This method does not perform
    // under/overflow checking, so this shouldn't be used in the general case (but is fine
    // if one is already make those checks before calling in to the rounding).
    public static int fastRound(final double d) {
        return (d > 0.0) ? (int) (d + 0.5d) : (int) (d - 0.5d);
    }

    public static double approximateLogSumLog(final double[] vals) {
        return approximateLogSumLog(vals, vals.length);
    }

    public static double approximateLogSumLog(final double[] vals, final int endIndex) {

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

    public static double approximateLogSumLog(final double a, final double b, final double c) {
        return approximateLogSumLog(a, approximateLogSumLog(b, c));
    }

    public static double approximateLogSumLog(final double a, final double b) {
        // this code works only when a <= b so we flip them if the order is opposite
        if (a > b) {
            return approximateLogSumLog(b, a);
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
    public static double logBinomialCoefficient(final int n, final int k) {
        if ( n < 0 ) {
            throw new IllegalArgumentException("n: Must have non-negative number of trials");
        }
        if ( k > n || k < 0 ) {
            throw new IllegalArgumentException("k: Must have non-negative number of successes, and no more successes than number of trials");
        }

        return logFactorial(n) - logFactorial(k) - logFactorial(n - k);
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
        return Math.exp(logBinomialProbability(n, k, Math.log(p)));
    }

    /**
     * binomial Probability(int, int, double) with log applied to result
     */
    public static double logBinomialProbability(final int n, final int k, final double logp) {
        if ( logp > 1e-18 )
            throw new IllegalArgumentException("logp: Log-probability must be 0 or less");
        double logOneMinusP = Math.log(1 - Math.exp(logp));
        return logBinomialCoefficient(n, k) + logp * k + logOneMinusP * (n - k);
    }

    public static double logSumLog(final double[] logValues, final int start) {
        return logSumLog(logValues, start, logValues.length);
    }

    public static double logSumLog(final double[] logValues) {
        return logSumLog(logValues, 0);
    }

    public static double logSumLog(final double[] logValues, final int start, final int finish) {
        if (start >= finish) {
            return Double.NEGATIVE_INFINITY;
        }
        final int maxElementIndex = maxElementIndex(logValues, start, finish);
        final double maxValue = logValues[maxElementIndex];
        if(maxValue == Double.NEGATIVE_INFINITY) {
            return maxValue;
        }
        double sum = 1.0;
        for (int i = start; i < finish; i++) {
            double curVal = logValues[i];
            double scaled_val = curVal - maxValue;
            if (i == maxElementIndex || curVal == Double.NEGATIVE_INFINITY) {
                continue;
            } else {
                sum += Math.exp(scaled_val);
            }
        }
        if ( Double.isNaN(sum) || sum == Double.POSITIVE_INFINITY ) {
            throw new IllegalArgumentException("log p: Values must be non-infinite and non-NAN");
        }
        return maxValue + (sum != 1.0 ? Math.log(sum) : 0.0);
    }

    /**
     * normalizes the log-probability array.  ASSUMES THAT ALL ARRAY ENTRIES ARE <= 0 (<= 1 IN REAL-SPACE).
     *
     * @param array the array to be normalized
     * @return a newly allocated array corresponding the normalized values in array
     */
    public static double[] normalizeFromLog(final double[] array) {
        return normalizeFromLog(array, false);
    }

    /**
     * normalizes the log-probability array.  ASSUMES THAT ALL ARRAY ENTRIES ARE <= 0 (<= 1 IN REAL-SPACE).
     *
     * @param array             the array to be normalized
     * @param takeLogOfOutput if true, the output will be transformed back into log units
     * @return a newly allocated array corresponding the normalized values in array, maybe log transformed
     */
    public static double[] normalizeFromLog(final double[] array, final boolean takeLogOfOutput) {
        return normalizeFromLog(array, takeLogOfOutput, false);
    }


    /**
     * See #normalizeFromLog but with the additional option to use an approximation that keeps the calculation always in log-space
     *
     * @param array
     * @param takeLogOfOutput
     * @param keepInLogSpace
     *
     * @return
     */
    public static double[] normalizeFromLog(final double[] array, final boolean takeLogOfOutput, final boolean keepInLogSpace) {
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
            normalized[i] = Math.exp(array[i] - maxValue);

        // normalize
        double sum = 0.0;
        for (int i = 0; i < array.length; i++)
            sum += normalized[i];
        for (int i = 0; i < array.length; i++) {
            double x = normalized[i] / sum;
            if (takeLogOfOutput) {
                x = Math.log(x);
                if ( x < LOG_P_OF_ZERO || Double.isInfinite(x) )
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
     * Checks that the result is a well-formed log probability
     *
     * @param result a supposedly well-formed log probability value.  By default allows
     *               -Infinity values, as log(0.0) == -Infinity.
     * @return true if result is really well formed
     */
    public static boolean goodLogProbability(final double result) {
        return goodLogProbability(result, true);
    }

    /**
     * Checks that the result is a well-formed log probability
     *
     * @param result a supposedly well-formed log probability value
     * @param allowNegativeInfinity should we consider a -Infinity value ok?
     * @return true if result is really well formed
     */
    public static boolean goodLogProbability(final double result, final boolean allowNegativeInfinity) {
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

    public static double logFactorial(final int x) {
        if (x >= LogFactorialCache.size() || x < 0)
            return Gamma.logGamma(x + 1);
        else
            return LogFactorialCache.get(x);
    }

    /**
     * Wrapper class so that the logFactorial array is only calculated if it's used
     */
    private static class LogFactorialCache {

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
                LogCache.ensureCacheContains(CACHE_SIZE);
                cache = new double[CACHE_SIZE];
                cache[0] = 0.0;
                for (int k = 1; k < cache.length; k++)
                    cache[k] = cache[k-1] + LogCache.get(k);
            }
        }

        private static double[] cache = null;
    }

    /**
     * Compute in a numerical correct way the quantity log(1-x)
     *
     * Uses the approximation log(1-x) = log(1/x - 1) + log(x) to avoid very quick underflow
     * in 1-x when x is very small
     *
     * @param x a positive double value between 0.0 and 1.0
     * @return an estimate of log(1-x)
     */
    public static double logOneMinusX(final double x) {
        if ( x == 1.0 )
            return Double.NEGATIVE_INFINITY;
        else if ( x == 0.0 )
            return 0.0;
        else {
            final double d = Math.log(1 / x - 1) + Math.log(x);
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
}
