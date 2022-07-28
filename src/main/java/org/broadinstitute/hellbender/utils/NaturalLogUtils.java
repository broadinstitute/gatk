package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.MathArrays;

import java.util.Collection;
import java.util.Collections;

public class NaturalLogUtils {
    public static final double LOG_ONE_HALF = FastMath.log(0.5);
    public static final double LOG_ONE_THIRD = FastMath.log(1.0/3);
    private static final double LOG1MEXP_THRESHOLD = Math.log(0.5);
    private static final double PHRED_TO_LOG_ERROR_PROB_FACTOR = -Math.log(10)/10;

    private static final double[] qualToLogProbCache = new IndexRange(0, QualityUtils.MAX_QUAL + 1)
            .mapToDouble(n -> log1mexp(qualToLogErrorProb((double)n)));

    private NaturalLogUtils() { }

    /**
     * normalizes the log-probability array in-place.
     *
     * @param array the array to be normalized
     * @return the normalized-in-place array
     */
    public static double[] normalizeFromLogToLinearSpace(final double[] array) {
        return normalizeLog(Utils.nonNull(array), false, true);
    }

    /**
     * normalizes the log-probability array in-place.
     *
     * @param array             the array to be normalized
     * @return the normalized-in-place array, maybe log transformed
     */
    public static double[] normalizeLog(final double[] array) {
        return normalizeLog(Utils.nonNull(array), true, true);
    }


    /**
     * See #normalizeFromLog but with the additional option to use an approximation that keeps the calculation always in log-space
     *
     * @param array
     * @param takeLogOfOutput
     * @param inPlace           if true, modify the input array in-place
     *
     * @return
     */
    public static double[] normalizeLog(final double[] array, final boolean takeLogOfOutput, final boolean inPlace) {
        final double logSum = logSumExp(Utils.nonNull(array));
        final double[] result = inPlace ? MathUtils.applyToArrayInPlace(array, x -> x - logSum) : MathUtils.applyToArray(array, x -> x - logSum);
        return takeLogOfOutput ? result : MathUtils.applyToArrayInPlace(result, Math::exp);
    }

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
    public static double logSumExp(final double... logValues) {
        Utils.nonNull(logValues);
        final int maxElementIndex = MathUtils.maxElementIndex(logValues);
        final double maxValue = logValues[maxElementIndex];
        if(maxValue == Double.NEGATIVE_INFINITY) {
            return maxValue;
        }
        double sum = 1.0;
        for (int i = 0; i < logValues.length; i++) {
            final double curVal = logValues[i];
            if (i == maxElementIndex || curVal == Double.NEGATIVE_INFINITY) {
                continue;
            } else {
                final double scaled_val = curVal - maxValue;
                sum += Math.exp(scaled_val);
            }
        }
        if ( Double.isNaN(sum) || sum == Double.POSITIVE_INFINITY ) {
            throw new IllegalArgumentException("logValues must be non-infinite and non-NAN");
        }
        return maxValue + (sum != 1.0 ? Math.log(sum) : 0.0);
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

    public static double[] posteriors(double[] logPriors, double[] logLikelihoods) {
        return normalizeFromLogToLinearSpace(MathArrays.ebeAdd(logPriors, logLikelihoods));
    }

    /**
     * Convert a phred-scaled quality score to its log10 probability of being wrong (Q30 => log(0.001))
     *
     * WARNING -- because this function takes a byte for maxQual, you must be careful in converting
     * integers to byte.  The appropriate way to do this is ((byte)(myInt & 0xFF))
     *
     * @param qual a phred-scaled quality score encoded as a byte
     * @return a probability (0.0-1.0)
     */
    public static double qualToLogErrorProb(final byte qual){
        return qualToLogErrorProb((double)(qual & 0xFF));
    }

    /**
     * Convert a phred-scaled quality score to its log10 probability of being wrong (Q30 => log10(0.001))
     *
     * This is the Phred-style conversion, *not* the Illumina-style conversion.
     *
     * The calculation is extremely efficient
     *
     * @param qual a phred-scaled quality score encoded as a double
     * @return log of probability (0.0-1.0)
     */
    public static double qualToLogErrorProb(final double qual) {
        Utils.validateArg( qual >= 0.0, () -> "qual must be >= 0.0 but got " + qual);
        return qual * PHRED_TO_LOG_ERROR_PROB_FACTOR;
    }

    public static double qualToLogProb(final byte qual) {
        return qualToLogProbCache[(int) qual & 0xff]; // Map: 127 -> 127; -128 -> 128; -1 -> 255; etc.
    }

    public static double logSumLog(final double a, final double b) {
        return a > b ? a + FastMath.log(1 + FastMath.exp(b - a)) : b + FastMath.log(1 + FastMath.exp(a - b));
    }
}
