package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import java.util.stream.IntStream;

/**
 * Extra MathUtils that should be moved to gatk-public
 * Created by davidben on 1/22/16.
 */
public class GATKProtectedMathUtils {

    /**
     * Computes ln(Sum_i(e^a_i)) trying to avoid underflow issues by using the log-sum-exp trick.
     *
     * <p>
     * This trick consists on shift all the log values by the maximum so that exponent values are
     * much larger (close to 1) before they are summed. Then the result is shifted back down by
     * the same amount in order to obtain the correct value.
     * </p>
     * @return any double value.
     */
    public static double naturalLogSumExp(final double ... values) {
        double max = MathUtils.arrayMax(Utils.nonNull(values));
        double sum = 0.0;
        for (int i = 0; i < values.length; ++i) {
            if (values[i] != Double.NEGATIVE_INFINITY) {
                sum += java.lang.Math.exp(values[i] - max);
            }
        }
        return max + Math.log(sum);
    }

    public static double interquartileRange(final double ... values) {
        final Percentile percentile = new Percentile();
        percentile.setData(values);
        return percentile.evaluate(75.0) - percentile.evaluate(25.0);
    }

    public static double mean(final double ... values) {
        Utils.nonNull(values);
        return MathUtils.mean(values, 0, values.length);
    }

    public static double[] rowMeans(final RealMatrix matrix) {
        Utils.nonNull(matrix);
        return IntStream.range(0, matrix.getRowDimension()).mapToDouble(r -> mean(matrix.getRow(r))).toArray();
    }

    public static double[] rowVariances(final RealMatrix matrix) {
        return IntStream.range(0, matrix.getRowDimension()).mapToDouble(r -> new Variance().evaluate(matrix.getRow(r))).toArray();
    }
}
