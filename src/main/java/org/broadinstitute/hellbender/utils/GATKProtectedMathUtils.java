package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import java.util.Collection;
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
     * Calculates the variances per column from a matrix.
     * @param matrix the input matrix.
     * @return never {@code null}, an array with as many positions as columns in {@code matrix}.
     * @throws IllegalArgumentException if {@code matrix} is {@code null}.
     */
    public static double[] columnVariances(final RealMatrix matrix) {
        Utils.nonNull(matrix);
        final Variance varianceEvaluator = new Variance();
        return IntStream.range(0, matrix.getColumnDimension())
                .mapToDouble(c -> varianceEvaluator.evaluate(matrix.getColumn(c)))
                .toArray();
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
        if (values.contains(null)) {
            throw new IllegalArgumentException("input values must not contain a null");
        }
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
}
