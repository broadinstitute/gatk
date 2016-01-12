package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.rank.Median;

import java.util.stream.IntStream;

/**
 * Static class for implementing some mtarix operations not in Apache, Spark, etc
 */
public class MatrixUtils {

    private MatrixUtils() {}

    /**
     * Return an array containing the median for each column in the given matrix.
     * @param m Not {@code null}.  Size MxN, where neither dimension is zero.
     * @return array of size N.  Never {@code null}
     */
    public static double[] getColumnMedians(final RealMatrix m) {
        Utils.nonNull(m, "Cannot calculate medians on a null matrix.");
        final Median medianCalculator = new Median();
        return IntStream.range(0, m.getColumnDimension()).boxed()
                .mapToDouble(i -> medianCalculator.evaluate(m.getColumn(i))).toArray();
    }

    /**
     * Return an array containing the median for each row in the given matrix.
     * @param m Not {@code null}.  Size MxN.
     * @return array of size M.  Never {@code null}
     */
    public static double[] getRowMedians(final RealMatrix m) {
        Utils.nonNull(m, "Cannot calculate medians on a null matrix.");
        final Median medianCalculator = new Median();
        return IntStream.range(0, m.getRowDimension()).boxed()
                .mapToDouble(i -> medianCalculator.evaluate(m.getRow(i))).toArray();
    }

    /**
     * Return an array containing the variance for each row in the given matrix.
     * @param m Not {@code null}.  Size MxN.
     * @return array of size M.  Never {@code null}
     */
    public static double[] getRowVariances(final RealMatrix m) {
        Utils.nonNull(m, "Cannot calculate medians on a null matrix.");
        final StandardDeviation std = new StandardDeviation();
        return IntStream.range(0, m.getRowDimension()).boxed()
                .mapToDouble(i -> std.evaluate(m.getColumn(i))).toArray();
    }

}
