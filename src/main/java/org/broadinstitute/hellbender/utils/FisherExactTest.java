package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.distribution.AbstractIntegerDistribution;
import org.apache.commons.math3.distribution.HypergeometricDistribution;

import java.util.Arrays;
import java.util.stream.DoubleStream;

/**
 * Implements the Fisher's exact test for 2x2 tables
 * assuming the null hypothesis of odd ratio of 1.
 */
public final class FisherExactTest {
    private static final double REL_ERR = 1 - 10e-7;

    /**
     * Computes the 2-sided pvalue of the Fisher's exact test on a normalized table that ensures that the sum of
     * all four entries is less than 2 * 200.
     */
    public static double twoSidedPValue(final int[][] normalizedTable) {
        Utils.nonNull(normalizedTable);
        Utils.validateArg(normalizedTable.length == 2, () -> "input must be 2x2 " + Arrays.deepToString(normalizedTable));
        Utils.validateArg(normalizedTable[0] != null && normalizedTable[0].length == 2, () -> "input must be 2x2 " + Arrays.deepToString(normalizedTable));
        Utils.validateArg(normalizedTable[1] != null && normalizedTable[1].length == 2, () -> "input must be 2x2 " + Arrays.deepToString(normalizedTable));

        //Note: this implementation follows the one in R base package
        final int[][] x= normalizedTable;
        final int m = x[0][0] + x[0][1];
        final int n = x[1][0] + x[1][1];
        final int k = x[0][0] + x[1][0];
        final int lo = Math.max(0, k - n);
        final int hi = Math.min(k, m);
        final IndexRange support = new IndexRange(lo, hi + 1);

        if (support.size() <= 1){     //special case, support has only one value
            return 1.0;
        }

        final AbstractIntegerDistribution dist = new HypergeometricDistribution(null, m+n, m, k);
        final double[] logds = support.mapToDouble(dist::logProbability);
        final double threshold = logds[x[0][0] - lo] * REL_ERR;
        final double[] log10ds = DoubleStream.of(logds).filter(d -> d <= threshold).map(MathUtils::logToLog10).toArray();
        final double pValue = MathUtils.sumLog10(log10ds);

        // min is necessary as numerical precision can result in pValue being slightly greater than 1.0
        return Math.min(pValue, 1.0);
    }
}
