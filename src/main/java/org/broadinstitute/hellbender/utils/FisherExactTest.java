package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.distribution.AbstractIntegerDistribution;
import org.apache.commons.math3.distribution.HypergeometricDistribution;

import java.util.Arrays;
import java.util.function.DoubleUnaryOperator;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import static java.lang.Math.*;
import static org.apache.commons.math3.util.MathArrays.ebeAdd;
import static org.broadinstitute.hellbender.utils.MathUtils.arrayMax;
import static org.broadinstitute.hellbender.utils.MathUtils.promote;
import static org.broadinstitute.hellbender.utils.MathUtils.sum;

/**
 * Implements the Fisher's exact test for 2x2 tables
 * assuming the null hypothesis of odd ratio of 1.
 */
public final class FisherExactTest {

    /**
     * Computes the 2-sided pvalue of the Fisher's exact test.
     */
    public static double twoSidedPValue(final int[][] normalizedTable) {
        Utils.nonNull(normalizedTable);
        Utils.validateArg(normalizedTable.length == 2, "input must be 2x2 " + Arrays.deepToString(normalizedTable));
        Utils.validateArg(normalizedTable[0] != null && normalizedTable[0].length == 2, "input must be 2x2 " + Arrays.deepToString(normalizedTable));
        Utils.validateArg(normalizedTable[1] != null && normalizedTable[1].length == 2, "input must be 2x2 " + Arrays.deepToString(normalizedTable));

        //Note: this implementation follows the one in R base package
        final int[][] x= normalizedTable;
        final int m = addExact(x[0][0], x[0][1]);
        final int n = addExact(x[1][0], x[1][1]);

        if (m+n == 0){     //special case, all entries zero
            return 1.0;
        }

        final int k = addExact(x[0][0], x[1][0]);
        final int lo = max(0, k - n);
        final int hi = min(k, m);
        final int[] support = range(lo, hi);
        final AbstractIntegerDistribution dist = new HypergeometricDistribution(null, m+n, m, k);
        final double[] logdc = IntStream.of(support).mapToDouble(i -> dist.logProbability(i)).toArray();
        final double oddsRatioAtNull = 1.0;
        final double[] ds = dnHyper(oddsRatioAtNull, logdc, support);
        final double relErr = 1 + pow(10, -7);

        final double threshold= ds[x[0][0] - lo] * relErr;
        final double pValue = DoubleStream.of(ds).filter(d -> d <= threshold).sum();

        // min is necessary as numerical precision can result in pValue being slightly greater than 1.0
        return min(pValue, 1.0);
    }

    private static double[] dnHyper(final double ncp, final double[] logdc, final int[] support){
        final double[] d1 = ebeAdd(logdc, apply(promote(support), d -> d * log(ncp)));
        final double maxD1 = arrayMax(d1);

        final double[] d2 = apply(apply(d1, d -> d - maxD1), d -> exp(d));
        final double sumD2 = sum(d2);

        return apply(d2,  d -> d/sumD2);
    }

    private static double[] apply(final double[] arr, final DoubleUnaryOperator op){
        return DoubleStream.of(arr).map(op).toArray();
    }

    private static int[] range(final int from, final int to){
        final int[] support = new int[to-from+1];
        for (int i = 0; i < support.length; i++){
            support[i] = from + i;
        }
        return support;
    }
}
