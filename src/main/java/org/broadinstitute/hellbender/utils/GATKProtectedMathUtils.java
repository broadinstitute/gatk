package org.broadinstitute.hellbender.utils;

/**
 * Extra MathUtils that should be moved to gatk-public
 * Created by davidben on 1/22/16.
 */
public class GATKProtectedMathUtils {

    /**
     * Compute log(e^a + e^b + e^c. . .) = log(e^M [e^(a-M) + e^(b-M) + e^(c-M). . .]) where M = max(a,b,c. . .)
     *
     * @return
     */

    public static double naturalLogSumNaturalLog(double[] values) {
        double max = MathUtils.arrayMax(values);
        double sum = 0.0;
        for (int i = 0; i < values.length; ++i)
            if (values[i] != Double.NEGATIVE_INFINITY) {
                sum += java.lang.Math.exp(values[i] - max);
            }
        return max + Math.log(sum);
    }
}
