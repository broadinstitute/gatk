package org.broadinstitute.hellbender.tools.sv;

import org.apache.commons.math3.distribution.TDistribution;
import org.broadinstitute.hellbender.utils.Utils;
import org.hipparchus.stat.descriptive.moment.StandardDeviation;

class TTestPowerCalculator {

    public static double power(double[] control, double[] treat, final double sigLevel) {
        Utils.nonNull(control);
        Utils.nonNull(treat);
        Utils.validateArg(sigLevel > 0 && sigLevel < 1, "significance level not on (0, 1)");

        // Check if both groups have more than 1 sample

            // Set treatment group to be the smaller
            if (treat.length > control.length) {
                final double[] tmp = control;
                control = treat;
                treat = tmp;
            }

            // Calculate standard deviation of control group
            final double controlSd = new StandardDeviation().evaluate(control) * (control.length - 1.0) / control.length;

            // Avoid division by zero
            if (controlSd == 0.0 || Double.isNaN(controlSd)) {
                return Double.NaN;
            }

            // Calculate effect size (Cohen's d)
            final double effectSize = 0.5 / controlSd;

            // Calculate power using two-sample t-test with unequal sample sizes
            return powerT2nTest(control.length, treat.length, sigLevel, effectSize);
    }

    /**
     * Calculate statistical power for two-sample t-test with unequal sample sizes
     * This is equivalent to R's pwr.t2n.test function, assuming "greater" alternative hypothesis
     *
     * @param n1 Sample size of group 1
     * @param n2 Sample size of group 2
     * @param sigLevel Significance level (alpha)
     * @param effectSize Effect size (Cohen's d)
     * @return Statistical power
     */
    private static double powerT2nTest(int n1, int n2, double sigLevel, double effectSize) {
        Utils.validateArg(sigLevel > 0 && sigLevel < 1, "significance level not on (0, 1)");
        Utils.validateArg(effectSize > 0, "effect size must be positive");

        // Not enough samples
        if (n1 < 2 || n2 < 2) {
            return Double.NaN;
        }

        // Calculate degrees of freedom
        double df = n1 + n2 - 2;

        // Calculate pooled standard error
        double se = Math.sqrt((1.0 / n1) + (1.0 / n2));

        // Calculate non-centrality parameter
        double ncp = effectSize / se;

        // Calculate critical value based on alternative hypothesis and significance level
        // Assume "greater than" alternative hypothesis
        final TDistribution tDistribution = new TDistribution(df);
        final double criticalValue = tDistribution.inverseCumulativeProbability(1 - sigLevel);

        // Calculate power using non-central t-distribution
        final double power = 1 - approximateNonCentralTCdf(criticalValue, df, ncp);

        return Math.max(0.0, Math.min(1.0, power)); // Ensure power is between 0 and 1
    }

    /**
     * Non-central t-distribution cumulative distribution function
     * Approximation for moderate non-centrality parameters
     *
     * @param t t-value
     * @param df Degrees of freedom
     * @param ncp Non-centrality parameter
     * @return Cumulative probability
     */
    private static double approximateNonCentralTCdf(double t, double df, double ncp) {
        // Simple approximation using shifted central t-distribution
        double shiftedT = t - ncp * Math.sqrt(df / (df + t * t));
        return new TDistribution(df).cumulativeProbability(shiftedT);
    }
}