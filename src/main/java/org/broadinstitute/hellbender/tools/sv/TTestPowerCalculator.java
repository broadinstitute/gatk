package org.broadinstitute.hellbender.tools.sv;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.TDistribution;
import org.broadinstitute.hellbender.utils.Utils;
import org.hipparchus.stat.descriptive.moment.StandardDeviation;

public class TTestPowerCalculator {

    public static double powerCalc(double[] control, double[] treat) {

        // Check if both groups have more than 1 sample
        if (control.length > 1 && treat.length > 1) {

            // Set treatment group to be the smaller
            if (treat.length > control.length) {
                final double[] tmp = control;
                control = treat;
                treat = tmp;
            }

            // Calculate standard deviation of control group
            double controlSd = new StandardDeviation().evaluate(control) * (control.length - 1.0) / control.length;

            // Avoid division by zero
            if (controlSd == 0.0 || Double.isNaN(controlSd)) {
                return Double.NaN;
            }

            // Calculate effect size (Cohen's d)
            double effectSize = 0.5 / controlSd;

            // Calculate power using two-sample t-test with unequal sample sizes
            double power = powerT2nTest(control.length, treat.length, 0.05,
                    "greater", effectSize);

            return power;
        } else {
            return Double.NaN;
        }
    }

    /**
     * Calculate statistical power for two-sample t-test with unequal sample sizes
     * This is equivalent to R's pwr.t2n.test function
     *
     * @param n1 Sample size of group 1
     * @param n2 Sample size of group 2
     * @param sigLevel Significance level (alpha)
     * @param alternative Type of alternative hypothesis
     * @param effectSize Effect size (Cohen's d)
     * @return Statistical power
     */
    public static double powerT2nTest(int n1, int n2, double sigLevel,
                                      String alternative, double effectSize) {
        Utils.validateArg(n1 > 1, "n1 less than 2");
        Utils.validateArg(n2 > 1, "n2 less than 2");
        Utils.validateArg(sigLevel > 0 && sigLevel < 1, "significance level not on (0, 1)");
        Utils.nonNull(alternative);
        Utils.validateArg(effectSize > 0, "effect size must be positive");

        // Calculate degrees of freedom
        double df = n1 + n2 - 2;

        // Calculate pooled standard error
        double se = Math.sqrt((1.0 / n1) + (1.0 / n2));

        // Calculate non-centrality parameter
        double ncp = effectSize / se;

        // Calculate critical value based on alternative hypothesis and significance level
        double criticalValue;
        final TDistribution tDistribution = new TDistribution(df);
        if ("greater".equals(alternative)) {
            criticalValue = tDistribution.inverseCumulativeProbability(1 - sigLevel);
        } else if ("less".equals(alternative)) {
            criticalValue = tDistribution.inverseCumulativeProbability(sigLevel);
        } else { // "two.sided"
            criticalValue = tDistribution.inverseCumulativeProbability(1 - sigLevel / 2);
        }

        // Calculate power using non-central t-distribution
        double power;
        if ("greater".equals(alternative)) {

            power = 1 - nonCentralTCdf(criticalValue, df, ncp);
        } else if ("less".equals(alternative)) {
            power = nonCentralTCdf(criticalValue, df, ncp);
        } else { // "two.sided"
            power = nonCentralTCdf(-criticalValue, df, ncp) +
                    (1 - nonCentralTCdf(criticalValue, df, ncp));
        }

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
    private static double nonCentralTCdf(double t, double df, double ncp) {
        // For small non-centrality parameter, use approximation
        if (Math.abs(ncp) < 0.1) {
            return centralTCdf(t, df);
        }

        // Simple approximation using shifted central t-distribution
        double shiftedT = t - ncp * Math.sqrt(df / (df + t * t));
        return centralTCdf(shiftedT, df);
    }

    /**
     * Central t-distribution cumulative distribution function
     * Approximation using normal distribution for large df
     *
     * @param t t-value
     * @param df Degrees of freedom
     * @return Cumulative probability
     */
    private static double centralTCdf(double t, double df) {
        if (df > 30) {
            return new NormalDistribution().cumulativeProbability(t);
        } else {
            // Simple approximation for small df
            return 0.5 + 0.5 * Math.signum(t) * Math.sqrt(t * t / (t * t + df));
        }
    }
}