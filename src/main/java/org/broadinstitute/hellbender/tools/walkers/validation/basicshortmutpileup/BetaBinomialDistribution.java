package org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup;

import org.apache.commons.math3.distribution.AbstractIntegerDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import static java.lang.Math.max;

/**
 * Beta-binomial using the Apache Math3 Framework.
 *
 * Like a binomial distribution, but the binomial p parameter is itself drawn from a beta distribution.
 *
 * https://en.wikipedia.org/wiki/Beta-binomial_distribution
 */
public class BetaBinomialDistribution extends AbstractIntegerDistribution {

    private final static long serialVersionUID = 2325;
    /**
     * Beta distribution shape parameter
     */
    private double alpha;

    /**
     * Beta distribution shape parameter
     */
    private double beta;

    /**
     * Number of trials
     */
    private int n;

    /**
     * @param rng a random number generator.  {@code null} is okay
     * @param alpha shape parameter for underlying beta distribution.  Must be greater than zero.
     * @param beta shape parameter for underlying beta distribution.  Must be greater than zero.
     * @param n number of trials.  Must be Zero or positive.
     */
    public BetaBinomialDistribution(RandomGenerator rng, double alpha, double beta, int n) {
        super(rng);

        ParamUtils.isPositive(alpha, "alpha must be greater than zero.");
        ParamUtils.isPositive(beta, "beta must be greater than zero.");
        ParamUtils.isPositiveOrZero(n, "number of trials must be greater than (or equal to) zero.");

        this.alpha = alpha;
        this.beta = beta;
        this.n = n;
    }

    /**
     * @param k number of successes.  Must be positive or zero.
     * @return the value of the pdf at k
     */
    @Override
    public double probability(int k) {
        return Math.exp(logProbability(k));
    }

    /**
     * @param k number of successes.  Must be positive or zero.
     * @return the value of the pdf at k
     */
    @Override
    public double logProbability(int k) {
        ParamUtils.isPositiveOrZero(k, "Number of successes must be greater than or equal to zero.");

        // nchoosek * beta(k+alpha, n-k+beta)/ beta(alpha, beta)
        // binomialcoefficient is "n choose k"
        return k > n ? Double.NEGATIVE_INFINITY :
                CombinatoricsUtils.binomialCoefficientLog(n,k) + Beta.logBeta(k+alpha, n - k + beta) - Beta.logBeta(alpha, beta);
    }

    /**
     * @param k number of successes.  Must be positive or zero.
     * @return the value of the cdf at k
     */
    @Override
    public double cumulativeProbability(int k) {

        ParamUtils.isPositiveOrZero(k, "Number of successes must be greater than or equal to zero.");

        // s l o w, but possibly not as slow as calculating a generalized hypergeometric function.
        double result = 0.0;
        for (int i = 0; i <= k; i ++) {
            result += probability(i);
        }
        return result;
    }

    @Override
    public double getNumericalMean() {
        return (n*alpha)/(alpha+beta);
    }

    @Override
    public double getNumericalVariance() {
        return (n*alpha*beta)*(alpha+beta+n) / ((alpha+beta)*(alpha+beta)*(alpha+beta+1));
    }

    @Override
    public int getSupportLowerBound() {
        return 0;
    }

    @Override
    public int getSupportUpperBound() {
        return n;
    }

    @Override
    public boolean isSupportConnected() {
        return true;
    }
}
