package org.broadinstitute.hellbender.utils.mcmc.univariatesamplers;

import com.google.common.primitives.Doubles;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Random;
import java.util.function.Function;


/**
 * Unit tests for {@link SliceSampler}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SliceSamplerUnitTest {
    private static final int RANDOM_SEED = 42;
    private static final RandomGenerator rng =
            RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));

    private static double relativeError(final double x, final double xTrue) {
        return Math.abs((x - xTrue) / xTrue);
    }

    /**
     * Test slice sampling of a normal distribution.  Checks that input mean and standard deviation are recovered
     * by 10000 samples to a relative error of 0.5% and 2%, respectively.
     */
    @Test
    public void testSliceSamplingOfNormalDistribution() {
        rng.setSeed(RANDOM_SEED);

        final double mean = 5.;
        final double standardDeviation = 0.75;
        final NormalDistribution normalDistribution = new NormalDistribution(mean, standardDeviation);
        final Function<Double, Double> normalLogPDF = normalDistribution::logDensity;

        final double xInitial = 1.;
        final double xMin = Double.NEGATIVE_INFINITY;
        final double xMax = Double.POSITIVE_INFINITY;
        final double width = 0.5;
        final int numSamples = 10000;
        final SliceSampler normalSampler = new SliceSampler(rng, normalLogPDF, xMin, xMax, width);
        final double[] samples = Doubles.toArray(normalSampler.sample(xInitial, numSamples));

        final double sampleMean = new Mean().evaluate(samples);
        final double sampleStandardDeviation = new StandardDeviation().evaluate(samples);
        Assert.assertEquals(relativeError(sampleMean, mean), 0., 0.005);
        Assert.assertEquals(relativeError(sampleStandardDeviation, standardDeviation), 0., 0.02);
    }

    /**
     * Test slice sampling of a monotonic beta distribution as an example of sampling of a bounded random variable.
     * Checks that input mean and variance are recovered by 10000 samples to a relative error of 0.5% and 2%,
     * respectively.
     */
    @Test
    public void testSliceSamplingOfMonotonicBetaDistribution() {
        rng.setSeed(RANDOM_SEED);

        final double alpha = 10.;
        final double beta = 1.;
        final BetaDistribution betaDistribution = new BetaDistribution(alpha, beta);
        final Function<Double, Double> betaLogPDF = betaDistribution::logDensity;

        final double xInitial = 0.5;
        final double xMin = 0.;
        final double xMax = 1.;
        final double width = 0.1;
        final int numSamples = 10000;
        final SliceSampler betaSampler = new SliceSampler(rng, betaLogPDF, xMin, xMax, width);
        final double[] samples = Doubles.toArray(betaSampler.sample(xInitial, numSamples));

        final double mean = betaDistribution.getNumericalMean();
        final double variance = betaDistribution.getNumericalVariance();
        final double sampleMean = new Mean().evaluate(samples);
        final double sampleVariance = new Variance().evaluate(samples);
        Assert.assertEquals(relativeError(sampleMean, mean), 0., 0.005);
        Assert.assertEquals(relativeError(sampleVariance, variance), 0., 0.02);
    }

    /**
     * Test slice sampling of a peaked beta distribution as an example of sampling of a bounded random variable.
     * Checks that input mean and variance are recovered by 10000 samples to a relative error of 0.5% and 2%,
     * respectively.
     */
    @Test
    public void testSliceSamplingOfPeakedBetaDistribution() {
        rng.setSeed(RANDOM_SEED);

        final double alpha = 10.;
        final double beta = 4.;
        final BetaDistribution betaDistribution = new BetaDistribution(alpha, beta);
        final Function<Double, Double> betaLogPDF = betaDistribution::logDensity;

        final double xInitial = 0.5;
        final double xMin = 0.;
        final double xMax = 1.;
        final double width = 0.1;
        final int numSamples = 10000;
        final SliceSampler betaSampler = new SliceSampler(rng, betaLogPDF, xMin, xMax, width);
        final double[] samples = Doubles.toArray(betaSampler.sample(xInitial, numSamples));

        final double mean = betaDistribution.getNumericalMean();
        final double variance = betaDistribution.getNumericalVariance();
        final double sampleMean = new Mean().evaluate(samples);
        final double sampleVariance = new Variance().evaluate(samples);
        Assert.assertEquals(relativeError(sampleMean, mean), 0., 0.005);
        Assert.assertEquals(relativeError(sampleVariance, variance), 0., 0.02);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInitialPointOutOfRange() {
        rng.setSeed(RANDOM_SEED);

        final double mean = 5.;
        final double standardDeviation = 0.75;
        final NormalDistribution normalDistribution = new NormalDistribution(mean, standardDeviation);
        final Function<Double, Double> normalLogPDF = normalDistribution::logDensity;

        final double xInitial = -10.;
        final double xMin = 0.;
        final double xMax = 1.;
        final double width = 0.5;
        final SliceSampler normalSampler = new SliceSampler(rng, normalLogPDF, xMin, xMax, width);
        normalSampler.sample(xInitial);
    }
}