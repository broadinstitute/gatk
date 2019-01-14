package org.broadinstitute.hellbender.utils.mcmc;

import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.function.BiFunction;
import java.util.function.Function;


/**
 * Unit tests for {@link MinibatchSliceSampler}.  Modified from tests in {@link SliceSamplerUnitTest}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class MinibatchSliceSamplerUnitTest {
    private static final int RANDOM_SEED = 1;
    private static final RandomGenerator rng =
            RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));

    private static double relativeError(final double x, final double xTrue) {
        return Math.abs((x - xTrue) / xTrue);
    }

    private static final int NUM_DATA_POINTS = 1000;
    private static final Function<Double, Double> UNIFORM_LOG_PRIOR = x -> 0.;
    private static final int MINIBATCH_SIZE = 100;
    private static final double APPROX_THRESHOLD = 1E-3;

    /**
     * Tests slice sampling of a normal posterior = uniform prior x normal likelihood from 1000 data points.
     * Checks that input mean and standard deviation are recovered from the posterior of the mean parameter
     * by 500 burn-in samples + 1000 samples to a relative error of 1% and 5%, respectively.
     */
    @Test
    public void testSliceSamplingOfNormalPosterior() {
        rng.setSeed(RANDOM_SEED);

        final double mean = 5.;
        final double standardDeviation = 0.75;
        final NormalDistribution normalDistribution = new NormalDistribution(rng, mean, standardDeviation);
        final BiFunction<Double, Double, Double> normalLogLikelihood =
                (d, x) -> new NormalDistribution(null, x, standardDeviation).logDensity(d);
        final List<Double> data = Doubles.asList(normalDistribution.sample(NUM_DATA_POINTS));

        final double xInitial = 1.;
        final double xMin = Double.NEGATIVE_INFINITY;
        final double xMax = Double.POSITIVE_INFINITY;
        final double width = 0.5;
        final int numBurnInSamples = 500;
        final int numSamples = 1500;
        final MinibatchSliceSampler<Double> normalSampler = new MinibatchSliceSampler<>(
                rng, data, UNIFORM_LOG_PRIOR, normalLogLikelihood,
                xMin, xMax, width, MINIBATCH_SIZE, APPROX_THRESHOLD);
        final double[] samples = Doubles.toArray(normalSampler.sample(xInitial, numSamples).subList(numBurnInSamples, numSamples));

        final double sampleMean = new Mean().evaluate(samples);
        final double sampleStandardDeviation = new StandardDeviation().evaluate(samples);
        Assert.assertEquals(relativeError(sampleMean, mean), 0., 0.01);
        Assert.assertEquals(relativeError(sampleStandardDeviation, standardDeviation / Math.sqrt(NUM_DATA_POINTS)), 0., 0.05);
    }

    /**
     * Tests slice sampling of a beta posterior = uniform prior x binomial likelihood from 1000 data points.
     * Checks that input mean and standard deviation are recovered from the posterior of the mean parameter
     * by 500 burn-in samples + 1000 samples to a relative error of 1% and 5%, respectively.
     */
    @Test
    public void testSliceSamplingOfBetaPosterior() {
        rng.setSeed(RANDOM_SEED);

        final int numTrials = 100;
        final double p = 0.9;
        final BinomialDistribution binomialDistribution = new BinomialDistribution(rng, numTrials, p);
        final BiFunction<Integer, Double, Double> binomialLogLikelihood =
                (d, x) -> new BinomialDistribution(null, numTrials, x).logProbability(d);
        final List<Integer> data = Ints.asList(binomialDistribution.sample(NUM_DATA_POINTS));
        final double numSuccessesTotal = data.stream().mapToDouble(x -> x).sum();
        final double alpha = 1. + numSuccessesTotal;
        final double beta = 1. + (numTrials * NUM_DATA_POINTS - numSuccessesTotal);
        final double variance = new BetaDistribution(null, alpha, beta).getNumericalVariance();

        final double xInitial = 0.5;
        final double xMin = 0.;
        final double xMax = 1.;
        final double width = 0.1;
        final int numBurnInSamples = 500;
        final int numSamples = 1500;
        final MinibatchSliceSampler<Integer> betaSampler = new MinibatchSliceSampler<>(
                rng, data, UNIFORM_LOG_PRIOR, binomialLogLikelihood,
                xMin, xMax, width, MINIBATCH_SIZE, APPROX_THRESHOLD);
        final double[] samples = Doubles.toArray(betaSampler.sample(xInitial, numSamples).subList(numBurnInSamples, numSamples));

        final double sampleMean = new Mean().evaluate(samples);
        final double sampleStandardDeviation = new StandardDeviation().evaluate(samples);
        Assert.assertEquals(relativeError(sampleMean, p), 0., 0.01);
        Assert.assertEquals(relativeError(sampleStandardDeviation, Math.sqrt(variance)), 0., 0.05);
    }

    /**
     * Tests slice sampling of a uniform prior with no data points.
     * Checks that mean and standard deviation are recovered from the samples
     * by 500 burn-in samples + 1000 samples to a relative error of 1% and 5%, respectively.
     */
    @Test
    public void testSliceSamplingOfUniformPriorWithNoData() {
        rng.setSeed(RANDOM_SEED);

        final double mean = 0.5;
        final double standardDeviation = 1. / Math.sqrt(12.);
        final List<Double> data = Collections.emptyList();

        final double xInitial = 0.5;
        final double xMin = 0.;
        final double xMax = 1.;
        final double width = 0.1;
        final int numBurnInSamples = 500;
        final int numSamples = 1500;
        final MinibatchSliceSampler<Double> uniformSampler = new MinibatchSliceSampler<>(
                rng, data, UNIFORM_LOG_PRIOR, (d, x) -> 0.,
                xMin, xMax, width, MINIBATCH_SIZE, APPROX_THRESHOLD);
        final double[] samples = Doubles.toArray(uniformSampler.sample(xInitial, numSamples).subList(numBurnInSamples, numSamples));

        final double sampleMean = new Mean().evaluate(samples);
        final double sampleStandardDeviation = new StandardDeviation().evaluate(samples);
        Assert.assertEquals(relativeError(sampleMean, mean), 0., 0.01);
        Assert.assertEquals(relativeError(sampleStandardDeviation, standardDeviation), 0., 0.05);
    }
}