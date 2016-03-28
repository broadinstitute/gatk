package org.broadinstitute.hellbender.utils.mcmc;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

/**
 * Tests for {@link PosteriorSummaryUtils}.  Tests that posterior mode and highest-posterior-density credible interval
 * are recovered for samples drawn from various posterior distributions.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class PosteriorSummaryUtilsUnitTest extends BaseTest {
    private static final int RANDOM_SEED = 42;
    private static final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));

    private static final List<Double> normalSamples = toList(new NormalDistribution(rng, 10., 1).sample(100));
    private static final List<Double> normalSamplesSmall = toList(new NormalDistribution(rng, 10., 1).sample(10));
    private static final List<Double> normalSamplesScaled = toList(new NormalDistribution(rng, 10., 0.01).sample(100));

    private static final List<Double> betaSamples = toList(new BetaDistribution(rng, 10, 4).sample(1000));
    private static final List<Double> betaSamplesEdgeMode = toList(new BetaDistribution(rng, 10, 1).sample(1000));
    private static final List<Double> betaSamplesMAF = Arrays.stream(ArrayUtils.addAll(
            new BetaDistribution(rng, 10, 6).sample(1000),
            new BetaDistribution(rng, 6, 10).sample(1000)))
            .boxed().filter(d -> d <= 0.5).collect(Collectors.toList());    //mimics ACNV minor-allele-fraction posterior

    private static final List<Double> identicalSamples = Collections.nCopies(1000, 1.);
    private static final List<Double> withNaNSamples = toList(new double[]{Double.NEGATIVE_INFINITY, 0., 1.});

    @DataProvider(name = "dataKernelDensityEstimation")
    public Object[][] dataHPD() {
        return new Object[][]{
                //samples, alpha, relative error tolerance, expected (calculated using Mathematica for Beta tests)
                {normalSamples, 0.32, 0.05, new PosteriorSummary(10., 9., 11.)},
                {normalSamples, 0.05, 0.05, new PosteriorSummary(10., 8., 12.)},
                {normalSamplesSmall, 0.32, 0.15, new PosteriorSummary(10., 9., 11.)},
                {normalSamplesSmall, 0.05, 0.15, new PosteriorSummary(10., 8., 12.)},
                {normalSamplesScaled, 0.32, 0.05, new PosteriorSummary(10., 9.99, 10.01)},
                {normalSamplesScaled, 0.05, 0.05, new PosteriorSummary(10., 9.98, 10.02)},
                {betaSamples, 0.32, 0.05, new PosteriorSummary(0.75, 0.621, 0.854)},
                {betaSamples, 0.05, 0.05, new PosteriorSummary(0.75, 0.487, 0.925)},
                {betaSamplesEdgeMode, 0.32, 0.05, new PosteriorSummary(1.0, 0.8923, 1.)},
                {betaSamplesEdgeMode, 0.05, 0.05, new PosteriorSummary(1.0, 0.7411, 1.)},
                {betaSamplesMAF, 0.32, 0.05, new PosteriorSummary(0.415, 0.312, 0.5)},
                {betaSamplesMAF, 0.05, 0.05, new PosteriorSummary(0.415, 0.191, 0.5)},
                {identicalSamples, 0.32, 0.01, new PosteriorSummary(1., 1., 1.)},
                {withNaNSamples, 0.32, 0.01, new PosteriorSummary(Double.NaN, Double.NaN, Double.NaN)},
        };
    }

    @Test(dataProvider = "dataKernelDensityEstimation")
    public void testCalculateHighestPosteriorDensitySummary(final List<Double> samples,
                                                            final double credibleIntervalAlpha,
                                                            final double relativeError,
                                                            final PosteriorSummary expected) {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final PosteriorSummary result =
                PosteriorSummaryUtils.calculateHighestPosteriorDensitySummary(samples, credibleIntervalAlpha, ctx);
        assertEquals(result, expected, relativeError);
    }

    @Test(dataProvider = "dataKernelDensityEstimation")
    public void testCalculatePosteriorMode(final List<Double> samples,
                                           final double credibleIntervalAlpha,
                                           final double relativeError,
                                           final PosteriorSummary expected) {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final double result = PosteriorSummaryUtils.calculatePosteriorMode(samples, ctx);
        Assert.assertTrue(withinRelativeError(result, expected.getCenter(), relativeError));

    }

    private static boolean withinRelativeError(final double x, final double xTrue, final double relativeError) {
        if (Double.isNaN(xTrue)) {
            return Double.isNaN(x);
        }
        return Math.abs(x - xTrue) < Math.abs(relativeError * xTrue);
    }

    private static List<Double> toList(double[] array) {
        return Arrays.stream(array).boxed().collect(Collectors.toList());   }

    private static void assertEquals(final PosteriorSummary result, final PosteriorSummary expected, final double relativeError) {
        Assert.assertTrue(withinRelativeError(result.getCenter(), expected.getCenter(), relativeError));
        Assert.assertTrue(withinRelativeError(result.getLower(), expected.getLower(), relativeError));
        Assert.assertTrue(withinRelativeError(result.getUpper(), expected.getUpper(), relativeError));
    }
}