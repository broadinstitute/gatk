package org.broadinstitute.hellbender.utils.mcmc.univariatesamplers;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.function.Function;

/**
 * Created by davidben on 11/24/15.
 */
public final class AdaptiveMetropolisSamplerUnitTest {
    private static final double INITIAL_STEP_SIZE = 0.1;
    private static final double INITIAL_GAUSSIAN_GUESS = 0;
    private static final double INITIAL_BETA_GUESS = 0.5;
    private static final int NUM_BURN_IN_STEPS = 1000;
    private static final int NUM_SAMPLES = 5000;

    private static final int RANDOM_SEED = 13;

    @Test
    public void testBeta() {
        final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));
        for (final double a : Arrays.asList(10, 20, 30)) {
            for (final double b : Arrays.asList(10, 20, 30)) {

                final double theoreticalMean = a / (a + b);
                final double theoreticalVariance = a*b/((a+b)*(a+b)*(a+b+1));

                //Note: this is the theoretical standard deviation of the sample mean given uncorrelated
                //samples.  The sample mean will have a greater variance here because samples are correlated.
                final double standardDeviationOfMean = Math.sqrt(theoreticalVariance / NUM_SAMPLES);

                final Function<Double, Double> logPDF = x -> (a - 1) * Math.log(x) + (b - 1) * Math.log(1 - x);
                final AdaptiveMetropolisSampler sampler = new AdaptiveMetropolisSampler(INITIAL_BETA_GUESS, INITIAL_STEP_SIZE, 0, 1);
                final List<Double> samples = sampler.sample(rng, logPDF, NUM_SAMPLES, NUM_BURN_IN_STEPS);

                final double sampleMean = samples.stream().mapToDouble(x -> x).average().getAsDouble();
                final double sampleMeanSquare = samples.stream().mapToDouble(x -> x*x).average().getAsDouble();
                final double sampleVariance = (sampleMeanSquare - sampleMean * sampleMean)*NUM_SAMPLES/(NUM_SAMPLES-1);

                Assert.assertEquals(sampleMean, theoreticalMean, 10 * standardDeviationOfMean);
                Assert.assertEquals(sampleVariance, theoreticalVariance, 10e-4);
            }
        }
    }

    @Test
    public void testGaussian() {
        final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));
        for (final double theoreticalMean : Arrays.asList(0)) {
            for (final double precision : Arrays.asList(1.0)) {
                final double variance = 1/precision;

                //Note: this is the theoretical standard deviation of the sample mean given uncorrelated
                //samples.  The sample mean will have a greater variance here because samples are correlated.
                final double standardDeviationOfMean = Math.sqrt(variance / NUM_SAMPLES);

                final Function<Double, Double> logPDF = x -> -(precision/2)*(x-theoreticalMean)*(x-theoreticalMean);
                final AdaptiveMetropolisSampler sampler = new AdaptiveMetropolisSampler(INITIAL_GAUSSIAN_GUESS, INITIAL_STEP_SIZE,
                        Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
                final List<Double> samples = sampler.sample(rng, logPDF, NUM_SAMPLES, NUM_BURN_IN_STEPS);

                final double sampleMean = samples.stream().mapToDouble(x -> x).average().getAsDouble();
                final double sampleMeanSquare = samples.stream().mapToDouble(x -> x*x).average().getAsDouble();
                final double sampleVariance = (sampleMeanSquare - sampleMean * sampleMean)*NUM_SAMPLES/(NUM_SAMPLES-1);

                Assert.assertEquals(sampleMean, theoreticalMean, 6 * standardDeviationOfMean);
                Assert.assertEquals(sampleVariance, variance, variance/10);
            }
        }
    }
}