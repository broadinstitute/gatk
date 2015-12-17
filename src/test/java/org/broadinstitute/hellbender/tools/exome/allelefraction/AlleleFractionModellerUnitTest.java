package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.commons.math3.special.Gamma;
import org.broadinstitute.hellbender.tools.exome.AllelicCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

import static java.lang.Math.log;

/**
 * Test the likelihood methods and MCMC inference of the {@link AlleleFractionModeller} <p>
 *
 * We subject the likelihood (with bias marginalized out) to several tests.  Recall that the exact likelihood is:
 * <p>
 * if indicator == ALT_MINOR:
 * <p>
 * log { [beta^alpha / Gamma(alpha)][(1-pi)/2] * int_{0 to infty} f^a * (1-f)^r * lambda^(alpha + r - 1) * exp(-beta*lambda)/(f + (1-f)*lambda)^n d lambda }
 * ,<p> if indicator == REF_MINOR same as ALT_MINOR but with f <--> 1 - f </p>
 * ,<p>if indicator == OUTLIER, log {pi * a!r!/(n+1)!} </p>
 * where alpha = mu*beta and n = a + r.
 * @author David Benjamin
 */
public final class AlleleFractionModellerUnitTest {
    private static final SimpleInterval DUMMY = new SimpleInterval("dummy", 1, 2);
    private static final double EPSILON = 1e-10;

    //if f is very close to 0 we have an analytic result for comparison
    @Test
    public void testHetLogLikelihoodMinorFractionNearZero() {
        final double pi = 0.01; //pi is just a prefactor so we don't need to test it thoroughly here
        for (final double f : Arrays.asList(1e-6, 1e-7, 1e-8)) {
            for (final double mean : Arrays.asList(0.9, 1.0, 1.1)) {
                for (final double variance : Arrays.asList(0.01, 0.005, 0.001)) {
                    final double alpha = mean * mean / variance;
                    final double beta = mean / variance;
                    final AlleleFractionState state = AlleleFractionModeller.makeSingleSegmentState(mean, variance, pi, f);
                    for (final int a : Arrays.asList(1, 2, 3)) {  //alt count
                        for (final int r : Arrays.asList(50, 100, 200)) { //ref count
                            final AllelicCount count = new AllelicCount(DUMMY, r, a);
                            final double actual = AlleleFractionModeller.hetLogLikelihood(state, 0, count, AlleleFractionIndicator.ALT_MINOR);
                            final double expected = a * log(beta) + Gamma.logGamma(alpha - a) - Gamma.logGamma(alpha)
                                    + log((1 - pi) / 2) + a * log(f / (1 - f));
                            Assert.assertEquals(actual, expected, 1e-3);
                        }
                    }
                }
            }
        }
    }

    // if f is very close to 1 we have an analytic result for comparison
    @Test
    public void testHetLogLikelihoodMinorFractionNearOne() {
        final double pi = 0.01; //pi is just a prefactor so we don't need to test it thoroughly here
        for (final double f : Arrays.asList(1 - 1e-6, 1 - 1e-7, 1 - 1e-8)) {
            for (final double mean : Arrays.asList(0.9, 1.0, 1.1)) {
                for (final double variance : Arrays.asList(0.01, 0.005, 0.001)) {
                    final double alpha = mean * mean / variance;
                    final double beta = mean / variance;
                    final AlleleFractionState state = AlleleFractionModeller.makeSingleSegmentState(mean, variance, pi, f);
                    for (final int a : Arrays.asList(1, 10, 20)) {  //alt count
                        for (final int r : Arrays.asList(1, 10, 20)) { //ref count
                            final AllelicCount count = new AllelicCount(DUMMY, r, a);
                            final double actual = AlleleFractionModeller.hetLogLikelihood(state, 0, count, AlleleFractionIndicator.ALT_MINOR);
                            final double expected = -r * log(beta) + Gamma.logGamma(alpha + r) - Gamma.logGamma(alpha)
                                    + log((1 - pi) / 2) - r * log(f / (1 - f));
                            Assert.assertEquals(actual, expected,1e-4);
                        }
                    }
                }
            }
        }
    }

    //if variance is tiny we can approximate lambda ~ mu in the tricky part of the integral to get an analytic result
    @Test
    public void testHetLogLikelihoodTightDistribution() {
        final double pi = 0.01; //pi is just a prefactor so we don't need to test it thoroughly here
        for (final double f : Arrays.asList(0.1, 0.2, 0.3)) {
            for (final double mean : Arrays.asList(0.9, 1.0, 1.1)) {
                for (final double variance : Arrays.asList(1e-6, 1e-7, 1e-8)) {
                    final AlleleFractionState state = AlleleFractionModeller.makeSingleSegmentState(mean, variance, pi, f);
                    for (final int a : Arrays.asList(1, 10, 20)) {  //alt count
                        for (final int r : Arrays.asList(1, 10, 20)) { //ref count
                            final AllelicCount count = new AllelicCount(DUMMY, r, a);
                            final double actual = AlleleFractionModeller.hetLogLikelihood(state, 0, count, AlleleFractionIndicator.ALT_MINOR);
                            final double expected = log((1 - pi) / 2) + a * log(f) + r * log(1-f) + r * log(mean) - (a+r) * log(f + (1-f)*mean);
                            Assert.assertEquals(actual, expected, 1e-3);
                        }
                    }
                }
            }
        }
    }

    //ALT_MINOR <--> REF_MINOR is equivalent to f <--> 1 - f
    @Test
    public void testRefMinor() {
        final double pi = 0.01; //pi is just a prefactor so we don't need to test it thoroughly here
        for (final double f : Arrays.asList(0.1, 0.2, 0.3)) {
            for (final double mean : Arrays.asList(0.9, 1.0, 1.1)) {
                for (final double variance : Arrays.asList(0.02, 0.01)) {
                    final AlleleFractionState altMinorState = AlleleFractionModeller.makeSingleSegmentState(mean, variance, pi, f);
                    final AlleleFractionState refMinorState = AlleleFractionModeller.makeSingleSegmentState(mean, variance, pi, 1 - f);
                    for (final int a : Arrays.asList(1, 10, 20)) {  //alt count
                        for (final int r : Arrays.asList(1, 10, 20)) { //ref count
                            final AllelicCount count = new AllelicCount(DUMMY, r, a);
                            final double altMinorLk = AlleleFractionModeller.hetLogLikelihood(altMinorState, 0, count, AlleleFractionIndicator.ALT_MINOR);
                            final double refMinorLk = AlleleFractionModeller.hetLogLikelihood(refMinorState, 0, count, AlleleFractionIndicator.REF_MINOR);
                            Assert.assertEquals(altMinorLk, refMinorLk, 1e-10);
                        }
                    }
                }
            }
        }
    }

    @Test
    public void testHetLogLikelihoodOutlierProbabilityDependence() {
        final AllelicCount count = new AllelicCount(DUMMY, 11, 37);
        final double f = 0.25;
        final double mean = 1.0;
        final double variance = 0.01;
        final double pi1 = 0.1;
        final double pi2 = 0.2;
        final double pi3 = 0.3;

        final AlleleFractionState state1 = AlleleFractionModeller.makeSingleSegmentState(mean, variance, pi1, f);
        final AlleleFractionState state2 = AlleleFractionModeller.makeSingleSegmentState(mean, variance, pi2, f);
        final AlleleFractionState state3 = AlleleFractionModeller.makeSingleSegmentState(mean, variance, pi3, f);

        final double lk1 = AlleleFractionModeller.hetLogLikelihood(state1, 0, count, AlleleFractionIndicator.ALT_MINOR);
        final double lk2 = AlleleFractionModeller.hetLogLikelihood(state2, 0, count, AlleleFractionIndicator.ALT_MINOR);
        final double lk3 = AlleleFractionModeller.hetLogLikelihood(state3, 0, count, AlleleFractionIndicator.ALT_MINOR);

        Assert.assertEquals(lk2 - lk1, log(1 - pi2) - log(1 - pi1), EPSILON);
        Assert.assertEquals(lk3 - lk2, log(1 - pi3) - log(1 - pi2), EPSILON);
    }

    //when read counts are very high the posterior on lambda ought to be very tight, so we can check moments analytically
    @Test
    public void testBiasPosteriorMoment() {
        final double pi = 0.01;

        for (final double f : Arrays.asList(0.1, 0.2, 0.3)) {
            for (final double meanBias : Arrays.asList(0.8, 0.9, 1.0, 1.1, 1.2)) {
                final double biasVariance = 0.01;  //prior average bias agrees with the MLE
                final AlleleFractionState state = AlleleFractionModeller.makeSingleSegmentState(meanBias, biasVariance, pi, f);
                for (final double n : Arrays.asList(1000, 2000, 3000)) {
                    final int a = (int) (n * f / (f + (1 - f) * meanBias));
                    final int r = (int) (n - a);
                    final AllelicCount count = new AllelicCount(DUMMY, r, a);
                    final double actualMean = AlleleFractionModeller.biasPosteriorMoment(state, 0, count, AlleleFractionIndicator.ALT_MINOR, 1);
                    final double actualMeanSquare = AlleleFractionModeller.biasPosteriorMoment(state, 0, count, AlleleFractionIndicator.ALT_MINOR, 2);
                    Assert.assertEquals(actualMean, meanBias, 1e-2);
                    Assert.assertEquals(actualMeanSquare, meanBias*meanBias, 5e-2);
                }
            }
        }
    }

    @Test
    public void testMCMC() {
        final int numSamples = 300;
        final int numBurnIn = 100;

        final double averageHetsPerSegment = 20;
        final int numSegments = 100;
        final int averageDepth = 50;

        final double meanBias = 1.1;
        final double biasVariance = 0.01;
        final double outlierProbability = 0.02;

        // note: the following tolerances could actually be made much smaller if we used more segments and/or
        // more hets -- most of the error is the sampling error of a finite simulated data set, not numerical error of MCMC
        final double minorFractionTolerance = 0.03;
        final double meanBiasTolerance = 0.03;
        final double biasVarianceTolerance = 0.04;
        final double outlierProbabilityTolerance = 0.03;
        final AlleleFractionSimulatedData simulatedData = new AlleleFractionSimulatedData(averageHetsPerSegment, numSegments,
                averageDepth, meanBias, biasVariance, outlierProbability);

        final AlleleFractionModeller model = new AlleleFractionModeller(simulatedData.getSegmentedModel());
        model.fitMCMC(numSamples, numBurnIn);


        final List<Double> meanBiasSamples = model.getmeanBiasSamples();
        Assert.assertEquals(meanBiasSamples.size(), numSamples - numBurnIn);

        final List<Double> biasVarianceSamples = model.getBiasVarianceSamples();
        Assert.assertEquals(biasVarianceSamples.size(), numSamples - numBurnIn);

        final List<Double> outlierProbabilitySamples = model.getOutlierProbabilitySamples();
        Assert.assertEquals(outlierProbabilitySamples.size(), numSamples - numBurnIn);

        final List<AlleleFractionState.MinorFractions> minorFractionsSamples = model.getMinorFractionsSamples();
        Assert.assertEquals(minorFractionsSamples.size(), numSamples - numBurnIn);
        for (final AlleleFractionState.MinorFractions sample : minorFractionsSamples) {
            Assert.assertEquals(sample.size(), numSegments);
        }

        final List<List<Double>> minorFractionsSamplesBySegment = model.getMinorFractionSamplesBySegment();

        final double mcmcMeanBias = meanBiasSamples.stream().mapToDouble(x -> x).average().getAsDouble();
        final double mcmcBiasVariance = biasVarianceSamples.stream().mapToDouble(x -> x).average().getAsDouble();
        final double mcmcOutlierProbabilityr = outlierProbabilitySamples.stream().mapToDouble(x -> x).average().getAsDouble();
        final List<Double> mcmcMinorFractions = minorFractionsSamplesBySegment
                .stream().map(list -> list.stream().mapToDouble(x -> x).average().getAsDouble())
                .collect(Collectors.toList());

        double totalSegmentError = 0.0;
        for (int segment = 0; segment < numSegments; segment++) {
            totalSegmentError += Math.abs(mcmcMinorFractions.get(segment) - simulatedData.getTrueState().minorFractionInSegment(segment));
        }

        Assert.assertEquals(mcmcMeanBias, meanBias, meanBiasTolerance);
        Assert.assertEquals(mcmcBiasVariance, biasVariance, biasVarianceTolerance);
        Assert.assertEquals(mcmcOutlierProbabilityr, outlierProbability, outlierProbabilityTolerance);
        Assert.assertEquals(totalSegmentError / numSegments, 0.0, minorFractionTolerance);
    }

    @Test
    public void testSamplers() {
        final int randomSeed = 15;
        final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(randomSeed));
        final int numSamples = 500;
        final double averageHetsPerSegment = 20;
        final int numSegments = 100;
        final int averageDepth = 50;
        final double meanBias = 1.1;
        final double biasVariance = 0.01;
        final double outlierProbability = 0.02;

        // note: the following tolerances could actually be made much smaller if we used more segments and/or
        // more hets -- most of the error is the sampling error of a finite simulated data set, not numerical error of MCMC
        final double meanBiasTolerance = 0.02;
        final double biasVarianceTolerance = 0.01;
        final double outlierProbabilityTolerance = 0.02;

        // as for the minor fraction tolerance, the "truth data" assigns some value of f for each segment
        // but we draw actual read counts for each het from a binomial (ignoring bias) with total number ~ d = average depth
        // and minor allele probability f.  this implies a variance ~ f(1-f)*d in the total minor read count for a het
        // and a variance f(1-f)*d*N, where N is average hets per segment, per segment.  We divide by the total number of reads,
        // roughly N*d, to get an empirical minor fraction of the segment.  This estimator thus has variance f(1-f)/(N*d).
        // averaging from 0 < f < 1/2 gives an average variance of 5/(24*N*d).  Assuming normality,
        // the absolute value of this sampling error (i.e. inherent to the simulated data and having nothing to do with
        // the MCMC or the model) has mean sqrt[5 / (12*pi*N*d)]
        final double MINOR_FRACTION_SAMPLING_ERROR = Math.sqrt(5.0 / (12.0 * 3.14 * averageHetsPerSegment * averageDepth));
        final double MINOR_FRACTION_TOLERANCE = MINOR_FRACTION_SAMPLING_ERROR;

        final AlleleFractionSimulatedData SIMULATED_DATA = new AlleleFractionSimulatedData(averageHetsPerSegment, numSegments, averageDepth, meanBias, biasVariance, outlierProbability);
        final AlleleFractionState INITIAL_STATE = SIMULATED_DATA.getTrueState();
        final AlleleFractionData DATA = new AlleleFractionData(SIMULATED_DATA.getSegmentedModel());

        final AlleleFractionModeller.MeanBiasSampler meanBiasSampler =
                new AlleleFractionModeller.MeanBiasSampler(INITIAL_STATE.meanBias());
        final AlleleFractionModeller.BiasVarianceSampler biasVarianceSampler =
                new AlleleFractionModeller.BiasVarianceSampler(INITIAL_STATE.biasVariance());
        final AlleleFractionModeller.OutlierProbabilitySampler outlierProbabilitySampler =
                new AlleleFractionModeller.OutlierProbabilitySampler(INITIAL_STATE.outlierProbability());
        final AlleleFractionModeller.MinorFractionsSampler minorFractionsSampler =
                new AlleleFractionModeller.MinorFractionsSampler(INITIAL_STATE.minorFractions());

        final AlleleFractionState state = INITIAL_STATE.copy(AlleleFractionState.class);
        final List<Double> meanBiasSamples = new ArrayList<>();
        final List<Double> biasVarianceSamples = new ArrayList<>();
        final List<Double> outlierProbabilitySamples = new ArrayList<>();
        final List<AlleleFractionState.MinorFractions> minorFractionsSamples = new ArrayList<>();
        for (int n = 0; n < numSamples; n++) {
            meanBiasSamples.add(meanBiasSampler.sample(rng, state, DATA));
            biasVarianceSamples.add(biasVarianceSampler.sample(rng, state, DATA));
            outlierProbabilitySamples.add(outlierProbabilitySampler.sample(rng, state, DATA));
            minorFractionsSamples.add(minorFractionsSampler.sample(rng, state, DATA));
        }

        final double estimatedMeanBias = meanBiasSamples.stream().mapToDouble(x -> x).average().getAsDouble();
        Assert.assertEquals(estimatedMeanBias, meanBias, meanBiasTolerance);

        final double estimatedBiasVariance = biasVarianceSamples.stream().mapToDouble(x -> x).average().getAsDouble();
        Assert.assertEquals(estimatedBiasVariance, biasVariance, biasVarianceTolerance);

        final double estimatedOutlierProbability = outlierProbabilitySamples.stream().mapToDouble(x -> x).average().getAsDouble();
        Assert.assertEquals(estimatedOutlierProbability, outlierProbability, outlierProbabilityTolerance);

        final List<List<Double>> samplesBySegment = new ArrayList<>();
        for (int segment = 0; segment < numSegments; segment++) {
            final int seg = segment;    //to allow use in the following lambda
            samplesBySegment.add(minorFractionsSamples.stream().map(s -> s.get(seg)).collect(Collectors.toList()));
        }

        final List<Double> estimatedMinorFractions = samplesBySegment
                .stream().map(list -> list.stream().mapToDouble(x -> x).average().getAsDouble())
                .collect(Collectors.toList());

        double totalSegmentError = 0.0;
        for (int segment = 0; segment < numSegments; segment++) {
            totalSegmentError += Math.abs(estimatedMinorFractions.get(segment) - INITIAL_STATE.minorFractionInSegment(segment));
        }

        Assert.assertEquals(totalSegmentError / numSegments, 0.0, MINOR_FRACTION_TOLERANCE);
    }
}