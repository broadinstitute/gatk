package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

/**
 * Tests the sampler classes of the allele-fraction model.
 *
 *  @author David Benjamin
 */
public class AlleleFractionSamplersUnitTest {
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
        final double MINOR_FRACTION_TOLERANCE = Math.sqrt(5.0 / (12.0 * 3.14 * averageHetsPerSegment * averageDepth));

        final AlleleFractionSimulatedData SIMULATED_DATA = new AlleleFractionSimulatedData(averageHetsPerSegment, numSegments, averageDepth, meanBias, biasVariance, outlierProbability);
        final AlleleFractionState INITIAL_STATE = SIMULATED_DATA.getTrueState();
        final AlleleFractionData DATA = new AlleleFractionData(SIMULATED_DATA.getSegmentedModel());

        final AlleleFractionSamplers.MeanBiasSampler meanBiasSampler =
                new AlleleFractionSamplers.MeanBiasSampler(INITIAL_STATE, 0.01);
        final AlleleFractionSamplers.BiasVarianceSampler biasVarianceSampler =
                new AlleleFractionSamplers.BiasVarianceSampler(INITIAL_STATE, 0.01);
        final AlleleFractionSamplers.OutlierProbabilitySampler outlierProbabilitySampler =
                new AlleleFractionSamplers.OutlierProbabilitySampler(INITIAL_STATE, 0.01);
        final AlleleFractionSamplers.MinorFractionsSampler minorFractionsSampler =
                new AlleleFractionSamplers.MinorFractionsSampler(INITIAL_STATE, Collections.nCopies(numSegments, 0.01));

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