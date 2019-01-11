package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.distribution.IntegerDistribution;
import org.broadinstitute.hellbender.tools.spark.utils.IntHistogram;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Random;

public class IntHistogramTest extends GATKBaseTest {
    private static final int MAX_TRACKED_VALUE = 2000;
    private static final float SIGNIFICANCE = .05f;

    @Test
    void testTrim() {
        final IntHistogram intHistogram = new IntHistogram(8);
        intHistogram.addObservation(0);
        intHistogram.addObservations(1, 3L);
        intHistogram.addObservations(2, 5L);
        intHistogram.addObservations(3, 3L);
        intHistogram.addObservation(4);
        intHistogram.addObservation(10);
        Assert.assertEquals(intHistogram.trim().getMaximumTrackedValue(), 3);
    }

    @Test
    void testTrimTo0() {
        final IntHistogram intHistogram = new IntHistogram(2);
        intHistogram.addObservations(0, 10L);
        intHistogram.addObservation(3);
        Assert.assertEquals(intHistogram.trim().getMaximumTrackedValue(), 0);
    }

    @Test
    void testTrimEarlyBailout() {
        final IntHistogram intHistogram = new IntHistogram(0);
        intHistogram.addObservations(0, 10L);
        intHistogram.addObservations(1, 10L);
        Assert.assertTrue(intHistogram.trim() == intHistogram);
    }

    // TODO: figure out what the data should really look like
    @Test
    public void testSmallInsertion() {
        final IntHistogram.CDF cdf = genNormalSample(480, 25, 10000000).getCDF();
        final IntHistogram sample = genNormalSample(480, 25, 25);
        for ( int idx = 75; idx != MAX_TRACKED_VALUE; ++idx )
            sample.addObservations(idx-50, sample.getNObservations(idx));
        Assert.assertTrue(cdf.isDifferentByKSStatistic(sample, .01f));
    }

    @Test
    public void testNormalDistribution() {
        final IntHistogram.CDF cdf = genNormalSample(480, 25, 10000).getCDF();
        Assert.assertFalse(cdf.isDifferentByKSStatistic(genNormalSample(480, 25, 100), SIGNIFICANCE));
        // these numbers were jiggled by hand until the difference was barely detectable
        Assert.assertTrue(cdf.isDifferentByKSStatistic(genNormalSample(489, 25, 100), SIGNIFICANCE));
        Assert.assertTrue(cdf.isDifferentByKSStatistic(genNormalSample(470, 25, 100), SIGNIFICANCE));
        Assert.assertTrue(cdf.isDifferentByKSStatistic(genNormalSample(480, 15, 100), SIGNIFICANCE));
        Assert.assertTrue(cdf.isDifferentByKSStatistic(genNormalSample(480, 42, 100), SIGNIFICANCE));
    }

    @Test
    public void testLogNormalDistribution() {
        final IntHistogram.CDF cdf = genNormalSample(480, 25, 10000).getCDF();
        Assert.assertFalse(cdf.isDifferentByKSStatistic(genLogNormalSample(480, 25, 100), SIGNIFICANCE));
        Assert.assertTrue(cdf.isDifferentByKSStatistic(genLogNormalSample(490, 25, 100), SIGNIFICANCE));
        Assert.assertTrue(cdf.isDifferentByKSStatistic(genLogNormalSample(471, 25, 100), SIGNIFICANCE));
        Assert.assertTrue(cdf.isDifferentByKSStatistic(genLogNormalSample(480, 15, 100), SIGNIFICANCE));
        Assert.assertTrue(cdf.isDifferentByKSStatistic(genLogNormalSample(480, 42, 100), SIGNIFICANCE));
    }

    @Test
    public void testBiModalDistribution() {
        final IntHistogram.CDF cdf = genNormalSample(480, 25, 10000).getCDF();
        final IntHistogram sample1 = genNormalSample(480, 25, 50);
        final IntHistogram sample2 = genNormalSample(500, 25, 50);
        sample1.addObservations(sample2);
        Assert.assertTrue(cdf.isDifferentByKSStatistic(sample1, SIGNIFICANCE));
    }

    @Test
    public void testEmpiricalDistributionWithoutSmoothingSampling() {
        final IntHistogram largeSample = genNormalSample(480, 25, 10000);
        final IntegerDistribution dist = largeSample.empiricalDistribution(0);
        final IntHistogram distSample = new IntHistogram(MAX_TRACKED_VALUE);
        Arrays.stream(dist.sample(1000)).forEach(distSample::addObservation);
        Assert.assertFalse(largeSample.getCDF().isDifferentByKSStatistic(distSample, SIGNIFICANCE));
    }

    @Test(dataProvider = "smoothingValues")
    public void testEmpiricalDistributionSmoothing(final int smoothing) {
        final IntHistogram largeSample = genNormalSample(480, 25, 10000);
        final IntegerDistribution dist = largeSample.empiricalDistribution(smoothing);
        final long smoothedNumberOfObservations = largeSample.getMaximumTrackedValue() * smoothing + largeSample.getTotalObservations();
        double cumulative = 0;
        double expectation = 0;
        double sqExpectation = 0;
        for (int i = 0; i <= largeSample.getMaximumTrackedValue(); i++) {
            final double distProb = dist.probability(i);
            Assert.assertEquals(distProb, (largeSample.getNObservations(i) + smoothing) / (double) smoothedNumberOfObservations, 0.0001);
            cumulative += distProb;
            Assert.assertEquals(dist.cumulativeProbability(i), cumulative, 0.00001);
            expectation += distProb * i;
            sqExpectation += i * distProb * i;
        }
        Assert.assertEquals(dist.getNumericalMean(), expectation, 0.00001);
        Assert.assertEquals(dist.getNumericalVariance(), sqExpectation - expectation * expectation, 0.00001);
    }

    @DataProvider
    public Object[][] smoothingValues() {
        return new Object[][] { { 0 }, { 1 }, { 2 }, { 13 }, {100 }};
    }

    public static IntHistogram genNormalSample( final int mean, final int stdDev, final int nSamples ) {
        Random random = new Random(47L);
        final IntHistogram histogram = new IntHistogram(MAX_TRACKED_VALUE);
        for ( int sample = 0; sample != nSamples; ++sample ) {
            final int val = (int)Math.round(mean + stdDev * random.nextGaussian());
            // this is just for testing -- no need to be super correct about the distribution
            if ( val < 0 ) { sample -= 1; continue; }
            histogram.addObservation(val);
        }
        return histogram;
    }

    public static IntHistogram genLogNormalSample( final int mean, final int stdDev, final int nSamples ) {
        Random random = new Random(47L);
        final double phi = Math.sqrt((double)stdDev * stdDev + (double)mean * mean);
        final double mu = Math.log((double)mean * mean / phi);
        final double sigma = Math.sqrt(Math.log(phi * phi / mean / mean));
        final IntHistogram histogram = new IntHistogram(MAX_TRACKED_VALUE);
        for ( int sample = 0; sample != nSamples; ++sample ) {
            final int val = (int)Math.round(Math.exp(mu + sigma * random.nextGaussian()));
            // this is just for testing -- no need to be super correct about the distribution
            if ( val < 0 ) { sample -= 1; continue; }
            histogram.addObservation(val);
        }
        return histogram;
    }
}
