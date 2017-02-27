package org.broadinstitute.hellbender.utils;

import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.Random;

public final class KolmogorovSmirnovCalculatorTest {
    private static final int ARR_LEN = 1000;
    private static Random random;

    @BeforeMethod
    public void reseedRandom() {
        random = new Random(47L);
    }

    @Test
    public void testNormalDistribution() {
        final KolmogorovSmirnovCalculator calc =
                new KolmogorovSmirnovCalculator(genNormalSample(480, 25, 10000), .05f);
        Assert.assertFalse(calc.isDifferent(genNormalSample(480, 25, 100)));
        // these numbers were jiggled by hand until the difference was barely detectable
        Assert.assertTrue(calc.isDifferent(genNormalSample(489, 25, 100)));
        Assert.assertTrue(calc.isDifferent(genNormalSample(473, 25, 100)));
        Assert.assertTrue(calc.isDifferent(genNormalSample(480, 15, 100)));
        Assert.assertTrue(calc.isDifferent(genNormalSample(480, 42, 100)));
    }

    @Test
    public void testLogNormalDistribution() {
        final KolmogorovSmirnovCalculator calc =
                new KolmogorovSmirnovCalculator(genLogNormalSample(480, 25, 10000), .05f);
        Assert.assertFalse(calc.isDifferent(genLogNormalSample(480, 25, 100)));
        Assert.assertTrue(calc.isDifferent(genLogNormalSample(489, 25, 100)));
        Assert.assertTrue(calc.isDifferent(genLogNormalSample(473, 25, 100)));
        Assert.assertTrue(calc.isDifferent(genLogNormalSample(480, 15, 100)));
        Assert.assertTrue(calc.isDifferent(genLogNormalSample(480, 42, 100)));
    }

    @Test
    public void testBiModalDistribution() {
        final KolmogorovSmirnovCalculator calc =
                new KolmogorovSmirnovCalculator(genNormalSample(480, 25, 10000), .05f);
        final long[] sample1 = genNormalSample(480, 25, 50);
        final long[] sample2 = genNormalSample(500, 25, 50);
        for ( int idx = 0; idx != ARR_LEN; ++idx ) {
            sample1[idx] += sample2[idx];
        }
        Assert.assertTrue(calc.isDifferent(sample1));
    }

    private static long[] genNormalSample( final int mean, final int stdDev, final int nSamples ) {
        final long[] histogram = new long[ARR_LEN];
        for ( int sample = 0; sample != nSamples; ++sample ) {
            int bin = (int)Math.round(mean + stdDev * random.nextGaussian());
            if ( bin >= ARR_LEN ) bin = ARR_LEN - 1;
            histogram[bin] += 1;
        }
        return histogram;
    }

    private static long[] genLogNormalSample( final int mean, final int stdDev, final int nSamples ) {
        final double phi = Math.sqrt((double)stdDev*stdDev + (double)mean*mean);
        final double mu = Math.log((double)mean*mean/phi);
        final double sigma = Math.sqrt(Math.log(phi*phi/mean/mean));
        final long[] histogram = new long[ARR_LEN];
        for ( int sample = 0; sample != nSamples; ++sample ) {
            int bin = (int)Math.round(Math.exp(mu + sigma * random.nextGaussian()));
            if ( bin >= ARR_LEN ) bin = ARR_LEN - 1;
            histogram[bin] += 1;
        }
        return histogram;
    }
}
