package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import com.google.common.primitives.Doubles;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class Mutect2FilteringEngineUnitTest {
    @DataProvider(name = "falsePositiveRateData")
    public Object[][] makeFalsePositiveRateData() {
        return new Object[][]{
                {new double[]{0.01, 0.01, 0.05, 0.01, 0.3, 0.4}, 0.1, 0.3},
                // If there are ties e.g. 0.9's below, set the threshold such that the expected false positive rate does not exceed the requested max false positive rate
                // Filtering variants with p > 0.9 gives you the expected false positive rate of 0.076 > 0.5, so the threshold must be 0.01
                {new double[]{0.01, 0.01, 0.9, 0.9, 0.9}, 0.05, 0.01},
                {new double[]{0.01, 0.01, 0.01, 0.01, 0.01}, 0.05, 1.0}, // The FPR never exceeds the max FPR so do not filter
                {new double[]{0.99, 0.99, 0.99, 1.0, 1.0}, 0.05, 0.0}, // Impossible to meet the FPR with this data so filter everything
                {new double[]{0.01, 0.01, 0.05, 0.01, 0.3, 0.4}, 0.0, 0.0},
        };
    }

    @Test(dataProvider = "falsePositiveRateData")
    public void testCalculateThreshold(final double[] posteriors,
                                       final double maxErrorRate,
                                       final double expectedThreshold) {
        final double threshold = ThresholdCalculator.calculateThresholdBasedOnFalseDiscoveryRate(Doubles.asList(posteriors), maxErrorRate);
        Assert.assertEquals(threshold, expectedThreshold);
    }
}