package org.broadinstitute.hellbender.tools.sv.aggregation;

import org.testng.Assert;
import org.testng.annotations.Test;

public class SimpleGMM1DTest {

    @Test
    public void testConstructorAndScoreValidation() {
        Assert.assertThrows(IllegalArgumentException.class, () -> new SimpleGMM1D(0, 1, 1, 1));
        Assert.assertThrows(IllegalArgumentException.class, () -> new SimpleGMM1D(1, 0, 1, 1));
        Assert.assertThrows(IllegalArgumentException.class, () -> new SimpleGMM1D(1, 1, 0, 1));

        final SimpleGMM1D gmm = new SimpleGMM1D(3, 3, 50, 13L);
        Assert.assertThrows(IllegalArgumentException.class, () -> gmm.fit(new double[]{}));
        Assert.assertThrows(IllegalArgumentException.class, () -> gmm.score(new double[]{0.0}));
    }

    @Test
    public void testFitAndScoreAreDeterministicForSeed() {
        final double[] train = new double[]{-2.1, -2.0, -1.9, -0.1, 0.0, 0.1, 1.9, 2.0, 2.1, -2.05, 0.05, 2.05};
        final double[] eval = new double[]{-2.0, 0.0, 2.0};

        final SimpleGMM1D gmm1 = new SimpleGMM1D(3, 10, 100, 7L);
        final SimpleGMM1D gmm2 = new SimpleGMM1D(3, 10, 100, 7L);
        gmm1.fit(train);
        gmm2.fit(train);

        Assert.assertEquals(gmm1.score(eval), gmm2.score(eval), 1e-12);
    }

    @Test
    public void testInDistributionScoresBetterThanOutliers() {
        final double[] train = new double[]{-2.2, -2.1, -2.0, -1.9, -0.2, -0.1, 0.0, 0.1, 1.9, 2.0, 2.1, 2.2};
        final SimpleGMM1D gmm = new SimpleGMM1D(3, 10, 100, 11L);
        gmm.fit(train);

        final double inDistributionScore = gmm.score(new double[]{-2.0, 0.0, 2.0});
        final double outlierScore = gmm.score(new double[]{8.0, 9.0, 10.0});

        Assert.assertTrue(inDistributionScore > outlierScore);
    }

    @Test
    public void testFitOnConstantDataProducesFiniteScore() {
        final double[] train = new double[]{0.25, 0.25, 0.25, 0.25, 0.25, 0.25};
        final SimpleGMM1D gmm = new SimpleGMM1D(2, 5, 50, 3L);
        gmm.fit(train);

        final double score = gmm.score(new double[]{0.25, 0.25});
        Assert.assertTrue(Double.isFinite(score));
    }
}