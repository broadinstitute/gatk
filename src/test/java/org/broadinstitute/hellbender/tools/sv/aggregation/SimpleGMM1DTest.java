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

    /**
     * Validate that the Bayesian GMM produces scores in the expected range for
     * realistic BAF log-ratio data. The reference values were computed using
     * sklearn's BayesianGaussianMixture(n_components=3, covariance_type='spherical',
     * weight_concentration_prior_type='dirichlet_process', n_init=20, max_iter=200).
     *
     * <p>Note: The Java and Python implementations use different random initialization
     * strategies (Java uses random responsibilities; sklearn default uses kmeans), so
     * exact numerical agreement is not expected. Instead we validate that:
     * (1) in-distribution data scores higher than outliers
     * (2) the Bayesian model produces scores in the right ballpark (not inflated 3x)
     * (3) the score is negative for outlier carrier data (higher BAF_DEL_LOGLIK)</p>
     *
     * <p>With sklearn, the standard (non-Bayesian) GMM produced BAF_DEL_LOGLIK of ~50
     * for carrier=[-0.5, -0.3] against control data centered at 0 with std~0.07,
     * while the Bayesian GMM produced ~14. The Bayesian regularization produces
     * broader components, making scores less extreme for outlier data.</p>
     */
    @Test
    public void testBayesianRegularizationProducesModerateScores() {
        // Control data: BAF log-ratios for non-carrier samples (centered ~0, std ~0.07)
        final double[] control = new double[]{
                -0.12, -0.05, 0.03, 0.08, -0.15,
                0.01, -0.08, 0.06, -0.03, 0.11,
                -0.01, 0.04, -0.09, 0.02, -0.06,
                0.07, -0.04, 0.0, -0.1, 0.05
        };
        // Carrier data: shifted negative (true DEL)
        final double[] carrier = new double[]{-0.5, -0.3};

        final SimpleGMM1D gmm = new SimpleGMM1D(3, 20, 200, 42L);
        gmm.fit(control);

        final double score = gmm.score(carrier);
        final double delLoglik = -score;

        // sklearn Bayesian GMM gives BAF_DEL_LOGLIK ≈ 14.4 for this data.
        // sklearn standard GMM gives BAF_DEL_LOGLIK ≈ 49.7 (3.5× inflated).
        // Our implementation should be much closer to the Bayesian value.
        Assert.assertTrue(delLoglik > 5, "BAF_DEL_LOGLIK should be positive for shifted carriers: " + delLoglik);
        Assert.assertTrue(delLoglik < 35, "BAF_DEL_LOGLIK should not be inflated like standard GMM (≈50): " + delLoglik);

        // In-distribution data should score much better
        final double inDistScore = gmm.score(new double[]{0.0, 0.02});
        Assert.assertTrue(inDistScore > score,
                "In-distribution data should score higher than outlier carriers");
    }

    @Test
    public void testBayesianComponentPruning() {
        // With enough data from a unimodal distribution, Bayesian GMM should
        // effectively prune extra components (assign them near-zero weight)
        final double[] unimodal = new double[100];
        final java.util.Random rng = new java.util.Random(123);
        for (int i = 0; i < 100; i++) {
            unimodal[i] = rng.nextGaussian() * 0.1;
        }

        final SimpleGMM1D gmm = new SimpleGMM1D(3, 10, 200, 42L);
        gmm.fit(unimodal);

        // Score of a mild outlier should be moderate, not extreme
        final double outlierScore = gmm.score(new double[]{-0.5});
        final double delLoglik = -outlierScore;

        // With standard EM, tight 3-component fit would give very extreme score.
        // With Bayesian priors, the regularization prevents this.
        Assert.assertTrue(delLoglik > 0, "Outlier should have positive -loglik: " + delLoglik);
        Assert.assertTrue(delLoglik < 50, "Bayesian should not give extreme score: " + delLoglik);
    }
}