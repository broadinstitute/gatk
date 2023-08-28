package org.broadinstitute.hellbender.utils.clustering;

import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;

/**
 * Unit tests for static helper methods; checked against python implementations in sklearn.
 */
public final class BayesianGaussianMixtureUtilsUnitTest {

    private static final double RELATIVE_SYMMETRY_THRESHOLD = 1E-6;
    private static final double ABSOLUTE_POSITIVITY_THRESHOLD = 1E-10;
    private static final double EPSILON = 1E-10;

    @Test
    public void testLogDirichletNorm() {
        final double epsilon = 1E-10;
        Assert.assertEquals(
                BayesianGaussianMixtureUtils.logDirichletNorm(new ArrayRealVector(new double[]{1., 2., 3., 4.})),
                10.31692083029347, epsilon);
        Assert.assertEquals(
                BayesianGaussianMixtureUtils.logDirichletNorm(new ArrayRealVector(new double[]{0.1, 0.2, 0.3, 0.4})),
                -5.669252286684849, epsilon);
    }

    @Test
    public void testLogWishartNorm() {
        final double epsilon = 1E-6;
        Assert.assertEquals(
                BayesianGaussianMixtureUtils.logWishartNorm(
                        new ArrayRealVector(new double[]{3., 4., 5., 6.}),
                        new ArrayRealVector(new double[]{0.1, 0.2, 0.3, 0.4}),
                        2).toArray(), // nFeatures must be <= all elements of degreesOfFreedom
                new double[]{-2.2586593 , -3.45180648, -5.25041877, -7.53671313}, epsilon);
    }

    @Test
    public void testComputePrecisionCholesky() {
        final double epsilon = 1E-6;
        final List<RealMatrix> result = BayesianGaussianMixtureUtils.computePrecisionCholesky(
                Arrays.asList(  // must be symmetric, positive semidefinite covariance matrices
                        new Array2DRowRealMatrix(new double[][]{
                                {2., 1.},
                                {1., 2.}}),
                        new Array2DRowRealMatrix(new double[][]{
                                {4., 3.},
                                {3., 4.}}),
                        new Array2DRowRealMatrix(new double[][]{
                                {6., 5.},
                                {5., 6.}})),
                RELATIVE_SYMMETRY_THRESHOLD,
                ABSOLUTE_POSITIVITY_THRESHOLD);
        final List<RealMatrix> expected =
                Arrays.asList(
                        new Array2DRowRealMatrix(new double[][]{
                                {0.70710678, -0.40824829},
                                {0., 0.81649658}}),
                        new Array2DRowRealMatrix(new double[][]{
                                {0.5, -0.56694671},
                                {0., 0.75592895}}),
                        new Array2DRowRealMatrix(new double[][]{
                                {0.40824829, -0.61545745},
                                {0., 0.73854895}}));
        IntStream.range(0, result.size()).forEach(k ->
                Assert.assertTrue(result.get(k).subtract(expected.get(k)).getNorm() < epsilon));
    }

    @Test
    public void testLogDetCholesky() {
        final double epsilon = 1E-6;
        final RealVector result = BayesianGaussianMixtureUtils.computeLogDetCholesky(
                Arrays.asList(  // must be symmetric, positive semidefinite covariance matrices
                        new Array2DRowRealMatrix(new double[][]{
                                {2., 1.},
                                {1., 2.}}),
                        new Array2DRowRealMatrix(new double[][]{
                                {4., 3.},
                                {3., 4.}}),
                        new Array2DRowRealMatrix(new double[][]{
                                {6., 5.},
                                {5., 6.}})));
        final RealVector expected = new ArrayRealVector(new double[]{1.38629436, 2.77258872, 3.58351894});
        Assert.assertEquals(result.toArray(), expected.toArray(), epsilon);
    }

    @Test
    public void testEstimateLogGaussianProb() {
        final double epsilon = 1E-6;
        // nSamples = 5, nFeatures = 2, nComponents = 3
        final RealMatrix X = new Array2DRowRealMatrix(new double[][]{
                {1., 2.},
                {2., 3.},
                {3., 4.},
                {4., 5.},
                {5., 6.}});
        final List<RealVector> means = Arrays.asList(
                new ArrayRealVector(new double[]{1., 2.}),
                new ArrayRealVector(new double[]{3., 4.}),
                new ArrayRealVector(new double[]{5., 6.}));
        final List<RealMatrix> precisionsChol = Arrays.asList(  // using expected from testComputePrecisionCholesky
                new Array2DRowRealMatrix(new double[][]{
                        {0.70710678, -0.40824829},
                        {0., 0.81649658}}),
                new Array2DRowRealMatrix(new double[][]{
                        {0.5, -0.56694671},
                        {0., 0.75592895}}),
                new Array2DRowRealMatrix(new double[][]{
                        {0.40824829, -0.61545745},
                        {0., 0.73854895}}));
        final RealMatrix result = BayesianGaussianMixtureUtils.estimateLogGaussianProb(X, means, precisionsChol);
        final RealMatrix expected = new Array2DRowRealMatrix(new double[][]{
                {-2.38718321, -3.38226071, -4.49137016},
                {-2.72051654, -2.95368928, -3.85500652},
                {-3.72051654, -2.81083214, -3.40046107},
                {-5.38718321, -2.95368928, -3.12773379},
                {-7.72051654, -3.38226071, -3.0368247}});
        Assert.assertTrue(result.subtract(expected).getNorm() < epsilon);
    }

    @Test
    public void testEstimateGaussianParameters() {
        final double epsilon = 1E-5;
        final double regCovar = 1E-10;
        // nSamples = 5, nFeatures = 2, nComponents = 3
        final RealMatrix X = new Array2DRowRealMatrix(new double[][]{       // scale second feature so that covariances are interesting
                {1., 20.},
                {2., 30.},
                {3., 40.},
                {4., 50.},
                {5., 60.}});
        final RealMatrix resp = new Array2DRowRealMatrix(new double[][]{    // randomly generated: np.random.seed(0); sklearn.preprocessing.normalize(np.random.random((5, 3)), axis=1, norm='l1')
                {0.29399155, 0.38311672, 0.32289173},
                {0.33750765, 0.26241723, 0.40007512},
                {0.1908342,  0.38890714, 0.42025866},
                {0.22501625, 0.46461061, 0.31037314},
                {0.36304264, 0.59155755, 0.04539982}});
        final Triple<RealVector, List<RealVector>, List<RealMatrix>> result = BayesianGaussianMixtureUtils.estimateGaussianParameters(X, resp, regCovar, EPSILON);
        final Triple<RealVector, List<RealVector>, List<RealMatrix>> expected = Triple.of(
                new ArrayRealVector(new double[]{1.41039229, 2.09060924, 1.49899847}),
                Arrays.asList(
                        new ArrayRealVector(new double[]{3.01815862, 40.18158617}),
                        new ArrayRealVector(new double[]{3.29612182, 42.96121824}),
                        new ArrayRealVector(new double[]{2.56992231, 35.6992231})),
                Arrays.asList(
                        new Array2DRowRealMatrix(new double[][]{
                                {2.26192176, 22.6192076},
                                {22.6192076, 226.19207698}}),
                        new Array2DRowRealMatrix(new double[][]{
                                {2.12493338, 21.24932385},
                                {21.24932385, 212.49323947}}),
                        new Array2DRowRealMatrix(new double[][]{
                                {1.27174977, 12.71748769},
                                {12.71748769, 127.17487785}})));
        Assert.assertEquals(result.getLeft().toArray(), expected.getLeft().toArray(), epsilon);    // effective number
        IntStream.range(0, result.getMiddle().size()).forEach(k ->                              // means
                Assert.assertEquals(result.getMiddle().get(k).toArray(), expected.getMiddle().get(k).toArray(), epsilon));
        IntStream.range(0, result.getRight().size()).forEach(k ->                               // covariances
                Assert.assertTrue(result.getRight().get(k).subtract(expected.getRight().get(k)).getNorm() < epsilon));
    }
}