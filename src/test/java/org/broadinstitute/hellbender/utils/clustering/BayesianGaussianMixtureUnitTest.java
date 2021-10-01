package org.broadinstitute.hellbender.utils.clustering;

import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;

public final class BayesianGaussianMixtureUnitTest extends GATKBaseTest {
    private static final File TEST_SUB_DIR = new File(packageRootTestDir, "utils/clustering");

    /**
     * import numpy as np  # np.__version__    = 1.19.5
     * import torch        # torch.__version__ = 1.8.0
     * import pyro         # pyro.__version__  = 1.6.0
     * import pyro.distributions as dist
     *
     * def generate_gmm_data(seed=1,
     *                       num_points=10000,
     *                       num_components=3,
     *                       num_features=5,
     *                       weight_concentration=1.,
     *                       mean_scale=10.,
     *                       variance_scale=1.,
     *                       lkj_concentration=1. # uniform distribution over correlation matrices, see https://pyro.ai/examples/lkj.html
     *                       ):
     *     pyro.set_rng_seed(seed)
     *     pi_k = pyro.sample('weight_k', dist.Dirichlet(weight_concentration * torch.ones(num_components)))
     *
     *     with pyro.plate('components', size=num_components, dim=-2):
     *         mu_ki = pyro.sample('mean_ki', dist.Normal(0., mean_scale * torch.ones(num_features)))
     *         theta_ki = pyro.sample('variance_ki', dist.HalfCauchy(variance_scale * torch.ones(num_features)))
     *         L_corr_kij = pyro.sample('lower_cholesky_correlation_kij', dist.LKJCholesky(num_features, lkj_concentration))
     *
     *     L_cov_kij = torch.bmm(theta_ki.sqrt().diag_embed(), L_corr_kij.squeeze(dim=-3))
     *
     *     with pyro.plate('data', size=num_points):
     *         z_n = pyro.sample('assignment', dist.Categorical(pi_k))
     *         X_ni = pyro.sample('obs_n', dist.MultivariateNormal(mu_ki[z_n], scale_tril=L_cov_kij[z_n]))
     *
     *     return (x.detach().numpy() for x in [pi_k, mu_ki, L_cov_kij, z_n, X_ni])
     *
     * pi_k, mu_ki, L_cov_kij, z_n, X_ni = generate_gmm_data()
     * cov_kij = np.einsum('kij,klj->kil', L_cov_kij, L_cov_kij)
     * np.savetxt('bayesian-gaussian-mixture-simulated-data-10k-samples-4-components-3-features.tsv', X_ni, delimiter='\t')
     */
    private static final File SIMULATED_DATA_FILE = new File(TEST_SUB_DIR, "bayesian-gaussian-mixture-simulated-data-10k-samples-4-components-3-features.tsv");

    private static double[][] readData(final File file) throws IOException {
        final List<String> lines = Files.readAllLines(file.toPath(), StandardCharsets.UTF_8);
        final double[][] data = new double[lines.size()][];
        for (int i = 0; i < lines.size(); i++) {
            data[i] = Arrays.stream(lines.get(i).split("\t")).mapToDouble(Double::parseDouble).toArray();
        }
        return data;
    }

    @Test
    public void testSimulatedData() throws IOException {
        final int nComponents = 4;
        final int nFeatures = 3;

        final double[][] data = readData(SIMULATED_DATA_FILE);

        final double[] meanPrior = new double[nFeatures];
        Arrays.fill(meanPrior, 0.);
        final double[][] covariancePrior = MatrixUtils.createRealIdentityMatrix(nFeatures).getData();

        final BayesianGaussianMixture bgmm = new BayesianGaussianMixture.Builder()
                .nComponents(nComponents)
                .tol(1E-3)
                .regCovar(1E-6)
                .maxIter(100)
                .nInit(1)
                .initMethod(BayesianGaussianMixture.InitMethod.TEST)
                .weightConcentrationPrior(1E-2)
                .meanPrecisionPrior(10.)
                .meanPrior(meanPrior)
                .degreesOfFreedomPrior(nFeatures)
                .covariancePrior(covariancePrior)
                .seed(1)
                .warmStart(true)
                .verboseInterval(1)
                .build();

        bgmm.fit(data);
        System.out.println("weights: " + bgmm.getWeights());
        System.out.println("meanPrecision: " + bgmm.getMeanPrecision());
        System.out.println("means: " + bgmm.getMeans());
        System.out.println("precisionsCholesky: " + bgmm.getPrecisionsCholesky());
        System.out.println("covariances: " + bgmm.getCovariances());
        System.out.println("degreesOfFreedom: " + bgmm.getDegreesOfFreedom());
    }

//    @Test
//    public void testSimulatedData1Mx3x10() throws IOException {
//        final int nComponents = 10;
//        final int nFeatures = 10;
//
//        final double[][] data = readData(new File("/home/slee/working/malariagen/issues/hyperhet/simulated-data.tsv"));
//
//        final double[] meanPrior = new double[nFeatures];
//        Arrays.fill(meanPrior, 0.);
//        final double[][] covariancePrior = MatrixUtils.createRealIdentityMatrix(nFeatures).getData();
//
//        final BayesianGaussianMixture bgmm = new BayesianGaussianMixture.Builder()
//                .nComponents(nComponents)
//                .tol(1E-3)
//                .regCovar(1E-6)
//                .maxIter(100)
//                .nInit(1)
//                .initMethod(BayesianGaussianMixture.InitMethod.TEST)
//                .weightConcentrationPrior(1E-2)
//                .meanPrecisionPrior(10.)
//                .meanPrior(meanPrior)
//                .degreesOfFreedomPrior(nFeatures)
//                .covariancePrior(covariancePrior)
//                .seed(1)
//                .warmStart(true)
//                .verboseInterval(1)
//                .build();
//
//        bgmm.fit(Arrays.copyOfRange(data, 0, 100000));
//        System.out.println("weights: " + bgmm.getWeights());
//        System.out.println("meanPrecision: " + bgmm.getMeanPrecision());
//        System.out.println("means: " + bgmm.getMeans());
//        System.out.println("precisionsCholesky: " + bgmm.getPrecisionsCholesky());
//        System.out.println("covariances: " + bgmm.getCovariances());
//        System.out.println("degreesOfFreedom: " + bgmm.getDegreesOfFreedom());
//
//        bgmm.fit(data);
//        System.out.println("weights: " + bgmm.getWeights());
//        System.out.println("meanPrecision: " + bgmm.getMeanPrecision());
//        System.out.println("means: " + bgmm.getMeans());
//        System.out.println("precisionsCholesky: " + bgmm.getPrecisionsCholesky());
//        System.out.println("covariances: " + bgmm.getCovariances());
//        System.out.println("degreesOfFreedom: " + bgmm.getDegreesOfFreedom());
//    }

    @Test
    public void testSimulatedDataWithWarmStart() throws IOException {

    }

    @Test
    public void testSimulatedDataWithManyExtraComponents() throws IOException {

    }

    @Test
    public void testSimulatedDataWithMultipleInitializations() throws IOException {

    }

    //******************************************************************************************************************
    // unit tests for static helper methods; checked against python implementations
    //******************************************************************************************************************

    @Test
    public void testLogDirichletNorm() {
        final double epsilon = 1E-10;
        Assert.assertEquals(
                BayesianGaussianMixture.logDirichletNorm(new ArrayRealVector(new double[]{1., 2., 3., 4.})),
                10.31692083029347, epsilon);
        Assert.assertEquals(
                BayesianGaussianMixture.logDirichletNorm(new ArrayRealVector(new double[]{0.1, 0.2, 0.3, 0.4})),
                -5.669252286684849, epsilon);
    }

    @Test
    public void testLogWishartNorm() {
        final double epsilon = 1E-6;
        Assert.assertEquals(
                BayesianGaussianMixture.logWishartNorm(
                        new ArrayRealVector(new double[]{3., 4., 5., 6.}),
                        new ArrayRealVector(new double[]{0.1, 0.2, 0.3, 0.4}),
                2).toArray(), // nFeatures must be <= all elements of degreesOfFreedom
                new double[]{-2.2586593 , -3.45180648, -5.25041877, -7.53671313}, epsilon);
    }

    @Test
    public void testComputePrecisionCholesky() {
        final double epsilon = 1E-6;
        final List<RealMatrix> result = BayesianGaussianMixture.computePrecisionCholesky(
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
        final RealVector result = BayesianGaussianMixture.computeLogDetCholesky(
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
        final RealMatrix result = BayesianGaussianMixture.estimateLogGaussianProb(X, means, precisionsChol);
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
        final Triple<RealVector, List<RealVector>, List<RealMatrix>> result = BayesianGaussianMixture.estimateGaussianParameters(X, resp, regCovar);
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