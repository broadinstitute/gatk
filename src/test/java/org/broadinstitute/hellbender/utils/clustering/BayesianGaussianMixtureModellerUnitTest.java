package org.broadinstitute.hellbender.utils.clustering;

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

public final class BayesianGaussianMixtureModellerUnitTest extends GATKBaseTest {

    private static final File TEST_SUB_DIR = new File(packageRootTestDir, "utils/clustering");

    private static final double ALLOWED_DELTA = 1E-6;

    /**
     * Simulated data was generated using the following Python code:
     *
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

    /**
     * Checks against expected results generated using the following Python code:
     *
     * import numpy as np # 1.23.1
     * from sklearn.mixture import BayesianGaussianMixture # 1.0.2
     *
     * class BayesianGaussianMixtureTestInit(BayesianGaussianMixture):
     *     def _initialize_parameters(self, X, random_state):
     *         n_samples, _ = X.shape
     *         resp = np.repeat(np.arange(self.n_components, dtype=np.float)[np.newaxis, :], n_samples, axis=0)
     *         resp /= resp.sum(axis=1)[:, np.newaxis]
     *         self._initialize(X, resp)
     *
     * n_components = 4
     * n_features = 3
     * X_ni = np.loadtxt('bayesian-gaussian-mixture-simulated-data-10k-samples-4-components-3-features.tsv')
     * bgmm = BayesianGaussianMixtureTestInit(n_components=n_components,
     *                                        tol=1E-3,
     *                                        reg_covar=1E-6,
     *                                        max_iter=100,
     *                                        n_init=1,
     *                                        weight_concentration_prior_type='dirichlet_distribution',
     *                                        weight_concentration_prior=1E-2,
     *                                        mean_precision_prior=10.,
     *                                        mean_prior=np.zeros(n_features),
     *                                        degrees_of_freedom_prior=n_features,
     *                                        covariance_type='full',
     *                                        covariance_prior=np.eye(n_features),
     *                                        random_state=1,
     *                                        warm_start=False)
     * bgmm.fit(X_ni)
     *
     * np.set_printoptions(suppress=True)
     * print('lower_bound_:\n', bgmm.lower_bound_)
     * print('weight_concentration_:\n', bgmm.weight_concentration_)
     * print('weights_:\n', bgmm.weights_)
     * print('mean_precision_:\n', bgmm.mean_precision_)
     * print('means_:\n', bgmm.means_)
     * print('precisions_cholesky_:\n', bgmm.precisions_cholesky_)
     * print('covariances_:\n', bgmm.covariances_)
     * print('degrees_of_freedom_:\n', bgmm.degrees_of_freedom_)
     */
    @Test
    public void testSimulatedData() throws IOException {
        final int nComponents = 4;
        final int nFeatures = 3;
        final int nSubsample = 1000;

        final double[][] data = readData(SIMULATED_DATA_FILE);

        final double[] meanPrior = new double[nFeatures];
        Arrays.fill(meanPrior, 0.);
        final double[][] covariancePrior = MatrixUtils.createRealIdentityMatrix(nFeatures).getData();

        final BayesianGaussianMixtureModeller bgmm = new BayesianGaussianMixtureModeller.Builder()
                .nComponents(nComponents)
                .tol(1E-3)
                .regCovar(1E-6)
                .maxIter(100)
                .nInit(1)
                .initMethod(BayesianGaussianMixtureModeller.InitMethod.TEST) // responsibilities proportional to component index
                .weightConcentrationPrior(1E-2)
                .meanPrecisionPrior(10.)
                .meanPrior(meanPrior)
                .degreesOfFreedomPrior((double) nFeatures)
                .covariancePrior(covariancePrior)
                .seed(1)
                .warmStart(false)
                .verboseInterval(1)
                .relativeSymmetryThreshold(1E-6)
                .absolutePositivityThreshold(1E-10)
                .epsilon(1E-10)
                .build();

        final double[][] dataSubsample = IntStream.range(0, nSubsample).boxed().map(i -> data[i]).toArray(double[][]::new);
        bgmm.fit(dataSubsample); // warmStart = false, so this fit doesn't matter
        bgmm.fit(data);

        final BayesianGaussianMixtureModelPosterior fit = bgmm.getBestFit();

        final double expectedLowerBound = -55411.488111726714;
        final RealVector expectedWeightConcentration = new ArrayRealVector(new double[]
                {0.01, 4489.975392, 1194.846932, 4315.207676});
        final RealVector expectedWeights = new ArrayRealVector(new double[]
                {0.000001, 0.448996, 0.119484, 0.431519});
        final RealVector expectedMeanPrecision = new ArrayRealVector(new double[]
                {10., 4499.965392, 1204.836932, 4325.197676});
        final List<RealVector> expectedMeans = Arrays.asList(
                new ArrayRealVector(new double[]{0., 0., -0.}),
                new ArrayRealVector(new double[]{-2.279021, -7.080966, -6.519024}),
                new ArrayRealVector(new double[]{4.166038, 2.445518, -4.251733}),
                new ArrayRealVector(new double[]{-10.270325, -5.136062, -8.918113}));
        final List<RealMatrix> expectedPrecisionsCholesky = Arrays.asList(
                new Array2DRowRealMatrix(new double[][]{
                        {1.732051, -0., 0.},
                        {0., 1.732051, 0.},
                        {0., 0., 1.732051}}),
                new Array2DRowRealMatrix(new double[][]{
                        {0.396528, -0.697808, -0.856619},
                        {0., 0.300155, 0.969691},
                        {0., 0., 1.155904}}),
                new Array2DRowRealMatrix(new double[][]{
                        {0.736986, -0.202604, 0.305607},
                        {0., 0.227671, -0.216366},
                        {0., 0., 1.014755}}),
                new Array2DRowRealMatrix(new double[][]{
                        {0.484877, -0.197117, -0.033798},
                        {0., 0.033003, 0.018767},
                        {0., 0., 0.66848}}));
        final List<RealMatrix> expectedCovariances = Arrays.asList(
                new Array2DRowRealMatrix(new double[][]{
                        {0.333333, 0., -0.},
                        {0., 0.333333, -0.},
                        {-0., -0., 0.333333}}),
                new Array2DRowRealMatrix(new double[][]{
                        {6.359926, 14.785738, -7.690568},
                        {14.785738, 45.473974, -27.19079},
                        {-7.690568, -27.19079, 17.859521}}),
                new Array2DRowRealMatrix(new double[][]{
                        {1.84112, 1.638408, -0.205136},
                        {1.638408, 20.750324, 3.930963},
                        {-0.205136, 3.930963, 1.871072}}),
                new Array2DRowRealMatrix(new double[][]{
                        {4.253403, 25.404238, -0.498157},
                        {25.404238, 1069.833413, -28.750485},
                        {-0.498157, -28.750485, 3.019773}}));
        final RealVector expectedDegreesOfFreedom = new ArrayRealVector(new double[]
                {3., 4492.965392, 1197.836932, 4318.197676});

        Assert.assertEquals(bgmm.getLowerBound(), expectedLowerBound, ALLOWED_DELTA);
        assertEqualUpToAllowedDelta(fit.getWeightConcentration(), expectedWeightConcentration, ALLOWED_DELTA);
        assertEqualUpToAllowedDelta(fit.getWeights(), expectedWeights, ALLOWED_DELTA);
        assertEqualUpToAllowedDelta(fit.getMeanPrecision(), expectedMeanPrecision, ALLOWED_DELTA);
        assertListRealVectorsEqualUpToAllowedDelta(fit.getMeans(), expectedMeans, ALLOWED_DELTA);
        assertListRealMatricesEqualUpToAllowedDelta(fit.getPrecisionsCholesky(), expectedPrecisionsCholesky, ALLOWED_DELTA);
        assertListRealMatricesEqualUpToAllowedDelta(fit.getCovariances(), expectedCovariances, ALLOWED_DELTA);
        assertEqualUpToAllowedDelta(fit.getDegreesOfFreedom(), expectedDegreesOfFreedom, ALLOWED_DELTA);
    }

    /**
     * Checks against expected results generated using the following Python code:
     *
     * import numpy as np # 1.23.1
     * from sklearn.mixture import BayesianGaussianMixture # 1.0.2
     *
     * class BayesianGaussianMixtureTestInit(BayesianGaussianMixture):
     *     def _initialize_parameters(self, X, random_state):
     *         n_samples, _ = X.shape
     *         resp = np.repeat(np.arange(self.n_components, dtype=np.float)[np.newaxis, :], n_samples, axis=0)
     *         resp /= resp.sum(axis=1)[:, np.newaxis]
     *         self._initialize(X, resp)
     *
     * n_components = 4
     * n_features = 3
     * n_subsample = 1000
     * X_ni = np.loadtxt('bayesian-gaussian-mixture-simulated-data-10k-samples-4-components-3-features.tsv')
     * bgmm = BayesianGaussianMixtureTestInit(n_components=n_components,
     *                                        tol=1E-3,
     *                                        reg_covar=1E-6,
     *                                        max_iter=100,
     *                                        n_init=1,
     *                                        weight_concentration_prior_type='dirichlet_distribution',
     *                                        weight_concentration_prior=1E-2,
     *                                        mean_precision_prior=10.,
     *                                        mean_prior=np.zeros(n_features),
     *                                        degrees_of_freedom_prior=n_features,
     *                                        covariance_type='full',
     *                                        covariance_prior=np.eye(n_features),
     *                                        random_state=1,
     *                                        warm_start=True)
     * bgmm.fit(X_ni[:n_subsample])
     * bgmm.fit(X_ni)
     *
     * np.set_printoptions(suppress=True)
     * print('lower_bound_:\n', bgmm.lower_bound_)
     * print('weight_concentration_:\n', bgmm.weight_concentration_)
     * print('weights_:\n', bgmm.weights_)
     * print('mean_precision_:\n', bgmm.mean_precision_)
     * print('means_:\n', bgmm.means_)
     * print('precisions_cholesky_:\n', bgmm.precisions_cholesky_)
     * print('covariances_:\n', bgmm.covariances_)
     * print('degrees_of_freedom_:\n', bgmm.degrees_of_freedom_)
     */
    @Test
    public void testSimulatedDataWithWarmStart() throws IOException {
        final int nComponents = 4;
        final int nFeatures = 3;
        final int nSubsample = 1000;

        final double[][] data = readData(SIMULATED_DATA_FILE);

        final double[] meanPrior = new double[nFeatures];
        Arrays.fill(meanPrior, 0.);
        final double[][] covariancePrior = MatrixUtils.createRealIdentityMatrix(nFeatures).getData();

        final BayesianGaussianMixtureModeller bgmm = new BayesianGaussianMixtureModeller.Builder()
                .nComponents(nComponents)
                .tol(1E-3)
                .regCovar(1E-6)
                .maxIter(100)
                .nInit(1)
                .initMethod(BayesianGaussianMixtureModeller.InitMethod.TEST) // responsibilities proportional to component index
                .weightConcentrationPrior(1E-2)
                .meanPrecisionPrior(10.)
                .meanPrior(meanPrior)
                .degreesOfFreedomPrior((double) nFeatures)
                .covariancePrior(covariancePrior)
                .seed(1)
                .warmStart(true)
                .verboseInterval(1)
                .relativeSymmetryThreshold(1E-6)
                .absolutePositivityThreshold(1E-10)
                .epsilon(1E-10)
                .build();

        final double[][] dataSubsample = IntStream.range(0, nSubsample).boxed().map(i -> data[i]).toArray(double[][]::new);
        bgmm.fit(dataSubsample);
        bgmm.fit(data);

        final BayesianGaussianMixtureModelPosterior fit = bgmm.getBestFit();

        final double expectedLowerBound =  -55411.488032617875;
        final RealVector expectedWeightConcentration = new ArrayRealVector(new double[]
                {0.01, 4490.195421, 1194.635073, 4315.199506});
        final RealVector expectedWeights = new ArrayRealVector(new double[]
                {0.000001, 0.449018, 0.119463, 0.431518});
        final RealVector expectedMeanPrecision = new ArrayRealVector(new double[]
                {10., 4500.185421, 1204.625073, 4325.189506});
        final List<RealVector> expectedMeans = Arrays.asList(
                new ArrayRealVector(new double[]{0., 0., -0.}),
                new ArrayRealVector(new double[]{-2.278808, -7.080753, -6.519043}),
                new ArrayRealVector(new double[]{4.166343, 2.446393, -4.251281}),
                new ArrayRealVector(new double[]{-10.270331, -5.136057, -8.918112}));
        final List<RealMatrix> expectedPrecisionsCholesky = Arrays.asList(
                new Array2DRowRealMatrix(new double[][]{
                        {1.732051, -0., 0.},
                        {0., 1.732051, 0.},
                        {0., 0., 1.732051}}),
                new Array2DRowRealMatrix(new double[][]{
                        {0.396503, -0.697686, -0.856625},
                        {0., 0.300133, 0.969693},
                        {0., 0., 1.155907}}),
                new Array2DRowRealMatrix(new double[][]{
                        {0.737046, -0.202458, 0.305988},
                        {0., 0.227665, -0.216383},
                        {0., 0., 1.015106}}),
                new Array2DRowRealMatrix(new double[][]{
                        {0.484878, -0.197118, -0.0338},
                        {0., 0.033003, 0.018767},
                        {0., 0., 0.668482}}));
        final List<RealMatrix> expectedCovariances = Arrays.asList(
                new Array2DRowRealMatrix(new double[][]{
                        {0.333333, 0., -0.},
                        {0., 0.333333, -0.},
                        {-0., -0., 0.333333}}),
                new Array2DRowRealMatrix(new double[][]{
                        {6.360724, 14.786082, -7.690245},
                        {14.786082, 45.472877, -27.189564},
                        {-7.690245, -27.189564, 17.858701}}),
                new Array2DRowRealMatrix(new double[][]{
                        {1.840817, 1.637007, -0.205938},
                        {1.637007, 20.749173, 3.929499},
                        {-0.205938, 3.929499, 1.870158}}),
                new Array2DRowRealMatrix(new double[][]{
                        {4.253385, 25.404296, -0.49815},
                        {25.404296, 1069.835392, -28.750524},
                        {-0.49815, -28.750524, 3.019765}}));
        final RealVector expectedDegreesOfFreedom = new ArrayRealVector(new double[]
                {3., 4493.185421, 1197.625073, 4318.189506});

        Assert.assertEquals(bgmm.getLowerBound(), expectedLowerBound, ALLOWED_DELTA);
        assertEqualUpToAllowedDelta(fit.getWeightConcentration(), expectedWeightConcentration, ALLOWED_DELTA);
        assertEqualUpToAllowedDelta(fit.getWeights(), expectedWeights, ALLOWED_DELTA);
        assertEqualUpToAllowedDelta(fit.getMeanPrecision(), expectedMeanPrecision, ALLOWED_DELTA);
        assertListRealVectorsEqualUpToAllowedDelta(fit.getMeans(), expectedMeans, ALLOWED_DELTA);
        assertListRealMatricesEqualUpToAllowedDelta(fit.getPrecisionsCholesky(), expectedPrecisionsCholesky, ALLOWED_DELTA);
        assertListRealMatricesEqualUpToAllowedDelta(fit.getCovariances(), expectedCovariances, ALLOWED_DELTA);
        assertEqualUpToAllowedDelta(fit.getDegreesOfFreedom(), expectedDegreesOfFreedom, ALLOWED_DELTA);
    }

    /**
     * Generated with Python code identical to that used in {@link #testSimulatedData}, but with n_components = 10.
     */
    @Test
    public void testSimulatedDataWithManyExtraComponents() throws IOException {
        final int nComponents = 10;
        final int nFeatures = 3;

        final double[][] data = readData(SIMULATED_DATA_FILE);

        final double[] meanPrior = new double[nFeatures];
        Arrays.fill(meanPrior, 0.);
        final double[][] covariancePrior = MatrixUtils.createRealIdentityMatrix(nFeatures).getData();

        final BayesianGaussianMixtureModeller bgmm = new BayesianGaussianMixtureModeller.Builder()
                .nComponents(nComponents)
                .tol(1E-3)
                .regCovar(1E-6)
                .maxIter(100)
                .nInit(1)
                .initMethod(BayesianGaussianMixtureModeller.InitMethod.TEST) // responsibilities proportional to component index
                .weightConcentrationPrior(1E-2)
                .meanPrecisionPrior(10.)
                .meanPrior(meanPrior)
                .degreesOfFreedomPrior((double) nFeatures)
                .covariancePrior(covariancePrior)
                .seed(1)
                .warmStart(false)
                .verboseInterval(1)
                .relativeSymmetryThreshold(1E-6)
                .absolutePositivityThreshold(1E-10)
                .epsilon(1E-10)
                .build();

        bgmm.fit(data);

        final BayesianGaussianMixtureModelPosterior fit = bgmm.getBestFit();

        final double expectedLowerBound = -48059.319701346176;
        final RealVector expectedWeightConcentration = new ArrayRealVector(new double[]
                {0.01, 1666.009231, 2880.492916, 1133.141938, 0.01, 0.01, 0.01, 0.01, 0.01, 4320.395915});
        final RealVector expectedWeights = new ArrayRealVector(new double[]
                {0.000001, 0.166599, 0.288046, 0.113313, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.432035});
        final RealVector expectedMeanPrecision = new ArrayRealVector(new double[]
                {10., 1675.999231, 2890.482916, 1143.131938, 10., 10., 10., 10., 10., 4330.385915});
        final List<RealVector> expectedMeans = Arrays.asList(
                new ArrayRealVector(new double[]{0., 0., -0.}),
                new ArrayRealVector(new double[]{-5.087788, -15.640294, -1.235359}),
                new ArrayRealVector(new double[]{-0.517346, -1.982549, -9.556166}),
                new ArrayRealVector(new double[]{4.217605, 2.696524, -4.138238}),
                new ArrayRealVector(new double[]{0., 0., -0.}),
                new ArrayRealVector(new double[]{0., 0., -0.}),
                new ArrayRealVector(new double[]{0., 0., -0.}),
                new ArrayRealVector(new double[]{0., 0., -0.}),
                new ArrayRealVector(new double[]{0., 0., -0.}),
                new ArrayRealVector(new double[]{-10.266077, -5.142953, -8.915529}));
        final List<RealMatrix> expectedPrecisionsCholesky = Arrays.asList(
                new Array2DRowRealMatrix(new double[][]{
                        {1.732051, -0., 0.},
                        {0., 1.732051, 0.},
                        {0., 0., 1.732051}}),
                new Array2DRowRealMatrix(new double[][]{
                        {1.660752, -1.497842, -0.835875},
                        {0., 1.017818, -0.250983},
                        {0., 0., 2.57486}}),
                new Array2DRowRealMatrix(new double[][]{
                        {0.61583, -0.098797, -1.169101},
                        {0., 0.736576, 0.669181},
                        {0., 0., 1.434016}}),
                new Array2DRowRealMatrix(new double[][]{
                        {0.731082, -0.177327, 0.370228},
                        {0., 0.226852, -0.213995},
                        {0., 0., 1.07647}}),
                new Array2DRowRealMatrix(new double[][]{
                        {1.732051, -0., 0.},
                        {0., 1.732051, 0.},
                        {0., 0., 1.732051}}),
                new Array2DRowRealMatrix(new double[][]{
                        {1.732051, -0., 0.},
                        {0., 1.732051, 0.},
                        {0., 0., 1.732051}}),
                new Array2DRowRealMatrix(new double[][]{
                        {1.732051, -0., 0.},
                        {0., 1.732051, 0.},
                        {0., 0., 1.732051}}),
                new Array2DRowRealMatrix(new double[][]{
                        {1.732051, -0., 0.},
                        {0., 1.732051, 0.},
                        {0., 0., 1.732051}}),
                new Array2DRowRealMatrix(new double[][]{
                        {1.732051, -0., 0.},
                        {0., 1.732051, 0.},
                        {0., 0., 1.732051}}),
                new Array2DRowRealMatrix(new double[][]{
                        { 0.484257, -0.196219, -0.035342},
                        {0., 0.033007, 0.018823},
                        {0., 0., 0.669047}}));
        final List<RealMatrix> expectedCovariances = Arrays.asList(
                new Array2DRowRealMatrix(new double[][]{
                        {0.333333, 0., -0.},
                        {0., 0.333333, -0.},
                        {-0., -0., 0.333333}}),
                new Array2DRowRealMatrix(new double[][]{
                        {0.362569, 0.533563, 0.169709},
                        {0.533563, 1.750496, 0.343839},
                        {0.169709, 0.343839, 0.23944}}),
                new Array2DRowRealMatrix(new double[][]{
                        {2.63681, 0.353674, 1.984654},
                        {0.353674, 1.890604, -0.593909},
                        {1.984654, -0.593909, 2.38145}}),
                new Array2DRowRealMatrix(new double[][]{
                        {1.870972, 1.462518, -0.35274},
                        {1.462518, 20.575177, 3.587216},
                        {-0.35274, 3.587216, 1.697404}}),
                new Array2DRowRealMatrix(new double[][]{
                        {0.333333, 0., -0.},
                        {0., 0.333333, -0.},
                        {-0., -0., 0.333333}}),
                new Array2DRowRealMatrix(new double[][]{
                        {0.333333, 0., -0.},
                        {0., 0.333333, -0.},
                        {-0., -0., 0.333333}}),
                new Array2DRowRealMatrix(new double[][]{
                        {0.333333, 0., -0.},
                        {0., 0.333333, -0.},
                        {-0., -0., 0.333333}}),
                new Array2DRowRealMatrix(new double[][]{
                        {0.333333, 0., -0.},
                        {0., 0.333333, -0.},
                        {-0., -0., 0.333333}}),
                new Array2DRowRealMatrix(new double[][]{
                        {0.333333, 0., -0.},
                        {0., 0.333333, -0.},
                        {-0., -0., 0.333333}}),
                new Array2DRowRealMatrix(new double[][]{
                        {4.26431, 25.350321, -0.487939},
                        {25.350321, 1068.585872, -28.724042},
                        {-0.487939, -28.724042, 3.016354}}));
        final RealVector expectedDegreesOfFreedom = new ArrayRealVector(new double[]
                {3., 1668.999231, 2883.482916, 1136.131938, 3., 3., 3., 3., 3., 4323.385915});

        Assert.assertEquals(bgmm.getLowerBound(), expectedLowerBound, ALLOWED_DELTA);
        assertEqualUpToAllowedDelta(fit.getWeightConcentration(), expectedWeightConcentration, ALLOWED_DELTA);
        assertEqualUpToAllowedDelta(fit.getWeights(), expectedWeights, ALLOWED_DELTA);
        assertEqualUpToAllowedDelta(fit.getMeanPrecision(), expectedMeanPrecision, ALLOWED_DELTA);
        assertListRealVectorsEqualUpToAllowedDelta(fit.getMeans(), expectedMeans, ALLOWED_DELTA);
        assertListRealMatricesEqualUpToAllowedDelta(fit.getPrecisionsCholesky(), expectedPrecisionsCholesky, ALLOWED_DELTA);
        assertListRealMatricesEqualUpToAllowedDelta(fit.getCovariances(), expectedCovariances, ALLOWED_DELTA);
        assertEqualUpToAllowedDelta(fit.getDegreesOfFreedom(), expectedDegreesOfFreedom, ALLOWED_DELTA);
    }

    private static void assertEqualUpToAllowedDelta(final RealVector a,
                                                    final RealVector b,
                                                    final double delta) {
        Assert.assertEquals(a.getDimension(), b.getDimension());
        IntStream.range(0, a.getDimension()).forEach(i -> Assert.assertEquals(a.getEntry(i), b.getEntry(i), delta));
    }

    private static void assertListRealVectorsEqualUpToAllowedDelta(final List<RealVector> a,
                                                                   final List<RealVector> b,
                                                                   final double delta) {
        Assert.assertEquals(a.size(), b.size());
        IntStream.range(0, a.size()).forEach(i -> assertEqualUpToAllowedDelta(a.get(i), b.get(i), delta));
    }

    private static void assertEqualUpToAllowedDelta(final RealMatrix a,
                                                    final RealMatrix b,
                                                    final double delta) {
        Assert.assertEquals(a.getRowDimension(), b.getRowDimension());
        Assert.assertEquals(a.getColumnDimension(), b.getColumnDimension());
        IntStream.range(0, a.getRowDimension()).forEach(
                i -> assertEqualUpToAllowedDelta(a.getRowVector(i), b.getRowVector(i), delta));
    }

    private static void assertListRealMatricesEqualUpToAllowedDelta(final List<RealMatrix> a,
                                                                    final List<RealMatrix> b,
                                                                    final double delta) {
        Assert.assertEquals(a.size(), b.size());
        IntStream.range(0, a.size()).forEach(i -> assertEqualUpToAllowedDelta(a.get(i), b.get(i), delta));
    }
}