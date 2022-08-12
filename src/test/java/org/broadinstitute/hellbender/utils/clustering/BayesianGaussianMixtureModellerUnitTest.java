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
     *                                        random_state=1)
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
                .initMethod(BayesianGaussianMixtureModeller.InitMethod.TEST)
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
                        {14.785738, 45.473974, -27.19079 },
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

    @Test
    public void testSimulatedDataWithWarmStart() throws IOException {

    }

    @Test
    public void testSimulatedDataWithManyExtraComponents() throws IOException {

    }

    @Test
    public void testSimulatedDataWithMultipleInitializations() throws IOException {

    }
}