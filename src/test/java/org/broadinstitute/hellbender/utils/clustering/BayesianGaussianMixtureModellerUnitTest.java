package org.broadinstitute.hellbender.utils.clustering;

import org.apache.commons.math3.linear.MatrixUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;

public final class BayesianGaussianMixtureModellerUnitTest extends GATKBaseTest {
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
        System.out.println(bgmm);
        System.out.println(bgmm.getLowerBound());
        final BayesianGaussianMixtureModelPosterior fit = bgmm.getBestFit();
        System.out.println("weights: " + fit.getWeights());
        System.out.println("meanPrecision: " + fit.getMeanPrecision());
        System.out.println("means: " + fit.getMeans());
        System.out.println("precisionsCholesky: " + fit.getPrecisionsCholesky());
        System.out.println("covariances: " + fit.getCovariances());
        System.out.println("degreesOfFreedom: " + fit.getDegreesOfFreedom());

        //ll -55411.48811
        //weights__:
        // [0.000001   0.44899574 0.11948422 0.43151904]
        //mean_precision_:
        // [  10.         4499.96539243 1204.83693171 4325.19767586]
        //means_:
        // [[  0.           0.          -0.        ]
        // [ -2.27902144  -7.08096566  -6.51902374]
        // [  4.16603771   2.44551832  -4.2517331 ]
        // [-10.27032453  -5.13606169  -8.91811286]]
        //precisions_cholesky_:
        // [[[ 1.73205081 -0.          0.        ]
        //  [ 0.          1.73205081  0.        ]
        //  [ 0.          0.          1.73205081]]
        //
        // [[ 0.3965281  -0.69780817 -0.85661924]
        //  [ 0.          0.30015465  0.9696907 ]
        //  [ 0.          0.          1.15590395]]
        //
        // [[ 0.73698552 -0.20260387  0.30560695]
        //  [ 0.          0.22767109 -0.21636638]
        //  [ 0.          0.          1.01475474]]
        //
        // [[ 0.48487715 -0.19711704 -0.03379824]
        //  [ 0.          0.03300309  0.01876718]
        //  [ 0.          0.          0.6684803 ]]]
        //covariances_:
        // [[[   0.33333333    0.           -0.        ]
        //  [   0.            0.33333333   -0.        ]
        //  [  -0.           -0.            0.33333333]]
        //
        // [[   6.35992584   14.78573847   -7.69056826]
        //  [  14.78573847   45.47397408  -27.19079008]
        //  [  -7.69056826  -27.19079008   17.85952135]]
        //
        // [[   1.84111996    1.63840757   -0.20513603]
        //  [   1.63840757   20.75032424    3.93096334]
        //  [  -0.20513603    3.93096334    1.87107189]]
        //
        // [[   4.25340328   25.4042382    -0.49815715]
        //  [  25.4042382  1069.83341318  -28.75048534]
        //  [  -0.49815715  -28.75048534    3.0197733 ]]]
        //degrees_of_freedom_:
        // [   3.         4492.96539243 1197.83693171 4318.19767586]
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
}