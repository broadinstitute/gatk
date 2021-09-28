package org.broadinstitute.hellbender.utils.clustering;

import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class BayesianGaussianMixtureUnitTest {
    @Test
    public void testSmallInputs() {
        final int nComponents = 3;
        final double tol = 1E-3;
        final double regCovar = 1E-6;
        final int maxIter = 10;
        final int nInit = 1;
        final BayesianGaussianMixture.InitMethod initMethod = BayesianGaussianMixture.InitMethod.TEST;
        final double weightConcentrationPrior = 0.01;
        final double meanPrecisionPrior = 10.;
        final double[] meanPrior = new double[]{0., 0.};
        final double degreesOfFreedomPrior = 1.;
        final double[][] covariancePrior = new double[][]{{1., 0.}, {0., 1.}};
        final int seed = 1;
        final boolean warmStart = true;
        final int verboseInterval = 1;

        final BayesianGaussianMixture bgmm = new BayesianGaussianMixture(
                nComponents, tol, regCovar, maxIter, nInit, initMethod, weightConcentrationPrior, meanPrecisionPrior,
                meanPrior, degreesOfFreedomPrior, covariancePrior, seed, warmStart, verboseInterval);

        final double[][] data = new double[][]{
                {1., 2.},
                {2., 3.},
                {3., 4.},
                {4., 5.},
                {5., 6.}
        };

        bgmm.fit(data);
    }
}