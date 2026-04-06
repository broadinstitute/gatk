package org.broadinstitute.hellbender.tools.sv.aggregation;

import org.apache.commons.math3.special.Gamma;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;
import java.util.Random;

/**
 * Variational Bayesian 1D Gaussian Mixture Model matching sklearn's
 * {@code BayesianGaussianMixture(n_components=K, covariance_type='spherical',
 * weight_concentration_prior_type='dirichlet_process')}.
 *
 * <p>This implements variational inference with:
 * <ul>
 *   <li>Dirichlet Process prior on component weights (stick-breaking representation)</li>
 *   <li>Gaussian prior on component means</li>
 *   <li>Gamma prior on component precisions (1D specialization of Wishart)</li>
 * </ul>
 *
 * <p>The Bayesian priors regularize the model by:
 * <ul>
 *   <li>Pruning unnecessary components (weights shrink toward zero)</li>
 *   <li>Preventing variance collapse (precision regularized toward data variance)</li>
 * </ul>
 *
 * <p>Usage:
 * <pre>
 *     SimpleGMM1D gmm = new SimpleGMM1D(3, 20, 200, seed);
 *     gmm.fit(data);
 *     double score = gmm.score(testData);
 * </pre>
 */
public class SimpleGMM1D {

    private static final double REG_COVAR = 1e-6;
    private static final double CONVERGENCE_TOL = 1e-3;

    private final int nComponents;
    private final int nInit;
    private final int maxIter;
    private final long seed;

    // Hyperparameters (set during fit from data)
    private double weightConcentrationPrior; // alpha_0
    private double meanPrecisionPrior;       // beta_0
    private double meanPrior;                // m_0
    private double degreesOfFreedomPrior;    // nu_0
    private double covariancePrior;          // W_0

    // Posterior parameters (fitted)
    private double[] gamma1;            // DP stick-breaking: 1 + N_k
    private double[] gamma2;            // DP stick-breaking: alpha_0 + cumsum(N_j, j>k)
    private double[] meanPrecision;     // beta_k = beta_0 + N_k
    private double[] means;             // m_k (posterior means)
    private double[] degreesOfFreedom;  // nu_k = nu_0 + N_k
    private double[] covariances;       // posterior covariance (normalized by nu_k)

    /**
     * @param nComponents number of Gaussian components
     * @param nInit       number of random restarts (best lower bound wins)
     * @param maxIter     maximum variational EM iterations per restart
     * @param seed        random seed
     */
    public SimpleGMM1D(final int nComponents, final int nInit, final int maxIter, final long seed) {
        Utils.validateArg(nComponents > 0, "nComponents must be positive");
        Utils.validateArg(nInit > 0, "nInit must be positive");
        Utils.validateArg(maxIter > 0, "maxIter must be positive");
        this.nComponents = nComponents;
        this.nInit = nInit;
        this.maxIter = maxIter;
        this.seed = seed;
    }

    /**
     * Fits the Bayesian GMM using variational EM with multiple random restarts.
     *
     * @param data array of 1D observations
     */
    public void fit(final double[] data) {
        Utils.nonNull(data);
        Utils.validateArg(data.length > 0, "data must be non-empty");

        final int n = data.length;
        final Random rng = new Random(seed);

        // Set hyperparameters (matching sklearn defaults)
        weightConcentrationPrior = 1.0 / nComponents;
        meanPrecisionPrior = 1.0;
        meanPrior = Arrays.stream(data).average().orElse(0);
        degreesOfFreedomPrior = 1.0; // n_features = 1 for 1D
        // Sample variance with ddof=1 (matching sklearn for spherical covariance)
        final double dataMean = meanPrior;
        covariancePrior = Arrays.stream(data)
                .map(x -> (x - dataMean) * (x - dataMean))
                .sum() / (n - 1);
        covariancePrior = Math.max(covariancePrior, REG_COVAR);

        double bestLowerBound = Double.NEGATIVE_INFINITY;

        for (int init = 0; init < nInit; init++) {
            // Initialize with random responsibilities
            final double[][] resp = new double[n][nComponents];
            for (int i = 0; i < n; i++) {
                double sum = 0;
                for (int k = 0; k < nComponents; k++) {
                    resp[i][k] = rng.nextDouble();
                    sum += resp[i][k];
                }
                for (int k = 0; k < nComponents; k++) {
                    resp[i][k] /= sum;
                }
            }

            // Posterior parameter arrays for this restart
            final double[] g1 = new double[nComponents];
            final double[] g2 = new double[nComponents];
            final double[] mp = new double[nComponents];
            final double[] m = new double[nComponents];
            final double[] dof = new double[nComponents];
            final double[] cov = new double[nComponents];

            // Initial M-step from random responsibilities
            mStep(data, resp, g1, g2, mp, m, dof, cov);

            // Variational EM loop
            double prevLB = Double.NEGATIVE_INFINITY;
            for (int iter = 0; iter < maxIter; iter++) {
                // E-step: compute responsibilities and lower bound proxy
                final double lb = eStep(data, g1, g2, mp, m, dof, cov, resp);

                // Check convergence (average gain per sample, matching sklearn)
                if (Math.abs((lb - prevLB) / n) < CONVERGENCE_TOL) {
                    break;
                }
                prevLB = lb;

                // M-step: update posterior parameters
                mStep(data, resp, g1, g2, mp, m, dof, cov);
            }

            // Use final lower bound for model selection
            final double finalLB = eStep(data, g1, g2, mp, m, dof, cov, resp);
            if (finalLB > bestLowerBound) {
                bestLowerBound = finalLB;
                this.gamma1 = g1.clone();
                this.gamma2 = g2.clone();
                this.meanPrecision = mp.clone();
                this.means = m.clone();
                this.degreesOfFreedom = dof.clone();
                this.covariances = cov.clone();
            }
        }
    }

    /**
     * Computes the average log-likelihood of test data under the fitted Bayesian model.
     * Matches sklearn's {@code BayesianGaussianMixture.score(X)}.
     *
     * @param data test observations
     * @return average log-likelihood per observation
     */
    public double score(final double[] data) {
        Utils.nonNull(gamma1, "Model must be fitted before scoring");
        Utils.nonNull(data);
        Utils.validateArg(data.length > 0, "data must be non-empty");

        final double[] logWeights = estimateLogWeights(gamma1, gamma2);
        double totalLogLik = 0;
        for (final double x : data) {
            final double[] logProbs = new double[nComponents];
            double maxLogProb = Double.NEGATIVE_INFINITY;
            for (int k = 0; k < nComponents; k++) {
                logProbs[k] = logWeights[k] + estimateLogProb(x, means[k], covariances[k],
                        degreesOfFreedom[k], meanPrecision[k]);
                maxLogProb = Math.max(maxLogProb, logProbs[k]);
            }
            double sumExp = 0;
            for (int k = 0; k < nComponents; k++) {
                sumExp += Math.exp(logProbs[k] - maxLogProb);
            }
            totalLogLik += maxLogProb + Math.log(sumExp);
        }
        return totalLogLik / data.length;
    }

    /**
     * M-step: update posterior parameters from current responsibilities.
     *
     * <p>Follows sklearn's {@code BayesianGaussianMixture._m_step}
     * for spherical covariance with Dirichlet Process weight prior.</p>
     */
    private void mStep(final double[] data, final double[][] resp,
                       final double[] g1, final double[] g2,
                       final double[] mp, final double[] m,
                       final double[] dof, final double[] cov) {
        final int n = data.length;

        // Compute sufficient statistics: nk, xk (weighted mean), sk (weighted variance)
        final double[] nk = new double[nComponents];
        final double[] xk = new double[nComponents];
        final double[] sk = new double[nComponents];

        for (int k = 0; k < nComponents; k++) {
            for (int i = 0; i < n; i++) {
                nk[k] += resp[i][k];
            }
            nk[k] = Math.max(nk[k], 1e-10);
            for (int i = 0; i < n; i++) {
                xk[k] += resp[i][k] * data[i];
            }
            xk[k] /= nk[k];
            for (int i = 0; i < n; i++) {
                final double diff = data[i] - xk[k];
                sk[k] += resp[i][k] * diff * diff;
            }
            sk[k] /= nk[k];
            sk[k] += REG_COVAR;
        }

        // Update weight concentration (Dirichlet Process stick-breaking)
        for (int k = 0; k < nComponents; k++) {
            g1[k] = 1.0 + nk[k];
        }
        g2[nComponents - 1] = weightConcentrationPrior;
        double cumSumRight = 0;
        for (int k = nComponents - 1; k >= 1; k--) {
            cumSumRight += nk[k];
            g2[k - 1] = weightConcentrationPrior + cumSumRight;
        }

        // Update mean parameters (Gaussian prior)
        for (int k = 0; k < nComponents; k++) {
            mp[k] = meanPrecisionPrior + nk[k];
            m[k] = (meanPrecisionPrior * meanPrior + nk[k] * xk[k]) / mp[k];
        }

        // Update precision parameters (Wishart/Gamma prior for spherical 1D)
        for (int k = 0; k < nComponents; k++) {
            dof[k] = degreesOfFreedomPrior + nk[k];
            final double diff = xk[k] - meanPrior;
            cov[k] = (covariancePrior + nk[k]
                    * (sk[k] + meanPrecisionPrior / mp[k] * diff * diff)) / dof[k];
            cov[k] = Math.max(cov[k], REG_COVAR);
        }
    }

    /**
     * E-step: compute responsibilities from current posterior parameters.
     * Returns the sum of log-normalizers (proxy for variational lower bound).
     */
    private double eStep(final double[] data,
                         final double[] g1, final double[] g2,
                         final double[] mp, final double[] m,
                         final double[] dof, final double[] cov,
                         final double[][] resp) {
        final double[] logWeights = estimateLogWeights(g1, g2);
        double totalLogNorm = 0;

        for (int i = 0; i < data.length; i++) {
            double maxLogProb = Double.NEGATIVE_INFINITY;
            final double[] logProbs = new double[nComponents];
            for (int k = 0; k < nComponents; k++) {
                logProbs[k] = logWeights[k]
                        + estimateLogProb(data[i], m[k], cov[k], dof[k], mp[k]);
                maxLogProb = Math.max(maxLogProb, logProbs[k]);
            }
            double sumExp = 0;
            for (int k = 0; k < nComponents; k++) {
                sumExp += Math.exp(logProbs[k] - maxLogProb);
            }
            final double logNorm = maxLogProb + Math.log(sumExp);
            totalLogNorm += logNorm;
            for (int k = 0; k < nComponents; k++) {
                resp[i][k] = Math.exp(logProbs[k] - logNorm);
            }
        }
        return totalLogNorm;
    }

    /**
     * Expected log weights under Dirichlet Process stick-breaking prior.
     *
     * <p>Matches sklearn's {@code _estimate_log_weights()} for dirichlet_process.</p>
     */
    private double[] estimateLogWeights(final double[] g1, final double[] g2) {
        final double[] logWeights = new double[nComponents];
        double cumLogStick = 0;
        for (int k = 0; k < nComponents; k++) {
            final double digammaSum = Gamma.digamma(g1[k] + g2[k]);
            logWeights[k] = Gamma.digamma(g1[k]) - digammaSum + cumLogStick;
            cumLogStick += Gamma.digamma(g2[k]) - digammaSum;
        }
        return logWeights;
    }

    /**
     * Expected log probability of data point x under component k.
     *
     * <p>Matches sklearn's {@code _estimate_log_prob()} for BayesianGMM
     * with spherical covariance, specialized to 1D. The formula is:</p>
     * <pre>
     *   log_gauss(x; m_k, cov_k) - 0.5*log(nu_k) + 0.5*(log(2) + digamma(nu_k/2) - 1/beta_k)
     * </pre>
     * <p>where the first term is a standard Gaussian log-pdf using the posterior
     * mean and covariance (already normalized by degrees of freedom), and the
     * remaining terms are Bayesian correction factors from the expected log-precision
     * (Wishart) and mean uncertainty.</p>
     */
    private static double estimateLogProb(final double x, final double mean,
                                          final double cov, final double dof,
                                          final double mp) {
        // Standard Gaussian log-pdf with posterior precision = 1/cov
        final double precision = 1.0 / cov;
        final double diff = x - mean;
        final double logGaussStd = -0.5 * (Math.log(2.0 * Math.PI)
                + diff * diff * precision - Math.log(precision));

        // Bayesian corrections from sklearn's _estimate_log_prob:
        // (1) -0.5 * log(nu_k): adjusts for covariance normalization
        // (2) +0.5 * E[log|Lambda_k|]: for 1D = log(2) + digamma(nu_k / 2)
        // (3) -0.5 / beta_k: mean uncertainty term
        final double logLambda = Math.log(2.0) + Gamma.digamma(0.5 * dof);
        return logGaussStd - 0.5 * Math.log(dof) + 0.5 * (logLambda - 1.0 / mp);
    }
}
