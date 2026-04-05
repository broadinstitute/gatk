package org.broadinstitute.hellbender.tools.sv.aggregation;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;
import java.util.Random;

/**
 * Simple 1D Gaussian Mixture Model fitted via Expectation-Maximization.
 *
 * <p>This implements a simplified Bayesian-style GMM for 1D data, matching the behavior of
 * sklearn's {@code BayesianGaussianMixture(n_components=K, covariance_type='spherical')}
 * as used in the v1.1 BAF test. The Bayesian aspect is approximated by adding a small
 * regularization to the covariance to prevent degenerate components.</p>
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

    // Fitted parameters
    private double[] weights;
    private double[] means;
    private double[] variances;

    /**
     * @param nComponents number of Gaussian components
     * @param nInit       number of random restarts (best log-likelihood wins)
     * @param maxIter     maximum EM iterations per restart
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
     * Fits the GMM on the given 1D data using EM with multiple random restarts.
     *
     * @param data array of 1D observations
     */
    public void fit(final double[] data) {
        Utils.nonNull(data);
        Utils.validateArg(data.length > 0, "data must be non-empty");

        final Random rng = new Random(seed);
        double bestLogLikelihood = Double.NEGATIVE_INFINITY;

        for (int init = 0; init < nInit; init++) {
            final double[] w = new double[nComponents];
            final double[] m = new double[nComponents];
            final double[] v = new double[nComponents];

            // Initialize means from random data points
            initializeParameters(data, w, m, v, rng);

            double prevLogLik = Double.NEGATIVE_INFINITY;
            for (int iter = 0; iter < maxIter; iter++) {
                // E-step: compute responsibilities
                final double[][] resp = new double[data.length][nComponents];
                final double logLik = eStep(data, w, m, v, resp);

                // Check convergence
                if (Math.abs(logLik - prevLogLik) < CONVERGENCE_TOL) {
                    break;
                }
                prevLogLik = logLik;

                // M-step: update parameters
                mStep(data, resp, w, m, v);
            }

            final double finalLogLik = computeLogLikelihood(data, w, m, v);
            if (finalLogLik > bestLogLikelihood) {
                bestLogLikelihood = finalLogLik;
                this.weights = w.clone();
                this.means = m.clone();
                this.variances = v.clone();
            }
        }
    }

    /**
     * Computes the average log-likelihood of test data under the fitted model.
     * This matches sklearn's {@code gmm.score(X)} which returns per-sample average log-likelihood.
     *
     * @param data test observations
     * @return average log-likelihood per observation
     */
    public double score(final double[] data) {
        Utils.nonNull(weights, "Model must be fitted before scoring");
        Utils.nonNull(data);
        Utils.validateArg(data.length > 0, "data must be non-empty");
        return computeLogLikelihood(data, weights, means, variances) / data.length;
    }

    private void initializeParameters(final double[] data, final double[] w, final double[] m,
                                      final double[] v, final Random rng) {
        final int n = data.length;
        // Uniform weights
        Arrays.fill(w, 1.0 / nComponents);

        // Pick random data points as initial means (with replacement)
        for (int k = 0; k < nComponents; k++) {
            m[k] = data[rng.nextInt(n)];
        }

        // Initialize variances to data variance
        final double mean = Arrays.stream(data).average().orElse(0);
        final double dataVar = Arrays.stream(data).map(x -> (x - mean) * (x - mean)).average().orElse(1.0);
        Arrays.fill(v, Math.max(dataVar, REG_COVAR));
    }

    /**
     * E-step: compute responsibilities and return total log-likelihood.
     */
    private double eStep(final double[] data, final double[] w, final double[] m,
                         final double[] v, final double[][] resp) {
        double totalLogLik = 0;
        for (int i = 0; i < data.length; i++) {
            double maxLogProb = Double.NEGATIVE_INFINITY;
            final double[] logProbs = new double[nComponents];
            for (int k = 0; k < nComponents; k++) {
                logProbs[k] = Math.log(w[k]) + logGaussian(data[i], m[k], v[k]);
                maxLogProb = Math.max(maxLogProb, logProbs[k]);
            }
            // Log-sum-exp for numerical stability
            double sumExp = 0;
            for (int k = 0; k < nComponents; k++) {
                sumExp += Math.exp(logProbs[k] - maxLogProb);
            }
            final double logNorm = maxLogProb + Math.log(sumExp);
            totalLogLik += logNorm;
            for (int k = 0; k < nComponents; k++) {
                resp[i][k] = Math.exp(logProbs[k] - logNorm);
            }
        }
        return totalLogLik;
    }

    /**
     * M-step: update weights, means, variances from responsibilities.
     */
    private void mStep(final double[] data, final double[][] resp,
                       final double[] w, final double[] m, final double[] v) {
        final int n = data.length;
        for (int k = 0; k < nComponents; k++) {
            double nk = 0;
            double sumX = 0;
            for (int i = 0; i < n; i++) {
                nk += resp[i][k];
                sumX += resp[i][k] * data[i];
            }
            // Avoid division by zero
            nk = Math.max(nk, 1e-10);
            w[k] = nk / n;
            m[k] = sumX / nk;
            double sumVar = 0;
            for (int i = 0; i < n; i++) {
                final double diff = data[i] - m[k];
                sumVar += resp[i][k] * diff * diff;
            }
            v[k] = Math.max(sumVar / nk, REG_COVAR);
        }
    }

    private double computeLogLikelihood(final double[] data, final double[] w,
                                        final double[] m, final double[] v) {
        double totalLogLik = 0;
        for (final double x : data) {
            double maxLogProb = Double.NEGATIVE_INFINITY;
            final double[] logProbs = new double[nComponents];
            for (int k = 0; k < nComponents; k++) {
                logProbs[k] = Math.log(w[k]) + logGaussian(x, m[k], v[k]);
                maxLogProb = Math.max(maxLogProb, logProbs[k]);
            }
            double sumExp = 0;
            for (int k = 0; k < nComponents; k++) {
                sumExp += Math.exp(logProbs[k] - maxLogProb);
            }
            totalLogLik += maxLogProb + Math.log(sumExp);
        }
        return totalLogLik;
    }

    /**
     * Log PDF of 1D Gaussian.
     */
    private static double logGaussian(final double x, final double mean, final double variance) {
        final double diff = x - mean;
        return -0.5 * (Math.log(2 * Math.PI * variance) + diff * diff / variance);
    }
}
