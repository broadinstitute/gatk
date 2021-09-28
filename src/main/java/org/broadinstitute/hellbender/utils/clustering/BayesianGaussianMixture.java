package org.broadinstitute.hellbender.utils.clustering;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.NonPositiveDefiniteMatrixException;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.commons.math3.special.Gamma;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.NaturalLogUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;

public final class BayesianGaussianMixture {

    private static final double LOG_2_PI = Math.log(2. * Math.PI);
    private static final double EPSILON = 1E-10;

    public enum InitMethod {
        K_MEANS, RANDOM, TEST
    }

    protected static final Logger logger = LogManager.getLogger(BayesianGaussianMixture.class);

    private final int nComponents;
    private final double tol;
    private final double regCovar;
    private final int maxIter;
    private final int nInit;
    private final InitMethod initMethod;
    private final double weightConcentrationPrior;
    private final double meanPrecisionPrior;
    private final RealVector meanPrior;
    private final double degreesOfFreedomPrior;
    private final RealMatrix covariancePrior;
    private final int seed;
    private final boolean warmStart;
    private final int verboseInterval;

    private final RandomGenerator rng;
    private boolean isConverged;
    private double lowerBound;

    private RealVector weightConcentration;
    private RealVector meanPrecision;
    private List<RealVector> means;
    private List<RealMatrix> precisionsCholesky;
    private List<RealMatrix> covariances;
    private RealVector degreesOfFreedom;

    private void logHeapUsage(final String message) {
        final int mb = 1024 * 1024;
        final Runtime runtime = Runtime.getRuntime();
        logger.info("Used memory [MB]: " + (runtime.totalMemory() - runtime.freeMemory()) / mb);
        System.gc();
        logger.info(message);
    }

    public BayesianGaussianMixture(final int nComponents,
                                   final double tol,
                                   final double regCovar,
                                   final int maxIter,
                                   final int nInit,
                                   final InitMethod initMethod,
                                   final double weightConcentrationPrior,
                                   final double meanPrecisionPrior,
                                   final double[] meanPrior,
                                   final double degreesOfFreedomPrior,
                                   final double[][] covariancePrior,
                                   final int seed,
                                   final boolean warmStart,
                                   final int verboseInterval) {
        // TODO: use a Builder to enable default parameter values
        this.nComponents = nComponents;
        this.tol = tol;
        this.regCovar = regCovar;
        this.maxIter = maxIter;
        this.nInit = nInit;
        this.initMethod = initMethod;
        this.weightConcentrationPrior = weightConcentrationPrior;
        this.meanPrecisionPrior = meanPrecisionPrior;
        this.meanPrior = new ArrayRealVector(meanPrior, true);
        this.degreesOfFreedomPrior = degreesOfFreedomPrior;
        this.covariancePrior = new Array2DRowRealMatrix(covariancePrior, true);
        this.seed = seed;
        this.warmStart = warmStart;
        this.verboseInterval = verboseInterval;

        rng = RandomGeneratorFactory.createRandomGenerator(new Random(seed));
        isConverged = false;
        lowerBound = Double.NEGATIVE_INFINITY;
    }

    /**
     * @param data matrix of data with dimensions (n_samples, n_features); defensive copy will not be made
     */
    public void fit(final double[][] data) {
//        TODO X = self._validate_data(X, dtype=[np.float64, np.float32], ensure_min_samples=2)

        final RealMatrix X = new Array2DRowRealMatrix(data, false); // we do not make a defensive copy

        final int nSamples = X.getRowDimension();
        if (nSamples < nComponents) {
            throw new IllegalArgumentException(
                    String.format("Number of data samples = %d is not greater than or equal to nComponents = %d",
                            nSamples, nComponents));
        }
        checkInitialParameters(X);

        // if we enable warm_start, we will have a unique initialisation
        final boolean doInit = !(warmStart && isConverged);
        final int nInit = doInit ? this.nInit : 1;

        double maxLowerBound = Double.NEGATIVE_INFINITY;
        isConverged = false;

        for (int init = 0; init < nInit; init++) {
            logger.info(String.format("Initialization %d...", init));

            if (doInit) {
                initializeParameters(X, seed);
            }

            double lowerBound = doInit ? Double.NEGATIVE_INFINITY : this.lowerBound;

            for (int nIter = 1; nIter <= maxIter; nIter++) {
                final double prevLowerBound = lowerBound;

                final Pair<Double, RealMatrix> logRespAndlogProbNorm = EStep(X);
                final RealMatrix logResp = logRespAndlogProbNorm.getRight();
                MStep(X, logResp);
                lowerBound = computeLowerBound(logResp);

                final double change = lowerBound - prevLowerBound;

                if (nIter % verboseInterval == 0) {
                    logger.info(String.format("Iteration %d, lower-bound change = %.4f...", nIter, change));
                }

                if (Math.abs(change) < tol) {
                    isConverged = true;
                    logger.info(String.format("Initialization %d converged after %d iterations, final lower-bound change = %.4f...", init, nIter, change));
                    break;
                }

                if (nInit == maxIter && !isConverged) {
                    logger.info(String.format("Initialization %d did not converge after %d iterations, final lower-bound change = %.4f...", init, nIter, change));
                }
            }

            if (lowerBound > maxLowerBound || maxLowerBound == Double.NEGATIVE_INFINITY) {
                maxLowerBound = lowerBound;
//                TODO best_params = self._get_parameters()
//                best_n_iter = n_iter
            }

            if (!isConverged) {
                logger.warn("No initializations converged. Try changing initialization parameters, increasing maxIter," +
                        "increasing tol, or checking for degenerate data.");
            }

//            TODO self._set_parameters(best_params)
//            self.n_iter_ = best_n_iter
            this.lowerBound = maxLowerBound;
        }

//        TODO add this for fitPredict
//        Always do a final e-step to guarantee that the labels returned by
//        fit_predict(X) are always consistent with fit(X).predict(X)
//        for any value of max_iter and tol (and any random_state).
//        _, log_resp = self._e_step(X)
//
//        return log_resp.argmax(axis=1)
    }

    public double[] scoreSamples(final double[][] data) {
        return null;
    }

    private void checkInitialParameters(final RealMatrix X) {
    }

    private void initializeParameters(final RealMatrix X,
                                      final int seed) {
        final int nSamples = X.getRowDimension();
        rng.setSeed(seed);

        if (initMethod == InitMethod.K_MEANS) {
            // TODO add implementation of K-means
        } else if (initMethod == InitMethod.RANDOM) {
            final RealMatrix resp = X.createMatrix(nSamples, nComponents);
            resp.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
                @Override
                public double visit(int sampleIndex, int componentIndex, double value) {
                    return rng.nextDouble();
                }
            });
            final RealVector respSumOverComponents = new ArrayRealVector(nSamples);
            IntStream.range(0, nSamples).forEach(
                    i -> respSumOverComponents.setEntry(i, Arrays.stream(resp.getRow(i)).sum()));
            IntStream.range(0, nSamples).forEach(
                    i -> resp.setRowVector(i, resp.getRowVector(i).mapDivide(respSumOverComponents.getEntry(i))));
            initialize(X, resp);
        } else if (initMethod == InitMethod.TEST) {
            final RealMatrix resp = X.createMatrix(nSamples, nComponents);
            resp.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
                @Override
                public double visit(int sampleIndex, int componentIndex, double value) {
                    return 1. / nComponents;
                }
            });
            initialize(X, resp);
        } else {
            throw new GATKException.ShouldNeverReachHereException("An initialization method must be implemented for each InitMethod.");
        }
    }

    private void initialize(final RealMatrix X,
                            final RealMatrix resp) {
        final Triple<RealVector, List<RealVector>, List<RealMatrix>> parameters = estimateGaussianParameters(X, resp, regCovar);
        final RealVector nk = parameters.getLeft();
        final List<RealVector> xk = parameters.getMiddle();
        final List<RealMatrix> sk = parameters.getRight();

        estimateWeights(nk);
        estimateMeans(nk, xk);
        estimatePrecisions(nk, xk, sk);
    }

    private Pair<Double, RealMatrix> EStep(final RealMatrix X) {
        final Pair<RealVector, RealMatrix> logProbNormAndLogProbResp = estimateLogProbResp(X);
        final RealVector logProbNorm = logProbNormAndLogProbResp.getLeft();
        final RealMatrix logProbResp = logProbNormAndLogProbResp.getRight();

        final double meanLogProbNorm = Arrays.stream(logProbNorm.toArray()).average().getAsDouble();
        return Pair.of(meanLogProbNorm, logProbResp);
    }

    private void MStep(final RealMatrix X,
                       final RealMatrix logResp) {
        final Triple<RealVector, List<RealVector>, List<RealMatrix>> parameters =
                estimateGaussianParameters(X, map(logResp, Math::exp), regCovar);
        final RealVector nk = parameters.getLeft();
        final List<RealVector> xk = parameters.getMiddle();
        final List<RealMatrix> sk = parameters.getRight();

        estimateWeights(nk);
        estimateMeans(nk, xk);
        estimatePrecisions(nk, xk, sk);
    }

    private void estimateWeights(final RealVector nk) {
        weightConcentration = nk.mapAdd(weightConcentrationPrior);
    }

    private void estimateMeans(final RealVector nk,
                               final List<RealVector> xk) {
        meanPrecision = nk.mapAdd(meanPrecisionPrior);
        final List<RealVector> means = new ArrayList<>(Collections.nCopies(nComponents, meanPrior));
        IntStream.range(0, nComponents).forEach(
                k -> means.set(k,
                        meanPrior.mapMultiply(meanPrecisionPrior)
                                .add(xk.get(k).mapMultiply(nk.getEntry(k)))
                                .mapDivide(meanPrecision.getEntry(k))));
        this.means = means;
    }

    private void estimatePrecisions(final RealVector nk,
                                    final List<RealVector> xk,
                                    final List<RealMatrix> sk) {
        estimateWishartFull(nk, xk, sk);
        precisionsCholesky = computePrecisionCholesky(covariances);
    }

    private void estimateWishartFull(final RealVector nk,
                                     final List<RealVector> xk,
                                     final List<RealMatrix> sk) {
        // Warning : in some versions of Bishop, there is a typo in the formula 10.63
        // `degrees_of_freedom_k = degrees_of_freedom_0 + N_k` is the correct formula
        degreesOfFreedom = nk.mapAdd(degreesOfFreedomPrior);

        covariances = new ArrayList<>(Collections.nCopies(nComponents, covariancePrior));
        for (int k = 0; k < nComponents; k++) {
            final RealVector diff = xk.get(k).subtract(meanPrior);
            final RealMatrix cov = covariancePrior
                    .add(sk.get(k).scalarMultiply(nk.getEntry(k)))
                    .add(diff.outerProduct(diff).scalarMultiply(nk.getEntry(k) * meanPrecisionPrior / meanPrecision.getEntry(k)))
                    .scalarMultiply(1. / degreesOfFreedom.getEntry(k)); // Contrary to the original Bishop book, we normalize the covariances
            covariances.set(k, cov);
        }
    }

    private RealMatrix estimateWeightedLogProb(final RealMatrix X) {
        final int nSamples = X.getRowDimension();
        final RealMatrix weightedLogProb = estimateLogProb(X);
        final RealVector logWeights = estimateLogWeights();
        IntStream.range(0, nSamples).forEach(
                i -> weightedLogProb.setRowVector(i, weightedLogProb.getRowVector(i).add(logWeights)));
        return weightedLogProb;
    }

    private Pair<RealVector, RealMatrix> estimateLogProbResp(final RealMatrix X) {
        final int nSamples = X.getRowDimension();
        final RealMatrix weightedLogProb = estimateWeightedLogProb(X);
        final RealVector logProbNorm = new ArrayRealVector(nSamples);
        IntStream.range(0, nSamples).forEach(
                i -> logProbNorm.setEntry(i, NaturalLogUtils.logSumExp(weightedLogProb.getRow(i))));
        final RealMatrix logResp = weightedLogProb.copy();
        IntStream.range(0, nSamples).forEach(
                i -> logResp.setRowVector(i, logResp.getRowVector(i).mapSubtract(logProbNorm.getEntry(i))));
        return Pair.of(logProbNorm, logResp);
    }

    private RealVector estimateLogWeights() {
        return weightConcentration.map(Gamma::digamma)
                .mapAdd(Gamma.digamma(sum(weightConcentration)));
    }

    private RealMatrix estimateLogProb(final RealMatrix X) {
        final int nSamples = X.getRowDimension();
        final int nFeatures = X.getColumnDimension();

        // We remove nFeatures * degreesOfFreedom.map(Math::log) because
        // the precision matrix is normalized
        final RealMatrix result = estimateLogGaussianProb(X, means, precisionsCholesky);
        final RealVector normalizingTerm = degreesOfFreedom.mapMultiply(0.5 * nFeatures);
        IntStream.range(0, nSamples).forEach(
                i -> result.setRowVector(i, result.getRowVector(i).subtract(normalizingTerm)));

        final RealVector digammaSumTerm = new ArrayRealVector(nComponents);
        IntStream.range(0, nComponents).forEach(
            k -> digammaSumTerm.setEntry(k,
                    IntStream.range(0, nFeatures)
                            .mapToDouble(i -> Gamma.digamma(0.5 * (degreesOfFreedom.getEntry(k) - i)))
                            .sum()));
        final RealVector logLambda = digammaSumTerm.mapAdd(nFeatures * MathUtils.LOG_2);

        final RealVector addToLogGaussTerm = logLambda.subtract(meanPrecision.map(x -> nFeatures / x)).mapMultiply(0.5);
        IntStream.range(0, nSamples).forEach(
                i -> result.setRowVector(i, result.getRowVector(i).add(addToLogGaussTerm)));

        return result;
    }

    private double computeLowerBound(final RealMatrix logResp) { // we remove an unused logProbNorm parameter from the python code
        // Contrary to the original formula, we have done some simplification
        // and removed all the constant terms.
        final int nFeatures = meanPrior.getDimension();

        // We removed 0.5 * nFeatures * degreesOfFreedom.map(Math::log)
        // because the precision matrix is normalized.
        final RealVector logDetPrecisionsChol = computeLogDetCholesky(precisionsCholesky);
        logDetPrecisionsChol.subtract(degreesOfFreedom.mapMultiply(0.5 * nFeatures));

        final double logWishart = sum(logWishartNorm(degreesOfFreedom, logDetPrecisionsChol, nFeatures));
        final double logNormWeight = logDirichletNorm(weightConcentration);

        return -sum(ebeMultiply(map(logResp, Math::exp), logResp)) - logWishart - logNormWeight
                - 0.5 * nFeatures * sum(meanPrecision.map(Math::log));
    }

    //******************************************************************************************************************
    // static helper methods (package-protected only for testing)
    //******************************************************************************************************************

    /**
     * See <a href="https://github.com/scikit-learn/scikit-learn/blob/1.0/sklearn/mixture/_bayesian_mixture.py#L20">here</a>.
     * @param dirichletConcentration parameters of the Dirichlet distribution. (nComponents, )
     * @return                       log normalization of the Dirichlet distribution.
     */
    @VisibleForTesting
    static double logDirichletNorm(final RealVector dirichletConcentration) {
        return Gamma.logGamma(sum(dirichletConcentration)) - sum(dirichletConcentration.map(Gamma::logGamma));
    }

    /**
     * See <a href="https://github.com/scikit-learn/scikit-learn/blob/1.0/sklearn/mixture/_bayesian_mixture.py#L38">here</a>.
     * @param degreesOfFreedom      number of degrees of freedom on the covariance Wishart distributions for all components. (nComponents, )
     * @param logDetPrecisionChol   determinants of the precision matrices for all components. (nComponents, )
     * @return                      log normalizations of the Wishart distributions for all components. (nComponents, )
     */
    @VisibleForTesting
    static RealVector logWishartNorm(final RealVector degreesOfFreedom,
                                     final RealVector logDetPrecisionChol,
                                     final int nFeatures) {
        final int nComponents = degreesOfFreedom.getDimension();
        final RealVector logGammaSumTerm = new ArrayRealVector(nComponents);
        IntStream.range(0, nComponents).forEach(
                k -> logGammaSumTerm.setEntry(k,
                        IntStream.range(0, nFeatures)
                                .mapToDouble(i -> Gamma.logGamma(0.5 * (degreesOfFreedom.getEntry(k) - i)))
                                .sum()));
        // To simplify the computation we have removed the Math.log(Math.PI) term
        return degreesOfFreedom.ebeMultiply(logDetPrecisionChol)
                .subtract(degreesOfFreedom.mapMultiply(0.5 * MathUtils.LOG_2 * nFeatures))
                .subtract(logGammaSumTerm);
    }

    /**
     * See <a href="https://github.com/scikit-learn/scikit-learn/blob/1.0/sklearn/mixture/_gaussian_mixture.py#L300">here</a>.
     * @param covariances           covariance matrices for all components. List of (nFeatures, nFeatures) with length nComponents
     * @return                      Cholesky decompositions of precision matrices for all components. List of (nFeatures, nFeatures) with length nComponents
     */
    @VisibleForTesting
    static List<RealMatrix> computePrecisionCholesky(final List<RealMatrix> covariances) {
        final int nComponents = covariances.size();
        final int nFeatures = covariances.get(0).getRowDimension();

        final List<RealMatrix> precisionsChol = new ArrayList<>(Collections.nCopies(nComponents, null));
        for (int k = 0; k < nComponents; k++) {
            final RealMatrix cov = covariances.get(k);
            try {
                final RealMatrix covChol = new CholeskyDecomposition(cov).getL();
                final RealMatrix precisionChol = cov.createMatrix(nFeatures, nFeatures);
                for (int l = 0; l < nFeatures; l++) {
                    final RealVector b = new ArrayRealVector(nFeatures);
                    b.setEntry(l, 1.);
                    MatrixUtils.solveLowerTriangularSystem(covChol, b);
                    precisionChol.setColumnVector(l, b);
                }
                precisionsChol.set(k, precisionChol.transpose());
            } catch (final NonPositiveDefiniteMatrixException e) {
                throw new UserException(
                        "Fitting the Bayesian Gaussian mixture model failed because some components have " +
                                "ill-defined empirical covariance (perhaps caused by singleton " +
                                "or collapsed samples). Try to decrease the number of components " +
                                "or increase the covariance-regularization parameter.");
            }
        }
        return precisionsChol;
    }

    /**
     * See <a href="https://github.com/scikit-learn/scikit-learn/blob/1.0/sklearn/mixture/_gaussian_mixture.py#L354">here</a>.
     * @param matrixChol            Cholesky decompositions of matrices for all components. List of (nFeatures, nFeatures) with length nComponents
     * @return                      log determinants of the Cholesky decompositions of matrices for all components. (nComponents, )
     */
    @VisibleForTesting
    static RealVector computeLogDetCholesky(final List<RealMatrix> matrixChol) { // we remove the nFeatures parameter from the python code (unnecessary when restricting covariance type to full)
        final int nComponents = matrixChol.size();
        final RealVector logDetChol = new ArrayRealVector(nComponents);
        IntStream.range(0, nComponents).forEach(
                k -> logDetChol.setEntry(k, sum(diag(matrixChol.get(k)).map(Math::log))));
        return logDetChol;
    }

    /**
     * See <a href="https://github.com/scikit-learn/scikit-learn/blob/1.0/sklearn/mixture/_gaussian_mixture.py#L394">here</a>.
     * @param X                     data. (nSamples, nFeatures)
     * @param means                 mean vectors for all components. List of (nFeatures, ) with length nComponents
     * @param precisionsChol        Cholesky decompositions of precision matrices for all components. List of (nFeatures, nFeatures) with length nComponents
     * @return                      log Gaussian probabilities. (nSamples, nComponents)
     */
    @VisibleForTesting
    static RealMatrix estimateLogGaussianProb(final RealMatrix X,
                                              final List<RealVector> means,
                                              final List<RealMatrix> precisionsChol) {
        final int nSamples = X.getRowDimension();
        final int nComponents = means.size();
        final int nFeatures = precisionsChol.get(0).getRowDimension();

        // det(precisionChol) is half of det(precision)
        final RealVector logDet = computeLogDetCholesky(precisionsChol);

        final RealMatrix logProb = X.createMatrix(nSamples, nComponents);
        for (int k = 0; k < nComponents; k++) {
            final RealVector mu = means.get(k);
            final RealMatrix precChol = precisionsChol.get(k);
            final RealVector muDotPrecChol = precChol.preMultiply(mu);
            final RealMatrix y = X.multiply(precChol);
            IntStream.range(0, nSamples).forEach(
                    i -> y.setRowVector(i, y.getRowVector(i).subtract(muDotPrecChol)));
            final RealVector ySquaredSum = new ArrayRealVector(nSamples);
            IntStream.range(0, nSamples).forEach(
                    i -> ySquaredSum.setEntry(i, y.getRowVector(i).dotProduct(y.getRowVector(i))));
            logProb.setColumnVector(k, ySquaredSum);
        }

        final RealMatrix result = logProb.scalarAdd(nFeatures * LOG_2_PI).scalarMultiply(-0.5);
        IntStream.range(0, nSamples).forEach(
                i -> result.setRowVector(i, result.getRowVector(i).add(logDet)));
        return result;
    }

    /**
     * See <a href="https://github.com/scikit-learn/scikit-learn/blob/1.0/sklearn/mixture/_gaussian_mixture.py#L260">here</a>.
     * @param X                     data. (nSamples, nFeatures)
     * @param resp                  responsibilities for each data sample in X. (nSamples, nComponents)
     * @param regCovar              regularization added to the diagonal of the covariance matrices.
     * @return                      triple of effective number for all components, mean vectors for all components, and covariance matrices for all components.
     *                              i.e., Triple of [(nComponents, ), List of (nFeatures, ) with length nComponents, List of (nFeatures, nFeatures) with length nComponents]
     */
    @VisibleForTesting
    static Triple<RealVector, List<RealVector>, List<RealMatrix>> estimateGaussianParameters(final RealMatrix X,
                                                                                             final RealMatrix resp,
                                                                                             final double regCovar) {
        final int nComponents = resp.getColumnDimension();

        final RealVector nk = new ArrayRealVector(nComponents);
        IntStream.range(0, nComponents).forEach(
                k -> nk.setEntry(k, sum(resp.getColumnVector(k)) + EPSILON));

        final List<RealVector> means = new ArrayList<>(Collections.nCopies(nComponents, null));
        final RealMatrix respTDotX = resp.transpose().multiply(X);
        IntStream.range(0, nComponents).forEach(
                k -> means.set(k, respTDotX.getRowVector(k).mapDivide(nk.getEntry(k))));

        final List<RealMatrix> covariances = estimateGaussianCovariances(X, resp, nk, means, regCovar);

        return Triple.of(nk, means, covariances);
    }


    /**
     * See <a href="https://github.com/scikit-learn/scikit-learn/blob/1.0/sklearn/mixture/_gaussian_mixture.py#L154">here</a>.
     * @param X                     data. (nSamples, nFeatures)
     * @param resp                  responsibilities for each data sample in X. (nSamples, nComponents)
     * @param nk                    effective number for all components. (nComponents, )
     * @param means                 mean vectors for all components. List of (nFeatures, ) with length nComponents
     * @param regCovar              regularization added to the diagonal of the covariance matrices.
     * @return                      covariance matrices for all components. List of (nFeatures, nFeatures) with length nComponents
     */
    private static List<RealMatrix> estimateGaussianCovariances(final RealMatrix X,
                                                                final RealMatrix resp,
                                                                final RealVector nk,
                                                                final List<RealVector> means,
                                                                final double regCovar) {
        final int nSamples = X.getRowDimension();
        final int nComponents = means.size();
        final int nFeatures = means.get(0).getDimension();

        final List<RealMatrix> covariances = new ArrayList<>(Collections.nCopies(nComponents, null));
        for (int k = 0; k < nComponents; k++) {
            final RealMatrix diff = X.copy();
            final RealVector mean = means.get(k);
            IntStream.range(0, nSamples).forEach(
                    i -> diff.setRowVector(i, diff.getRowVector(i).subtract(mean)));

            final RealVector respComponent = resp.getColumnVector(k);
            final RealMatrix respComponentTimesDiffT = diff.transpose();
            IntStream.range(0, nFeatures).forEach(
                    i -> respComponentTimesDiffT.setRowVector(i, respComponentTimesDiffT.getRowVector(i).ebeMultiply(respComponent)));
            final RealMatrix cov = respComponentTimesDiffT.multiply(diff).scalarMultiply(1. / nk.getEntry(k));

            IntStream.range(0, nFeatures).forEach(
                    i -> cov.addToEntry(i, i, regCovar));

            covariances.set(k, cov);
        }
        return covariances;
    }

    private static double sum(final RealVector v) {
        return Arrays.stream(v.toArray()).sum();
    }

    private static double sum(final RealMatrix m) {
        double sum = 0.;
        for (int i = 0; i < m.getRowDimension(); i++) {
            for (int j = 0; j < m.getColumnDimension(); j++) {
                sum += m.getEntry(i, j);
            }
        }
        return sum;
    }

    private static RealVector diag(final RealMatrix m) { // m should be square matrix
        final int dim = m.getRowDimension();
        final RealVector diag = new ArrayRealVector(dim);
        IntStream.range(0, dim).forEach(i -> diag.setEntry(i, m.getEntry(i, i)));
        return diag;
    }

    private static RealMatrix map(final RealMatrix m,
                                  final UnivariateFunction function) {
        final RealMatrix result = m.copy();
        result.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int i, int j, double value) {
                return function.value(value);
            }
        });
        return result;
    }

    private static RealMatrix ebeMultiply(final RealMatrix m,
                                          final RealMatrix n) { // m and n should have same dimensions
        final RealMatrix result = m.copy();
        result.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int i, int j, double value) {
                return m.getEntry(i, j) * n.getEntry(i, j);
            }
        });
        return result;
    }

    private void checkParameters(final RealMatrix X) {}

    private void checkWeightsParameters() {}

    private void checkMeansParameters(final RealMatrix X) {}

    private void checkPrecisionParameters(final RealMatrix X) {}

    private void checkCovariancePriorParameter(final RealMatrix X) {}

    private void checkPrecisionPositivity(final RealVector precision) {}

    private void checkPrecisionMatrix(final RealMatrix precision) {}
}