package org.broadinstitute.hellbender.utils.clustering;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.NonPositiveDefiniteMatrixException;
import org.apache.commons.math3.linear.NonSquareMatrixException;
import org.apache.commons.math3.linear.NonSymmetricMatrixException;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.NaturalLogUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public final class BayesianGaussianMixture {

    private static final double LOG_2_PI = Math.log(2. * Math.PI);
    private static final double EPSILON = 1E-10;
    private static final double RELATIVE_SYMMETRY_THRESHOLD = 1E-6;
    private static final double ABSOLUTE_POSITIVITY_THRESHOLD = 1E-10;

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
    private RealVector meanPrior;             // not final; if not specified via Builder and set to null, will then be set to the mean of X (over samples, i.e., over rows) upon fit(X)
    private Double degreesOfFreedomPrior;     // not final; if not specified via Builder and set to null, will then be set to nFeatures upon fit(X)
    private RealMatrix covariancePrior;       // not final; if not specified via Builder and set to null, will then be set to the covariance of X upon fix(X)
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

    private RealVector bestWeightConcentration;
    private RealVector bestMeanPrecision;
    private List<RealVector> bestMeans;
    private List<RealMatrix> bestPrecisionsCholesky;
    private List<RealMatrix> bestCovariances;
    private RealVector bestDegreesOfFreedom;

    /**
     * We use a {@link Builder} to enable the setting of default values, which can be
     * done trivially in python. Our implementation of parameter validation thus slightly differs from that in sklearn.
     */
    private BayesianGaussianMixture(final int nComponents,
                                    final double tol,
                                    final double regCovar,
                                    final int maxIter,
                                    final int nInit,
                                    final InitMethod initMethod,
                                    final double weightConcentrationPrior,
                                    final double meanPrecisionPrior,
                                    final RealVector meanPrior,
                                    final Double degreesOfFreedomPrior,
                                    final RealMatrix covariancePrior,
                                    final int seed,
                                    final boolean warmStart,
                                    final int verboseInterval) {
        this.nComponents = nComponents;
        this.tol = tol;
        this.regCovar = regCovar;
        this.maxIter = maxIter;
        this.nInit = nInit;
        this.verboseInterval = verboseInterval;
        this.weightConcentrationPrior = weightConcentrationPrior;
        this.meanPrecisionPrior = meanPrecisionPrior;
        this.meanPrior = meanPrior;
        this.degreesOfFreedomPrior = degreesOfFreedomPrior;
        this.covariancePrior = covariancePrior;
        this.initMethod = initMethod;
        this.seed = seed;
        this.warmStart = warmStart;

        rng = RandomGeneratorFactory.createRandomGenerator(new Random(seed));
        isConverged = false;
        lowerBound = Double.NEGATIVE_INFINITY;

        logger.info(toString());
    }

    /**
     * @param data matrix of data with dimensions (n_samples, n_features); to minimize memory requirements, a defensive copy will not be made
     */
    public void fit(final double[][] data) {
        Utils.validateArg(data.length >= 2, "Data must contain at least 2 samples.");

        final RealMatrix X = new Array2DRowRealMatrix(data, false); // we do not make a defensive copy

        final int nSamples = X.getRowDimension();
        if (nSamples < nComponents) {
            throw new IllegalArgumentException(
                    String.format("Number of data samples = %d is not greater than or equal to nComponents = %d",
                            nSamples, nComponents));
        }
        checkPriorParameters(X);

        final boolean doInit = !(warmStart && isConverged);
        final int nInit = doInit ? this.nInit : 1;

        double maxLowerBound = Double.NEGATIVE_INFINITY;
        isConverged = false;

        logHeapUsage("Starting loop...");
        for (int init = 0; init < nInit; init++) {
            logHeapUsage("Starting initialization...");
            logger.info(String.format("Initialization %d...", init));

            if (doInit) {
                initializeParameters(X, seed);
            }

            double lowerBound = doInit ? Double.NEGATIVE_INFINITY : this.lowerBound;

            logHeapUsage("Starting iteration...");
            for (int nIter = 1; nIter <= maxIter; nIter++) {
                final double prevLowerBound = lowerBound;

                logHeapUsage("Starting EStep...");
                final RealMatrix logResp = EStep(X);

                logHeapUsage("Starting MStep...");
                MStep(X, logResp);

                logHeapUsage("Starting computeLowerBound...");
                lowerBound = computeLowerBound(logResp);


                final double change = lowerBound - prevLowerBound;

                if (nIter % verboseInterval == 0) {
                    logger.info(String.format("Iteration %d, lower bound = %.5f, lower-bound change = %.5f...",
                            nIter, lowerBound, change));
                }

                if (Math.abs(change) < tol) {
                    isConverged = true;
                    logger.info(String.format("Initialization %d converged after %d iterations, final lower bound = %.5f, final lower-bound change = %.5f...",
                            init, nIter, lowerBound, change));
                    break;
                }

                if (nInit == maxIter && !isConverged) {
                    logger.info(String.format("Initialization %d did not converge after %d iterations, final lower bound = %.5f, final lower-bound change = %.5f...",
                            init, nIter, lowerBound, change));
                }

                if (lowerBound > maxLowerBound || maxLowerBound == Double.NEGATIVE_INFINITY) {
                    maxLowerBound = lowerBound;
                    logger.info(String.format("New maximum lower bound = %.5f found with initialization %d...",
                            maxLowerBound, init));
                    bestWeightConcentration = weightConcentration.copy();
                    bestMeanPrecision = meanPrecision.copy();
                    bestMeans = means.stream().map(RealVector::copy).collect(Collectors.toList());
                    bestPrecisionsCholesky = precisionsCholesky.stream().map(RealMatrix::copy).collect(Collectors.toList());
                    bestCovariances = covariances.stream().map(RealMatrix::copy).collect(Collectors.toList());
                    bestDegreesOfFreedom = degreesOfFreedom.copy();
                }
            }

            if (!isConverged) {
                logger.warn("No initializations converged. Try changing initialization parameters, increasing maxIter, " +
                        "increasing tol, or checking for degenerate data.");
            }

            // set values to the best found over all initializations
            this.lowerBound = maxLowerBound;
            weightConcentration = bestWeightConcentration.copy();
            meanPrecision = bestMeanPrecision.copy();
            means = bestMeans.stream().map(RealVector::copy).collect(Collectors.toList());
            precisionsCholesky = bestPrecisionsCholesky.stream().map(RealMatrix::copy).collect(Collectors.toList());
            covariances = bestCovariances.stream().map(RealMatrix::copy).collect(Collectors.toList());
            degreesOfFreedom = bestDegreesOfFreedom.copy();
        }

//        TODO add this for fitPredict
//        (original sklearn comment) Always do a final e-step to guarantee that the labels returned by
//        fit_predict(X) are always consistent with fit(X).predict(X)
//        for any value of max_iter and tol (and any random_state).
//        _, log_resp = self._e_step(X)
//
//        return log_resp.argmax(axis=1)
    }

    public double[] scoreSamples(final double[][] data) {
        return null;
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
                    return componentIndex;
                }
            });
            final RealVector respSumOverComponents = new ArrayRealVector(nSamples);
            IntStream.range(0, nSamples).forEach(
                    i -> respSumOverComponents.setEntry(i, Arrays.stream(resp.getRow(i)).sum()));
            IntStream.range(0, nSamples).forEach(
                    i -> resp.setRowVector(i, resp.getRowVector(i).mapDivide(respSumOverComponents.getEntry(i))));
            initialize(X, resp);
        } else {
            throw new GATKException.ShouldNeverReachHereException("An initialization method must be implemented for each InitMethod.");
        }
    }

    private void initialize(final RealMatrix X,
                            final RealMatrix resp) { // unlike the sklearn implementation, we simply reuse the MStep method;
        MStep(X, map(resp, FastMath::log));          // we thus compute an extra log(resp), which fortunately does not seem to cause numerical discrepancies in tests on typical data
    }                                                // (although one might be wary of diverging at all from sklearn at this initialization step...)

    private RealMatrix EStep(final RealMatrix X) {  // we trivially collapse out the sklearn _estimate_log_prob_resp method and remove some unused returns
        final int nSamples = X.getRowDimension();
        final RealMatrix weightedLogProb = estimateWeightedLogProb(X);
        final RealVector logProbNorm = new ArrayRealVector(nSamples);
        IntStream.range(0, nSamples).forEach(
                i -> logProbNorm.setEntry(i, NaturalLogUtils.logSumExp(weightedLogProb.getRow(i))));
        final RealMatrix logResp = weightedLogProb.copy();
        IntStream.range(0, nSamples).forEach(
                i -> logResp.setRowVector(i, logResp.getRowVector(i).mapSubtract(logProbNorm.getEntry(i))));
        return logResp;
    }

    private void MStep(final RealMatrix X,
                       final RealMatrix logResp) {
        final Triple<RealVector, List<RealVector>, List<RealMatrix>> parameters =
                estimateGaussianParameters(X, map(logResp, FastMath::exp), regCovar);
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
        // (original sklearn comment) Warning : in some versions of Bishop, there is a typo in the formula 10.63
        // `degrees_of_freedom_k = degrees_of_freedom_0 + N_k` is the correct formula
        degreesOfFreedom = nk.mapAdd(degreesOfFreedomPrior);

        covariances = new ArrayList<>(Collections.nCopies(nComponents, covariancePrior));
        for (int k = 0; k < nComponents; k++) {
            final RealVector diff = xk.get(k).subtract(meanPrior);
            final RealMatrix cov = covariancePrior
                    .add(sk.get(k).scalarMultiply(nk.getEntry(k)))
                    .add(diff.outerProduct(diff).scalarMultiply(nk.getEntry(k) * meanPrecisionPrior / meanPrecision.getEntry(k)))
                    .scalarMultiply(1. / degreesOfFreedom.getEntry(k)); // (original sklearn comment) Contrary to the original Bishop book, we normalize the covariances
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

    private RealVector estimateLogWeights() {
        return weightConcentration.map(Gamma::digamma).mapSubtract(Gamma.digamma(sum(weightConcentration)));
    }

    private RealMatrix estimateLogProb(final RealMatrix X) {
        final int nSamples = X.getRowDimension();
        final int nFeatures = X.getColumnDimension();

        // (original sklearn comment) We remove nFeatures * degreesOfFreedom.map(FastMath::log) because
        // the precision matrix is normalized
        final RealMatrix result = estimateLogGaussianProb(X, means, precisionsCholesky);
        final RealVector normalizingTerm = degreesOfFreedom.map(FastMath::log).mapMultiply(0.5 * nFeatures);
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
        // (original sklearn comment) Contrary to the original formula, we have done some simplification
        // and removed all the constant terms.
        final int nFeatures = meanPrior.getDimension();

        // (original sklearn comment) We removed 0.5 * nFeatures * degreesOfFreedom.map(FastMath::log)
        // because the precision matrix is normalized.
        final RealVector logDetPrecisionsChol = computeLogDetCholesky(precisionsCholesky);
        logDetPrecisionsChol.combineToSelf(1., -0.5 * nFeatures, degreesOfFreedom.map(FastMath::log));

        final double logWishart = sum(logWishartNorm(degreesOfFreedom, logDetPrecisionsChol, nFeatures));
        final double logNormWeight = logDirichletNorm(weightConcentration);

        return -sum(ebeMultiply(map(logResp, FastMath::exp), logResp))
                - logWishart
                - logNormWeight
                - 0.5 * nFeatures * sum(meanPrecision.map(FastMath::log));
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
                .add(degreesOfFreedom.mapMultiply(0.5 * MathUtils.LOG_2 * nFeatures))
                .add(logGammaSumTerm)
                .mapMultiply(-1.);
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
                final RealMatrix covChol = new CholeskyDecomposition(
                        cov, RELATIVE_SYMMETRY_THRESHOLD, ABSOLUTE_POSITIVITY_THRESHOLD).getL();
                final RealMatrix precisionChol = cov.createMatrix(nFeatures, nFeatures);
                for (int l = 0; l < nFeatures; l++) {
                    final RealVector b = new ArrayRealVector(nFeatures);
                    b.setEntry(l, 1.);
                    MatrixUtils.solveLowerTriangularSystem(covChol, b);
                    precisionChol.setColumnVector(l, b);
                }
                precisionsChol.set(k, precisionChol.transpose());
            } catch (final NonSymmetricMatrixException | NonPositiveDefiniteMatrixException e) {
                throw new UserException(
                        "Numerical issues were encountered, causing at least one component to be assigned an invalid covariance " +
                                "(i.e., the covariance was not symmetric and/or positive semidefinite). " +
                                "Try to adjust the covariance prior or increase the covariance-regularization parameter.", e);
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
                k -> logDetChol.setEntry(k, sum(diag(matrixChol.get(k)).map(FastMath::log))));
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

    //******************************************************************************************************************
    // static RealVector and RealMatrix helper methods
    //******************************************************************************************************************

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

    //******************************************************************************************************************
    // parameter validation methods that require the data X
    //******************************************************************************************************************

    private void checkPriorParameters(final RealMatrix X) {
        checkMeanParameters(X);
        checkPrecisionParameters(X);
        checkCovarianceParameters(X);
    }

    private void checkMeanParameters(final RealMatrix X) {
        final int nSamples = X.getRowDimension();
        final int nFeatures = X.getColumnDimension();
        if (meanPrior == null) { // if not specified through the Builder, set to the mean of X (over samples, i.e., over rows) by default
            meanPrior = new ArrayRealVector(nComponents);
            IntStream.range(0, nFeatures).forEach(
                    i -> meanPrior.setEntry(i, sum(X.getColumnVector(i)) / nSamples));
        } else {
            Utils.validateArg(meanPrior.getDimension() == nFeatures,
                    String.format("meanPrior and the data should correspond to the same number of features, " +
                            "but meanPrior has dimension %d and the data has %d features.", meanPrior.getDimension(), nFeatures));
        }
    }

    private void checkPrecisionParameters(final RealMatrix X) {
        final int nFeatures = X.getColumnDimension();
        if (degreesOfFreedomPrior == null) { // if not specified through the Builder, set to nFeatures by default
            degreesOfFreedomPrior = (double) nFeatures;
        } else {
            Utils.validateArg(degreesOfFreedomPrior > nFeatures - 1.,
                    String.format("degreesOfFreedomPrior = %.5f must be strictly greater than the number of features minus 1 = %d.",
                            degreesOfFreedomPrior, nFeatures - 1));
        }
    }

    private void checkCovarianceParameters(final RealMatrix X) {
        final int nFeatures = X.getColumnDimension();
        if (covariancePrior == null) { // if not specified through the Builder, set to the covariance of X by default
            covariancePrior = new Covariance(X.transpose()).getCovarianceMatrix();
        } else {
            Builder.checkCovariance(covariancePrior);
            Utils.validateArg(covariancePrior.getRowDimension() == nFeatures,
                    String.format("covariancePrior and the data should correspond to the same number of features, " +
                            "but covariancePrior has dimensions %d x %d and the data has %d features.",
                            covariancePrior.getRowDimension(), covariancePrior.getColumnDimension(), nFeatures));
        }
    }

    private static void logHeapUsage(final String message) {
        final int mb = 1024 * 1024;
        final Runtime runtime = Runtime.getRuntime();
        logger.debug("Used memory [MB]: " + (runtime.totalMemory() - runtime.freeMemory()) / mb);
        logger.debug(message);
    }

    public RealVector getWeightConcentration() {
        return weightConcentration.copy();
    }

    public RealVector getMeanPrecision() {
        return meanPrecision.copy();
    }

    public List<RealVector> getMeans() {
        return means.stream().map(RealVector::copy).collect(Collectors.toList());
    }

    public List<RealMatrix> getPrecisionsCholesky() {
        return precisionsCholesky.stream().map(RealMatrix::copy).collect(Collectors.toList());
    }

    public List<RealMatrix> getCovariances() {
        return covariances.stream().map(RealMatrix::copy).collect(Collectors.toList());
    }

    public RealVector getDegreesOfFreedom() {
        return degreesOfFreedom.copy();
    }

    @Override
    public String toString() {
        return "BayesianGaussianMixture{" +
                "nComponents=" + nComponents +
                ", tol=" + tol +
                ", regCovar=" + regCovar +
                ", maxIter=" + maxIter +
                ", nInit=" + nInit +
                ", initMethod=" + initMethod +
                ", weightConcentrationPrior=" + weightConcentrationPrior +
                ", meanPrecisionPrior=" + meanPrecisionPrior +
                ", meanPrior=" + meanPrior +
                ", degreesOfFreedomPrior=" + degreesOfFreedomPrior +
                ", covariancePrior=" + covariancePrior +
                ", seed=" + seed +
                ", warmStart=" + warmStart +
                ", verboseInterval=" + verboseInterval +
                ", isConverged=" + isConverged +
                ", lowerBound=" + lowerBound +
                ", weightConcentration=" + weightConcentration +
                ", meanPrecision=" + meanPrecision +
                ", means=" + means +
                ", precisionsCholesky=" + precisionsCholesky +
                ", covariances=" + covariances +
                ", degreesOfFreedom=" + degreesOfFreedom +
                '}';
    }

    /**
     * This builder will set defaults and perform validation of basic parameters that do not require the data X to do so.
     */
    public static final class Builder {
        private int nComponents = 1;
        private double tol = 1E-3;
        private double regCovar = 1E-6;
        private int maxIter = 100;
        private int nInit = 1;
        private InitMethod initMethod = InitMethod.K_MEANS;
                                                                // some prior parameters require the data X to construct defaults and/or fully validate
        private Double weightConcentrationPrior = null;         // if null, will be set to 1. / nComponents upon build()
        private double meanPrecisionPrior = 1;
        private RealVector meanPrior = null;                    // if null, will be set to the mean of X (over samples, i.e., over rows) upon fix(X)
        private Double degreesOfFreedomPrior = null;            // if null, will be set to nFeatures (i.e., X.getColumnDimension()) upon fix(X)
        private RealMatrix covariancePrior;                     // if null, will be set to the covariance of X upon fix(X)

        private int seed = 0;                                   // defaults to None (i.e., uses the global random state) in the sklearn implementation; we disallow this to avoid reproducibility headaches
        private boolean warmStart = false;
        private int verboseInterval = 10;

        public Builder() {
        }

        public Builder nComponents(final int nComponents) {
            Utils.validateArg(nComponents >= 1, "nComponents must be >= 1.");
            this.nComponents = nComponents;
            return this;
        }

        public Builder tol(final double tol) {
            Utils.validateArg(tol > 0., "tol must be > 0.");
            this.tol = tol;
            return this;
        }

        public Builder regCovar(final double regCovar) {
            Utils.validateArg(regCovar >= 0., "regCovar must be >= 0.");
            this.regCovar = regCovar;
            return this;
        }

        public Builder maxIter(final int maxIter) {
            Utils.validateArg(maxIter >= 1, "maxIter must be >= 1.");
            this.maxIter = maxIter;
            return this;
        }

        public Builder nInit(final int nInit) {
            Utils.validateArg(nInit >= 1, "nInit must be >= 1.");
            this.nInit = nInit;
            return this;
        }

        public Builder initMethod(final InitMethod initMethod) {
            this.initMethod = initMethod;
            return this;
        }

        public Builder weightConcentrationPrior(final double weightConcentrationPrior) {
            Utils.validateArg(weightConcentrationPrior > 0., "weightConcentrationPrior must be > 0.");
            this.weightConcentrationPrior = weightConcentrationPrior;
            return this;
        }

        public Builder meanPrecisionPrior(final double meanPrecisionPrior) {
            Utils.validateArg(meanPrecisionPrior > 0., "meanPrecisionPrior must be > 0.");
            this.meanPrecisionPrior = meanPrecisionPrior;
            return this;
        }

        public Builder meanPrior(final double[] meanPrior) {
            Utils.nonNull(meanPrior, "meanPrior must be non-null. " +
                    "To set this to the empirical mean of the data, simply do not chain this method onto the Builder.");
            this.meanPrior = new ArrayRealVector(meanPrior, true);
            return this;
        }

        public Builder degreesOfFreedomPrior(final double degreesOfFreedomPrior) {
            Utils.validateArg(degreesOfFreedomPrior > 0., "degreesOfFreedomPrior must be > 0.");
            this.degreesOfFreedomPrior = degreesOfFreedomPrior;
            return this;
        }

        public Builder covariancePrior(final double[][] covariancePrior) {
            Utils.nonNull(covariancePrior, "covariancePrior must be non-null. " +
                    "To set this the empirical covariance of the data, simply do not chain this method onto the Builder.");
            this.covariancePrior = new Array2DRowRealMatrix(covariancePrior, true);
            checkCovariance(this.covariancePrior);
            return this;
        }

        public Builder seed(final int seed) {
            this.seed = seed;
            return this;
        }

        public Builder warmStart(final boolean warmStart) {
            this.warmStart = warmStart;
            return this;
        }

        public Builder verboseInterval(final int verboseInterval) {
            Utils.validateArg(verboseInterval >= 1, "verboseInterval must be >= 1.");
            this.verboseInterval = verboseInterval;
            return this;
        }

        public BayesianGaussianMixture build() {
            weightConcentrationPrior = weightConcentrationPrior == null ? 1. / nComponents : weightConcentrationPrior;
            return new BayesianGaussianMixture(
                    nComponents,
                    tol,
                    regCovar,
                    maxIter,
                    nInit,
                    initMethod,
                    weightConcentrationPrior,
                    meanPrecisionPrior,
                    meanPrior,
                    degreesOfFreedomPrior,
                    covariancePrior,
                    seed,
                    warmStart,
                    verboseInterval);
        }

        private static void checkCovariance(final RealMatrix covariance) {
            try {
                new CholeskyDecomposition(covariance, BayesianGaussianMixture.RELATIVE_SYMMETRY_THRESHOLD, BayesianGaussianMixture.ABSOLUTE_POSITIVITY_THRESHOLD);
            } catch (final NonSquareMatrixException | NonSymmetricMatrixException | NonPositiveDefiniteMatrixException e) {
                throw new IllegalArgumentException("Covariance matrix must be square, symmetric, and positive semidefinite.", e);
            }
        }
    }
}