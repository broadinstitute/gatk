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
import org.apache.commons.math3.ml.clustering.CentroidCluster;
import org.apache.commons.math3.ml.clustering.Clusterable;
import org.apache.commons.math3.ml.clustering.DoublePoint;
import org.apache.commons.math3.ml.clustering.KMeansPlusPlusClusterer;
import org.apache.commons.math3.ml.distance.EuclideanDistance;
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
        K_MEANS_PLUS_PLUS, RANDOM, TEST
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

    private boolean isConverged;
    private double lowerBound;

    private RealVector weightConcentration;
    private RealVector meanPrecision;
    private List<RealVector> means;
    private List<RealMatrix> precisionsCholesky;
    private List<RealMatrix> covariances;
    private RealVector degreesOfFreedom;

    private int bestInit;
    private RealVector bestWeightConcentration;
    private RealVector bestMeanPrecision;
    private List<RealVector> bestMeans;
    private List<RealMatrix> bestPrecisionsCholesky;
    private List<RealMatrix> bestCovariances;
    private RealVector bestDegreesOfFreedom;
    private boolean isFitted = false;

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

        isConverged = false;
        lowerBound = Double.NEGATIVE_INFINITY;

        logger.info(toString());
    }

    /**
     * TODO
     * @param data      double[][] of data with dimensions (nSamples, nFeatures); to minimize memory requirements, a defensive copy will not be made
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

        final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(seed));

        logHeapUsage("Starting loop...");
        for (int init = 0; init < nInit; init++) {
            logHeapUsage("Starting initialization...");
            logger.info(String.format("Initialization %d...", init));

            if (doInit) {
                initializeParameters(X, rng);
            }

            double lowerBound = doInit ? Double.NEGATIVE_INFINITY : this.lowerBound;

            logHeapUsage("Starting iteration...");
            for (int nIter = 1; nIter <= maxIter; nIter++) {
                final double prevLowerBound = lowerBound;

                final RealMatrix logResp = EStep(X);
                MStep(X, logResp);
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

                if (nIter == maxIter && !isConverged) {
                    logger.info(String.format("Initialization %d did not converge after %d iterations, final lower bound = %.5f, final lower-bound change = %.5f...",
                            init, nIter, lowerBound, change));
                }
            }

            if (lowerBound > maxLowerBound || maxLowerBound == Double.NEGATIVE_INFINITY) {
                maxLowerBound = lowerBound;
                logger.info(String.format("New maximum lower bound = %.5f found with initialization %d...",
                        maxLowerBound, init));
                bestInit = init;
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
        logger.info(String.format("Fit complete. Maximum lower bound = %.5f found with initialization %d.",
                maxLowerBound, bestInit));
        weightConcentration = bestWeightConcentration.copy();
        meanPrecision = bestMeanPrecision.copy();
        means = bestMeans.stream().map(RealVector::copy).collect(Collectors.toList());
        precisionsCholesky = bestPrecisionsCholesky.stream().map(RealMatrix::copy).collect(Collectors.toList());
        covariances = bestCovariances.stream().map(RealMatrix::copy).collect(Collectors.toList());
        degreesOfFreedom = bestDegreesOfFreedom.copy();
        isFitted = true;
    }

    /**
     * TODO
     * In contrast to the sklearn implementation, our fitPredict calls fit (and not vice versa) to avoid an unused E step.
     * @param data      double[][] of data with dimensions (nSamples, nFeatures); to minimize memory requirements, a defensive copy will not be made
     * @return          int[] with dimension (nSamples, ) giving component assignments for each data sample
     */
    public int[] fitPredict(final double[][] data) {
        fit(data);

        final RealMatrix X = new Array2DRowRealMatrix(data, false); // we do not make a defensive copy
        final int nSamples = X.getRowDimension();

        // (original sklearn comment) Always do a final E step to guarantee that the labels returned by
        // fit_predict(X) are always consistent with fit(X).predict(X)
        // for any value of maxIter and tol (and any random seed).
        final RealMatrix logResp = EStep(X);

        return IntStream.range(0, nSamples).map(    // logResp.argmax(axis=1)
                i -> logResp.getRowVector(i).getMaxIndex())
                .toArray();
    }

    /**
     * TODO decide which log-likelihood score to report, Gaussian mixture or exact (i.e., Bishop 10.81-82)?
     * @param data      double[][] of data with dimensions (nSamples, nFeatures); to minimize memory requirements, a defensive copy will not be made
     * @return
     */
    public double[] scoreSamples(final double[][] data) {
        if (!isFitted) {
            throw new UnsupportedOperationException("Cannot score samples before model has been fit.");
        }
        final RealMatrix X = new Array2DRowRealMatrix(data);
        final int nSamples = X.getRowDimension();
        final RealMatrix weightedLogProb = estimateWeightedLogProb(X);
        final RealVector logProbNorm = new ArrayRealVector(                     // logProbNorm = logsumexp(weightedLogProb, axis=1)
                IntStream.range(0, nSamples).mapToDouble(
                        i -> NaturalLogUtils.logSumExp(weightedLogProb.getRow(i))).toArray());
        return logProbNorm.toArray();
    }

    /**
     * @param X         data. (nSamples, nFeatures)
     * @param rng       random number generator, which should be reset to this.seed at start of fit
     */
    private void initializeParameters(final RealMatrix X,
                                      final RandomGenerator rng) {
        final int nSamples = X.getRowDimension();

        if (initMethod == InitMethod.K_MEANS_PLUS_PLUS) {
            final List<IndexedDoublePoint> pointsX = IntStream.range(0, nSamples).boxed()
                    .map(i -> new IndexedDoublePoint(X.getRow(i), i))
                    .collect(Collectors.toList());

            final KMeansPlusPlusClusterer<IndexedDoublePoint> kMeansPlusPlusClusterer =
                    new KMeansPlusPlusClusterer<>(nComponents, -1, new EuclideanDistance(), rng);
            final List<CentroidCluster<IndexedDoublePoint>> centroids = kMeansPlusPlusClusterer.cluster(pointsX);
            final RealMatrix resp = X.createMatrix(nSamples, nComponents);
            IntStream.range(0, nComponents)
                    .forEach(k -> centroids.get(k).getPoints()
                            .forEach(p -> resp.setEntry(p.index, k, 1)));
            initialize(X, resp);
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

    /**
     * This could be replaced with MStep(X, map(resp, FastMath::log)), but we retain it to match the sklearn implementation.
     * @param resp      responsibilities. (nSamples, nComponents)
     */
    private void initialize(final RealMatrix X,
                            final RealMatrix resp) {
        final Triple<RealVector, List<RealVector>, List<RealMatrix>> parameters =
                estimateGaussianParameters(X, resp, regCovar);
        final RealVector nk = parameters.getLeft();
        final List<RealVector> xk = parameters.getMiddle();
        final List<RealMatrix> sk = parameters.getRight();

        estimateWeights(nk);
        estimateMeans(nk, xk);
        estimatePrecisions(nk, xk, sk);
    }

    /**
     * We trivially collapse out the sklearn _estimate_log_prob_resp method and remove some unused variables.
     * @param X         data. (nSamples, nFeatures)
     * @return          log responsibilities. (nSamples, nComponents)
     */
    private RealMatrix EStep(final RealMatrix X) {
        final int nSamples = X.getRowDimension();
        final RealMatrix result = estimateWeightedLogProb(X);                   // result = weightedLogProb; we will update result in place
        final RealVector logProbNorm = new ArrayRealVector(                     // logProbNorm = logsumexp(weightedLogProb, axis=1)
                IntStream.range(0, nSamples).mapToDouble(
                        i -> NaturalLogUtils.logSumExp(result.getRow(i))).toArray());
        IntStream.range(0, nSamples).forEach(                                // result = result - logProbNorm[:, np.newaxis]
                i -> result.setRowVector(i, result.getRowVector(i).mapSubtract(logProbNorm.getEntry(i))));
        return result;
    }

    /**
     * @param logResp   log responsibilities. (nSamples, nComponents)
     */
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

    /**
     * @param nk        (nComponents, )
     */
    private void estimateWeights(final RealVector nk) {
        weightConcentration = nk.mapAdd(weightConcentrationPrior);
    }

    /**
     * @param nk        (nComponents, )
     * @param xk        List with length nComponents of (nFeature, ) vectors
     */
    private void estimateMeans(final RealVector nk,
                               final List<RealVector> xk) {
        meanPrecision = nk.mapAdd(meanPrecisionPrior);
        final List<RealVector> means = new ArrayList<>(Collections.nCopies(nComponents, null));
        IntStream.range(0, nComponents).forEach(            // (meanPrecisionPrior * meanPrior + nk[:, np.newaxis] * xk) / meanPrecision[:, np.newaxis]
                k -> means.set(k,
                        meanPrior.mapMultiply(meanPrecisionPrior)
                                .add(xk.get(k).mapMultiply(nk.getEntry(k)))
                                .mapDivide(meanPrecision.getEntry(k))));
        this.means = means;
    }

    /**
     * @param nk        (nComponents, )
     * @param xk        List with length nComponents of (nFeature, ) vectors
     * @param sk        List with length nComponents of (nFeature, nFeature) matrices
     */
    private void estimatePrecisions(final RealVector nk,
                                    final List<RealVector> xk,
                                    final List<RealMatrix> sk) {
        estimateWishartFull(nk, xk, sk);
        precisionsCholesky = computePrecisionCholesky(covariances);
    }

    /**
     * @param nk        (nComponents, )
     * @param xk        List with length nComponents of (nFeature, ) vectors
     * @param sk        List with length nComponents of (nFeature, nFeature) matrices
     */
    private void estimateWishartFull(final RealVector nk,
                                     final List<RealVector> xk,
                                     final List<RealMatrix> sk) {
        // (original sklearn comment) Warning : in some versions of Bishop, there is a typo in the formula 10.63
        // `degrees_of_freedom_k = degrees_of_freedom_0 + N_k` is the correct formula
        degreesOfFreedom = nk.mapAdd(degreesOfFreedomPrior);
        covariances = new ArrayList<>(Collections.nCopies(nComponents, null));
        for (int k = 0; k < nComponents; k++) {
            final RealVector diff = xk.get(k).subtract(meanPrior);
            final RealMatrix cov = covariancePrior      // (covariancePrior + nk[k] * sk[k] + nk[k] * meanPrecisionPrior / meanPrecision_[k] * np.outer(diff, diff)) / degreesOfFreedom [:, np.newaxis, np.newaxis]
                    .add(sk.get(k).scalarMultiply(nk.getEntry(k)))
                    .add(diff.outerProduct(diff).scalarMultiply(nk.getEntry(k) * meanPrecisionPrior / meanPrecision.getEntry(k)))
                    .scalarMultiply(1. / degreesOfFreedom.getEntry(k)); // (original sklearn comment) Contrary to the original Bishop book, we normalize the covariances
            covariances.set(k, cov);
        }
    }

    /**
     * @param X         data. (nSamples, nFeatures)
     * @return          (nSamples, nComponents)
     */
    private RealMatrix estimateWeightedLogProb(final RealMatrix X) {
        final int nSamples = X.getRowDimension();
        final RealMatrix result = estimateLogProb(X);                  // result = logProb; we will update result in place
        final RealVector logWeights = estimateLogWeights();
        IntStream.range(0, nSamples).forEach(                       // result = result + logWeights
                i -> result.setRowVector(i, result.getRowVector(i).add(logWeights)));
        return result;
    }

    private RealVector estimateLogWeights() {
        return weightConcentration.map(Gamma::digamma).mapSubtract(Gamma.digamma(sum(weightConcentration)));
    }

    /**
     * @param X         data. (nSamples, nFeatures)
     * @return          (nSamples, nComponents)
     */
    private RealMatrix estimateLogProb(final RealMatrix X) {
        final int nSamples = X.getRowDimension();
        final int nFeatures = X.getColumnDimension();

        // (original sklearn comment) We remove nFeatures * degreesOfFreedom.map(FastMath::log) because
        // the precision matrix is normalized
        final RealMatrix result = estimateLogGaussianProb(X, means, precisionsCholesky);                           // result = estimateLogGaussianProb; we will update result in place
        final RealVector normalizingTerm = degreesOfFreedom.map(FastMath::log).mapMultiply(0.5 * nFeatures);
        IntStream.range(0, nSamples).forEach(                                                                   // result = result - 0.5 * nFeatures * np.log(degreesOfFreedom)
                i -> result.setRowVector(i, result.getRowVector(i).subtract(normalizingTerm)));

        final RealVector digammaSumTerm = new ArrayRealVector(                                                     // digammaSumTerm = np.sum(digamma(0.5 * (degreesOfFreedom - np.arange(0, nFeatures)[:, np.newaxis])), axis=0)
                IntStream.range(0, nComponents).mapToDouble(
                        k -> IntStream.range(0, nFeatures)
                                .mapToDouble(i -> Gamma.digamma(0.5 * (degreesOfFreedom.getEntry(k) - i)))
                                .sum())
                        .toArray());
        final RealVector logLambda = digammaSumTerm.mapAdd(nFeatures * MathUtils.LOG_2);                           // logLambda = digammaSumTerm + nFeatures * np.log(2.)

        final RealVector finalTerm = logLambda.subtract(meanPrecision.map(x -> nFeatures / x)).mapMultiply(0.5);   // finalTerm = 0.5 * (log_lambda - nFeatures / meanPrecision)
        IntStream.range(0, nSamples).forEach(                                                                   // result = result + finalTerm
                i -> result.setRowVector(i, result.getRowVector(i).add(finalTerm)));

        return result;
    }

    /**
     * We remove an unused log_prob_norm parameter from the sklearn implementation.
     * @param logResp       log responsibilities. (nSamples, nComponents)
     * @return              evidence lower bound
     */
    private double computeLowerBound(final RealMatrix logResp) {
        // (original sklearn comment) Contrary to the original formula, we have done some simplification
        // and removed all the constant terms.
        final int nFeatures = meanPrior.getDimension();

        // (original sklearn comment) We remove 0.5 * nFeatures * degreesOfFreedom.map(FastMath::log)
        // because the precision matrix is normalized.
        final RealVector logDetPrecisionsChol = computeLogDetCholesky(precisionsCholesky);
        logDetPrecisionsChol.combineToSelf(1., -0.5 * nFeatures, degreesOfFreedom.map(FastMath::log));      // update logDetPrecisionsChol in place

        final double logWishart = sum(logWishartNorm(degreesOfFreedom, logDetPrecisionsChol, nFeatures));
        final double logNormWeight = logDirichletNorm(weightConcentration);

        return -sum(ebeMultiply(map(logResp, FastMath::exp), logResp))
                - logWishart
                - logNormWeight
                - 0.5 * nFeatures * sum(meanPrecision.map(FastMath::log));
    }

    //******************************************************************************************************************
    // static helper methods (package-protected only for testing);
    // some of these are imported from other classes in the sklearn.mixture package, so we provide links for all methods
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
     * @param degreesOfFreedom      number of degrees of freedom of the covariance Wishart distributions for all components. (nComponents, )
     * @param logDetPrecisionsChol  determinants of the precision matrices for all components. (nComponents, )
     * @return                      log normalizations of the Wishart distributions for all components. (nComponents, )
     */
    @VisibleForTesting
    static RealVector logWishartNorm(final RealVector degreesOfFreedom,
                                     final RealVector logDetPrecisionsChol,
                                     final int nFeatures) {
        final int nComponents = degreesOfFreedom.getDimension();
        final RealVector logGammaSumTerm = new ArrayRealVector(         // logGammaSumTerm = np.sum(gammaln(0.5 * (degreesOfFreedom - np.arange(nFeatures)[:, np.newaxis])), axis=0)
                IntStream.range(0, nComponents).mapToDouble(
                        k -> IntStream.range(0, nFeatures)
                                .mapToDouble(i -> Gamma.logGamma(0.5 * (degreesOfFreedom.getEntry(k) - i)))
                                .sum())
                        .toArray());
        // (original sklearn comment) To simplify the computation we have removed the Math.log(Math.PI) term
        return degreesOfFreedom.ebeMultiply(logDetPrecisionsChol)        // -(degreesOfFreedom * logDetPrecisionsChol + degreesOfFreedom * 0.5 * np.log(2.) * nFeatures + logGammaSumTerm)
                .add(degreesOfFreedom.mapMultiply(0.5 * MathUtils.LOG_2 * nFeatures))
                .add(logGammaSumTerm)
                .mapMultiply(-1.);
    }

    /**
     * See <a href="https://github.com/scikit-learn/scikit-learn/blob/1.0/sklearn/mixture/_gaussian_mixture.py#L300">here</a>.
     * @param covariances           covariance matrices for all components. List with length nComponents of (nFeatures, nFeatures) matrices
     * @return                      Cholesky decompositions of precision matrices for all components. List with length nComponents of (nFeatures, nFeatures) matrices
     */
    @VisibleForTesting
    static List<RealMatrix> computePrecisionCholesky(final List<RealMatrix> covariances) {
        final int nComponents = covariances.size();
        final int nFeatures = covariances.get(0).getRowDimension();

        final List<RealMatrix> precisionsChol = new ArrayList<>(Collections.nCopies(nComponents, null));
        for (int k = 0; k < nComponents; k++) {
            final RealMatrix cov = covariances.get(k);
            try {
                final RealMatrix covChol = new CholeskyDecomposition(       // covChol = np.linalg.cholesky(cov, lower=True)
                        cov, RELATIVE_SYMMETRY_THRESHOLD, ABSOLUTE_POSITIVITY_THRESHOLD).getL();
                final RealMatrix precisionChol = cov.createMatrix(nFeatures, nFeatures);
                for (int l = 0; l < nFeatures; l++) {
                    final RealVector b = new ArrayRealVector(nFeatures);
                    b.setEntry(l, 1.);
                    MatrixUtils.solveLowerTriangularSystem(covChol, b);
                    precisionChol.setColumnVector(l, b);
                }
                precisionsChol.set(k, precisionChol.transpose());           // with the above loop, results in: precisionChol = np.linalg.solve_triangular(covChol, np.eye(nFeatures), lower=True).T
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
     * We remove an nFeatures parameter as it is unnecessary when restricting covariance type to full.
     * @param matrixChol            Cholesky decompositions of matrices for all components. List with length nComponents of (nFeatures, nFeatures) matrices
     * @return                      log determinants of the Cholesky decompositions of matrices for all components. (nComponents, )
     */
    @VisibleForTesting
    static RealVector computeLogDetCholesky(final List<RealMatrix> matrixChol) {
        return new ArrayRealVector(         // np.sum(np.log(matrixChol.reshape(nComponents, -1)[:, :: n_features + 1]), axis=1)
                matrixChol.stream().mapToDouble(
                        m -> sum(diag(m).map(FastMath::log)))
                        .toArray());
    }

    /**
     * See <a href="https://github.com/scikit-learn/scikit-learn/blob/1.0/sklearn/mixture/_gaussian_mixture.py#L394">here</a>.
     * @param X                     data. (nSamples, nFeatures)
     * @param means                 mean vectors for all components. List with length nComponents of (nFeatures, ) vectors
     * @param precisionsChol        Cholesky decompositions of precision matrices for all components. List with length nComponents of (nFeatures, nFeatures) matrices
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
            IntStream.range(0, nSamples).forEach(                       // y = np.dot(X, precChol) - np.dot(mu, precChol)
                    i -> y.setRowVector(i, y.getRowVector(i).subtract(muDotPrecChol)));
            final RealVector ySquaredSum = new ArrayRealVector(
                    IntStream.range(0, nSamples).mapToDouble(
                            i -> y.getRowVector(i).dotProduct(y.getRowVector(i)))
                            .toArray());
            logProb.setColumnVector(k, ySquaredSum);                       // logProb[:, k] = np.sum(np.square(y), axis=1)
        }

        final RealMatrix result = logProb.scalarAdd(nFeatures * LOG_2_PI).scalarMultiply(-0.5);
        IntStream.range(0, nSamples).forEach(                           // result = -0.5 * (nFeatures * np.log(2. * np.pi) + logProb) + logDet
                i -> result.setRowVector(i, result.getRowVector(i).add(logDet)));
        return result;
    }

    /**
     * See <a href="https://github.com/scikit-learn/scikit-learn/blob/1.0/sklearn/mixture/_gaussian_mixture.py#L260">here</a>.
     * @param X                     data. (nSamples, nFeatures)
     * @param resp                  responsibilities for each data sample in X. (nSamples, nComponents)
     * @param regCovar              regularization added to the diagonal of the covariance matrices.
     * @return                      triple of effective number for all components, mean vectors for all components, and covariance matrices for all components.
     *                              i.e., Triple of [(nComponents, ), List with length nComponents of (nFeatures, ) vectors, List with length nComponents of (nFeatures, nFeatures) matrices]
     */
    @VisibleForTesting
    static Triple<RealVector, List<RealVector>, List<RealMatrix>> estimateGaussianParameters(final RealMatrix X,
                                                                                             final RealMatrix resp,
                                                                                             final double regCovar) {
        final int nComponents = resp.getColumnDimension();

        final RealVector nk = new ArrayRealVector(
                IntStream.range(0, nComponents).mapToDouble(
                        k -> sum(resp.getColumnVector(k)) + EPSILON)
                        .toArray());

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
     * @param means                 mean vectors for all components. List with length nComponents of (nFeatures, ) vectors
     * @param regCovar              regularization added to the diagonal of the covariance matrices.
     * @return                      covariance matrices for all components. List with length nComponents of (nFeatures, nFeatures) matrices
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
            IntStream.range(0, nSamples).forEach(                                                                // diff = X - means[k]
                    i -> diff.setRowVector(i, diff.getRowVector(i).subtract(mean)));

            final RealVector respComponent = resp.getColumnVector(k);
            final RealMatrix respComponentTimesDiffT = diff.transpose();
            IntStream.range(0, nFeatures).forEach(                                                               // respComponentTimesDiffT = resp[:, k] * diff.T
                    i -> respComponentTimesDiffT.setRowVector(i, respComponentTimesDiffT.getRowVector(i).ebeMultiply(respComponent)));
            final RealMatrix cov = respComponentTimesDiffT.multiply(diff).scalarMultiply(1. / nk.getEntry(k));      // covariances[k] = np.dot(respComponentTimesDiffT, diff) / nk[k]

            IntStream.range(0, nFeatures).forEach(                                                               // covariances[k].flat[:: nFeatures + 1] += regCovar
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

    private static RealVector diag(final RealMatrix m) { // m is assumed to be a square matrix
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
                                          final RealMatrix n) { // m and n are assumed to have the same dimensions
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
    // parameter validation methods that require the data X (nSamples, nComponents)
    //******************************************************************************************************************

    /**
     * See <a href="https://github.com/scikit-learn/scikit-learn/blob/1.0/sklearn/mixture/_bayesian_mixture.py#L381">here</a>
     * for this and the following methods.
     */
    private void checkPriorParameters(final RealMatrix X) {
        checkMeanParameters(X);
        checkPrecisionParameters(X);
        checkCovarianceParameters(X);
    }

    private void checkMeanParameters(final RealMatrix X) {
        final int nSamples = X.getRowDimension();
        final int nFeatures = X.getColumnDimension();
        if (meanPrior == null) { // if not specified through the Builder, set to the mean of X (over samples, i.e., over rows) by default
            meanPrior = new ArrayRealVector(
                    IntStream.range(0, nFeatures).mapToDouble(
                            i ->  sum(X.getColumnVector(i)) / nSamples)
                            .toArray());
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
            checkCovariance(covariancePrior);
            Utils.validateArg(covariancePrior.getRowDimension() == nFeatures,
                    String.format("covariancePrior and the data should correspond to the same number of features, " +
                            "but covariancePrior has dimensions %d x %d and the data has %d features.",
                            covariancePrior.getRowDimension(), covariancePrior.getColumnDimension(), nFeatures));
        }
    }

    private static void checkCovariance(final RealMatrix covariance) {
        try {
            new CholeskyDecomposition(covariance, BayesianGaussianMixture.RELATIVE_SYMMETRY_THRESHOLD, BayesianGaussianMixture.ABSOLUTE_POSITIVITY_THRESHOLD);
        } catch (final NonSquareMatrixException | NonSymmetricMatrixException | NonPositiveDefiniteMatrixException e) {
            throw new IllegalArgumentException("Covariance matrix must be square, symmetric, and positive semidefinite.", e);
        }
    }

    //******************************************************************************************************************
    // getters, toString, and miscellaneous
    //******************************************************************************************************************

    public RealVector getWeights() {
        return weightConcentration.copy().mapDivideToSelf(sum(weightConcentration));
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
                ", seed=" + seed +
                ", warmStart=" + warmStart +
                ", verboseInterval=" + verboseInterval +
                '}';
    }

    private static void logHeapUsage(final String message) {
        final int mb = 1024 * 1024;
        final Runtime runtime = Runtime.getRuntime();
        logger.debug("Used memory [MB]: " + (runtime.totalMemory() - runtime.freeMemory()) / mb);
        logger.debug(message);
    }

    private static final class IndexedDoublePoint implements Clusterable {
        final DoublePoint point;
        final int index;

        private IndexedDoublePoint(final double[] point,
                                   final int index) {
            this.point = new DoublePoint(point);
            this.index = index;
        }

        @Override
        public double[] getPoint() {
            return point.getPoint();
        }
    }

    //******************************************************************************************************************
    // Builder
    //******************************************************************************************************************

    /**
     * This builder will set defaults and perform validation of basic parameters that do not require the data X to do so.
     */
    public static final class Builder {
        private int nComponents = 1;
        private double tol = 1E-3;
        private double regCovar = 1E-6;
        private int maxIter = 100;
        private int nInit = 1;
        private InitMethod initMethod = InitMethod.K_MEANS_PLUS_PLUS;
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
    }
}