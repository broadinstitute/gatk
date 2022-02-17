package org.broadinstitute.hellbender.utils.clustering;

import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
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
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.NaturalLogUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;

public final class BayesianGaussianMixtureModeller implements Serializable {

    private static final long serialVersionUID = 1L;

    public enum InitMethod {
        K_MEANS_PLUS_PLUS, RANDOM, TEST
    }

    private static final Logger logger = LogManager.getLogger(BayesianGaussianMixtureModeller.class);

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
    private RealMatrix covariancePrior;       // not final; if not specified via Builder and set to null, will then be set to the covariance of X upon fit(X)
    private final int seed;
    private final boolean warmStart;
    private final int verboseInterval;
    private final double relativeSymmetryThreshold;
    private final double absolutePositivityThreshold;
    private final double epsilon;

    private boolean isConverged;
    private boolean isFitAvailable;
    private double lowerBound;
    private int bestInit;

    private RealVector weightConcentration;
    private RealVector meanPrecision;
    private List<RealVector> means;
    private List<RealMatrix> precisionsCholesky;
    private List<RealMatrix> covariances;
    private RealVector degreesOfFreedom;

    private BayesianGaussianMixtureModelPosterior bestFit;

    /**
     * We use a {@link Builder} to enable the setting of default values, which can be
     * done trivially in python. Our implementation of parameter validation thus slightly differs from that in sklearn.
     */
    private BayesianGaussianMixtureModeller(final int nComponents,
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
                                            final int verboseInterval,
                                            final double relativeSymmetryThreshold,
                                            final double absolutePositivityThreshold,
                                            final double epsilon) {
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
        this.relativeSymmetryThreshold = relativeSymmetryThreshold;
        this.absolutePositivityThreshold = absolutePositivityThreshold;
        this.epsilon = epsilon;

        isConverged = false;
        isFitAvailable = false;
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

        BayesianGaussianMixtureUtils.logHeapUsage(logger, "Starting loop...");
        for (int init = 0; init < nInit; init++) {
            BayesianGaussianMixtureUtils.logHeapUsage(logger, "Starting initialization...");
            logger.info(String.format("Initialization %d...", init));

            if (doInit) {
                initializeParameters(X, rng);
            }

            double lowerBound = doInit ? Double.NEGATIVE_INFINITY : this.lowerBound;

            BayesianGaussianMixtureUtils.logHeapUsage(logger, "Starting iteration...");


            for (int nIter = 1; nIter <= maxIter; nIter++) {
                final double prevLowerBound = lowerBound;
                final RealVector prevWeights = weightConcentration.copy().mapDivideToSelf(BayesianGaussianMixtureUtils.sum(weightConcentration));

                final RealMatrix logResp = EStep(X);
                MStep(X, logResp);
                lowerBound = computeLowerBound(logResp);

                final RealVector weights = weightConcentration.copy().mapDivideToSelf(BayesianGaussianMixtureUtils.sum(weightConcentration));
                final double totalAbsoluteWeightsChange = IntStream.range(0, nComponents)
                        .mapToDouble(k -> Math.abs(prevWeights.getEntry(k) - weights.getEntry(k)))
                        .sum();
                prevWeights.combineToSelf(0, 1, weights);

                final double lowerBoundChange = lowerBound - prevLowerBound;

                if (nIter % verboseInterval == 0) {
                    logger.info(String.format("Iteration %d, lower bound = %.5f, lower-bound change = %.5f, total absolute weights change = %.5f...",
                            nIter, lowerBound, lowerBoundChange, totalAbsoluteWeightsChange));
                }

                if (Math.abs(lowerBoundChange) < tol) {
                    isConverged = true;
                    logger.info(String.format("Initialization %d converged after %d iterations, final lower bound = %.5f, final lower-bound change = %.5f, total absolute weights change = %.5f...",
                            init, nIter, lowerBound, lowerBoundChange, totalAbsoluteWeightsChange));
                    break;
                }

                if (nIter == maxIter && !isConverged) {
                    logger.info(String.format("Initialization %d did not converge after %d iterations, final lower bound = %.5f, final lower-bound change = %.5f, total absolute weights change = %.5f...",
                            init, nIter, lowerBound, lowerBoundChange, totalAbsoluteWeightsChange));
                }
            }

            if (lowerBound > maxLowerBound || maxLowerBound == Double.NEGATIVE_INFINITY) {
                maxLowerBound = lowerBound;
                logger.info(String.format("New maximum lower bound = %.5f found with initialization %d...",
                        maxLowerBound, init));
                bestInit = init;
                bestFit = new BayesianGaussianMixtureModelPosterior(
                        weightConcentration,
                        meanPrecision,
                        means,
                        precisionsCholesky,
                        covariances,
                        degreesOfFreedom);
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
        setCurrentAndBestFits(bestFit);
        isFitAvailable = true;
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
        if (!isFitAvailable) {
            throw new UnsupportedOperationException("Cannot score samples before model has been fit or set.");
        }
        final RealMatrix X = new Array2DRowRealMatrix(data);
        final int nSamples = X.getRowDimension();
        final RealMatrix weightedLogProb = estimateWeightedLogProb(X);
        return IntStream.range(0, nSamples).mapToDouble(                                // logsumexp(weightedLogProb, axis=1)
                i -> NaturalLogUtils.logSumExp(weightedLogProb.getRow(i))).toArray();
    }

    /**
     * @param X         data. (nSamples, nFeatures)
     * @param rng       random number generator, which should be reset to this.seed at start of fit
     */
    private void initializeParameters(final RealMatrix X,
                                      final RandomGenerator rng) {
        final int nSamples = X.getRowDimension();

        if (initMethod == InitMethod.K_MEANS_PLUS_PLUS) {
            final RealMatrix resp = BayesianGaussianMixtureUtils.calculateRespKMeansPlusPlus(X, nComponents, rng);
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
                BayesianGaussianMixtureUtils.estimateGaussianParameters(X, resp, regCovar, epsilon);
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
                BayesianGaussianMixtureUtils.estimateGaussianParameters(X, BayesianGaussianMixtureUtils.map(logResp, FastMath::exp), regCovar, epsilon);
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
        precisionsCholesky = BayesianGaussianMixtureUtils.computePrecisionCholesky(
                covariances, relativeSymmetryThreshold, absolutePositivityThreshold);
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
        return weightConcentration.map(Gamma::digamma).mapSubtract(Gamma.digamma(BayesianGaussianMixtureUtils.sum(weightConcentration)));
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
        final RealMatrix result = BayesianGaussianMixtureUtils.estimateLogGaussianProb(X, means, precisionsCholesky);                           // result = estimateLogGaussianProb; we will update result in place
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
        final RealVector logDetPrecisionsChol = BayesianGaussianMixtureUtils.computeLogDetCholesky(precisionsCholesky);
        logDetPrecisionsChol.combineToSelf(1., -0.5 * nFeatures, degreesOfFreedom.map(FastMath::log));      // update logDetPrecisionsChol in place

        final double logWishart = BayesianGaussianMixtureUtils.sum(BayesianGaussianMixtureUtils.logWishartNorm(degreesOfFreedom, logDetPrecisionsChol, nFeatures));
        final double logNormWeight = BayesianGaussianMixtureUtils.logDirichletNorm(weightConcentration);

        return -BayesianGaussianMixtureUtils.sum(BayesianGaussianMixtureUtils.ebeMultiply(BayesianGaussianMixtureUtils.map(logResp, FastMath::exp), logResp))
                - logWishart
                - logNormWeight
                - 0.5 * nFeatures * BayesianGaussianMixtureUtils.sum(meanPrecision.map(FastMath::log));
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
                            i -> BayesianGaussianMixtureUtils.sum(X.getColumnVector(i)) / nSamples)
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
            covariancePrior = new Covariance(X).getCovarianceMatrix();      // note Covariance assumes columns correspond to features; np.cov assumes rows correspond to features
        } else {
            BayesianGaussianMixtureUtils.checkCovariance(covariancePrior, relativeSymmetryThreshold, absolutePositivityThreshold);
            Utils.validateArg(covariancePrior.getRowDimension() == nFeatures,
                    String.format("covariancePrior and the data should correspond to the same number of features, " +
                            "but covariancePrior has dimensions %d x %d and the data has %d features.",
                            covariancePrior.getRowDimension(), covariancePrior.getColumnDimension(), nFeatures));
        }
    }

    //******************************************************************************************************************
    // getters and toString
    //******************************************************************************************************************

    public double getLowerBound() {
        if (lowerBound == Double.NEGATIVE_INFINITY) {
            throw new UnsupportedOperationException("Lower bound has not yet been determined. Call a fitting method first.");
        }
        return lowerBound;
    }

    public BayesianGaussianMixtureModelPosterior getBestFit() {
        if (!isFitAvailable) {
            throw new UnsupportedOperationException("Best fit has not yet been determined. Call a fitting method or set the fits first.");
        }
        return bestFit;
    }

    /**
     * Should only be used to set fits in preparation for calling {@link #scoreSamples}.
     * To be consistent with sklearn, setting fits will not affect the result of calling {@link #fit} or {@link #fitPredict},
     * even if {@link #warmStart} is true; the fits will essentially be overwritten by the new ones.
     */
    public void setCurrentAndBestFits(final BayesianGaussianMixtureModelPosterior fit) {
        weightConcentration = fit.getWeightConcentration();
        meanPrecision = fit.getMeanPrecision();
        means = fit.getMeans();
        precisionsCholesky = fit.getPrecisionsCholesky();
        covariances = fit.getCovariances();
        degreesOfFreedom = fit.getDegreesOfFreedom();
        bestFit = fit;
        isFitAvailable = true;
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
                ", relativeSymmetryThreshold=" + relativeSymmetryThreshold +
                ", absolutePositivityThreshold=" + absolutePositivityThreshold +
                ", epsilon=" + epsilon +
                '}';
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
        private RealVector meanPrior = null;                    // if null, will be set to the mean of X (over samples, i.e., over rows) upon fit(X)
        private Double degreesOfFreedomPrior = null;            // if null, will be set to nFeatures (i.e., X.getColumnDimension()) upon fit(X)
        private RealMatrix covariancePrior;                     // if null, will be set to the covariance of X upon fit(X)

        private int seed = 0;                                   // defaults to None (i.e., uses the global random state) in the sklearn implementation; we disallow this to avoid reproducibility headaches
        private boolean warmStart = false;
        private int verboseInterval = 10;

        private double relativeSymmetryThreshold = 1E-6;
        private double absolutePositivityThreshold = 1E-10;
        private double epsilon = 1E-10;

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

        public Builder weightConcentrationPrior(final Double weightConcentrationPrior) {
            if (weightConcentrationPrior != null) {
                Utils.validateArg(weightConcentrationPrior > 0., "weightConcentrationPrior must be > 0.");
            }
            this.weightConcentrationPrior = weightConcentrationPrior;
            return this;
        }

        public Builder meanPrecisionPrior(final double meanPrecisionPrior) {
            Utils.validateArg(meanPrecisionPrior > 0., "meanPrecisionPrior must be > 0.");
            this.meanPrecisionPrior = meanPrecisionPrior;
            return this;
        }

        public Builder meanPrior(final double[] meanPrior) {
            if (meanPrior != null) {
                this.meanPrior = new ArrayRealVector(meanPrior, true);
            } else {
                this.meanPrior = null;
            }
            return this;
        }

        public Builder degreesOfFreedomPrior(final Double degreesOfFreedomPrior) {
            if (degreesOfFreedomPrior != null) {
                Utils.validateArg(degreesOfFreedomPrior > 0., "degreesOfFreedomPrior must be > 0.");
            }
            this.degreesOfFreedomPrior = degreesOfFreedomPrior;
            return this;
        }

        public Builder covariancePrior(final double[][] covariancePrior) {
            if (covariancePrior != null) {
                this.covariancePrior = new Array2DRowRealMatrix(covariancePrior, true);
                BayesianGaussianMixtureUtils.checkCovariance(this.covariancePrior, relativeSymmetryThreshold, absolutePositivityThreshold);
            } else {
                this.covariancePrior = null;
            }
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

        public Builder relativeSymmetryThreshold(final double relativeSymmetryThreshold) {
            Utils.validateArg(relativeSymmetryThreshold >= 0., "relativeSymmetryThreshold must be >= 0.");
            this.relativeSymmetryThreshold = relativeSymmetryThreshold;
            return this;
        }

        public Builder absolutePositivityThreshold(final double absolutePositivityThreshold) {
            Utils.validateArg(absolutePositivityThreshold >= 0., "absolutePositivityThreshold must be >= 0.");
            this.absolutePositivityThreshold = absolutePositivityThreshold;
            return this;
        }

        public Builder epsilon(final double epsilon) {
            Utils.validateArg(epsilon >= 0., "absolutePositivityThreshold must be >= 0.");
            this.epsilon = epsilon;
            return this;
        }

        public BayesianGaussianMixtureModeller build() {
            weightConcentrationPrior = weightConcentrationPrior == null ? 1. / nComponents : weightConcentrationPrior;
            return new BayesianGaussianMixtureModeller(
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
                    verboseInterval,
                    relativeSymmetryThreshold,
                    absolutePositivityThreshold,
                    epsilon);
        }
    }
}