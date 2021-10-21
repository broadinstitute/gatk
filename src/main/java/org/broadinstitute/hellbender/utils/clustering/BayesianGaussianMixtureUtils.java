package org.broadinstitute.hellbender.utils.clustering;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.NonPositiveDefiniteMatrixException;
import org.apache.commons.math3.linear.NonSquareMatrixException;
import org.apache.commons.math3.linear.NonSymmetricMatrixException;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.IntStream;

final class BayesianGaussianMixtureUtils {

    private static final double LOG_2_PI = Math.log(2. * Math.PI);

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
     * @param covariances                   covariance matrices for all components. List with length nComponents of (nFeatures, nFeatures) matrices
     * @param relativeSymmetryThreshold     threshold above which off-diagonal elements are considered too different and matrix not symmetric
     * @param absolutePositivityThreshold   threshold below which diagonal elements are considered null and matrix not positive definite
     * @return                              Cholesky decompositions of precision matrices for all components. List with length nComponents of (nFeatures, nFeatures) matrices
     */
    @VisibleForTesting
    static List<RealMatrix> computePrecisionCholesky(final List<RealMatrix> covariances,
                                                     final double relativeSymmetryThreshold,
                                                     final double absolutePositivityThreshold) {
        final int nComponents = covariances.size();
        final int nFeatures = covariances.get(0).getRowDimension();

        final List<RealMatrix> precisionsChol = new ArrayList<>(Collections.nCopies(nComponents, null));
        for (int k = 0; k < nComponents; k++) {
            final RealMatrix cov = covariances.get(k);
            try {
                final RealMatrix covChol = new CholeskyDecomposition(       // covChol = np.linalg.cholesky(cov, lower=True)
                        cov, relativeSymmetryThreshold, absolutePositivityThreshold).getL();
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
     * @param epsilon               epsilon added to effective number; prevents divide-by-zero caused by empty components
     * @return                      triple of effective number for all components, mean vectors for all components, and covariance matrices for all components.
     *                              i.e., Triple of [(nComponents, ), List with length nComponents of (nFeatures, ) vectors, List with length nComponents of (nFeatures, nFeatures) matrices]
     */
    @VisibleForTesting
    static Triple<RealVector, List<RealVector>, List<RealMatrix>> estimateGaussianParameters(final RealMatrix X,
                                                                                             final RealMatrix resp,
                                                                                             final double regCovar,
                                                                                             final double epsilon) {
        final int nComponents = resp.getColumnDimension();

        final RealVector nk = new ArrayRealVector(
                IntStream.range(0, nComponents).mapToDouble(
                        k -> sum(resp.getColumnVector(k)) + epsilon)
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
    static List<RealMatrix> estimateGaussianCovariances(final RealMatrix X,
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

    static double sum(final RealVector v) {
        return Arrays.stream(v.toArray()).sum();
    }

    static double sum(final RealMatrix m) {
        double sum = 0.;
        for (int i = 0; i < m.getRowDimension(); i++) {
            for (int j = 0; j < m.getColumnDimension(); j++) {
                sum += m.getEntry(i, j);
            }
        }
        return sum;
    }

    static RealVector diag(final RealMatrix m) { // m is assumed to be a square matrix
        final int dim = m.getRowDimension();
        final RealVector diag = new ArrayRealVector(dim);
        IntStream.range(0, dim).forEach(i -> diag.setEntry(i, m.getEntry(i, i)));
        return diag;
    }

    static RealMatrix map(final RealMatrix m,
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

    static RealMatrix ebeMultiply(final RealMatrix m,
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

    static void checkCovariance(final RealMatrix covariance,
                                final double relativeSymmetryThreshold,
                                final double absolutePositivityThreshold) {
        try {
            new CholeskyDecomposition(covariance, relativeSymmetryThreshold, absolutePositivityThreshold);
        } catch (final NonSquareMatrixException | NonSymmetricMatrixException | NonPositiveDefiniteMatrixException e) {
            throw new IllegalArgumentException("Covariance matrix must be square, symmetric, and positive semidefinite.", e);
        }
    }

    static void logHeapUsage(final Logger logger,
                             final String message) {
        final int mb = 1024 * 1024;
        final Runtime runtime = Runtime.getRuntime();
        logger.debug("Used memory [MB]: " + (runtime.totalMemory() - runtime.freeMemory()) / mb);
        logger.debug(message);
    }
}
