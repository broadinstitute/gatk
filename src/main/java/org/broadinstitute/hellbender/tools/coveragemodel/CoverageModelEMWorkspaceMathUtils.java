package org.broadinstitute.hellbender.tools.coveragemodel;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.util.Precision;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.coveragemodel.nd4jutils.Nd4jApacheAdapterUtils;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.ops.transforms.Transforms;

import javax.annotation.Nonnull;

/**
 * Common math utilities for {@link CoverageModelEMWorkspace}
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class CoverageModelEMWorkspaceMathUtils {

    /**
     * Inverts a square INDArray matrix using Apache LU decomposition
     *
     * @param mat matrix to be inverted
     * @return inverted matrix as an {@link INDArray}
     */
    public static INDArray minv(@Nonnull final INDArray mat) {
        if (mat.isScalar()) {
            return Nd4j.onesLike(mat).divi(mat);
        }
        if (!mat.isSquare()) {
            throw new IllegalArgumentException("Invalid array: must be square matrix");
        }
        final RealMatrix rm = Nd4jApacheAdapterUtils.convertINDArrayToApacheMatrix(mat);
        final RealMatrix rmInverse = new LUDecomposition(rm).getSolver().getInverse();
        return Nd4jApacheAdapterUtils.convertApacheMatrixToINDArray(rmInverse);
    }

    /**
     * Solves a linear system using Apache commons methods [mat].[x] = [vec]
     *
     * @param mat the coefficients matrix (must be square and full-rank)
     * @param vec the right hand side vector
     * @return solution of the linear system
     */
    public static INDArray linsolve(@Nonnull final INDArray mat, @Nonnull final INDArray vec) {
        if (mat.isScalar()) {
            return vec.div(mat.getDouble(0));
        }
        if (!mat.isSquare()) {
            throw new IllegalArgumentException("invalid array: must be a square matrix");
        }
        final RealVector sol = new LUDecomposition(Nd4jApacheAdapterUtils.convertINDArrayToApacheMatrix(mat))
                .getSolver().solve(Nd4jApacheAdapterUtils.convertINDArrayToApacheVector(vec));
        return Nd4j.create(sol.toArray(), vec.shape());
    }

    /**
     * Return the norm-infinity of {@code arr}
     * @param arr the array to be normed
     * @return norm-infinity
     */
    public static double getINDArrayNormInfinity(@Nonnull final INDArray arr) {
        return Transforms.abs(arr, true).maxNumber().doubleValue();
    }

    /**
     * Takes a square symmetric real matrix [M] and finds an orthogonal transformation [U] such that
     * [U] [M] [U]^T is diagonal, and diagonal entries are sorted in descending order.
     *
     * @param matrix a symmetric matrix
     * @param symmetrize enforce symmetry
     * @param logger a logger instance
     * @return [U]
     */
    public static ImmutablePair<double[], RealMatrix> getOrthogonalizerAndSorterTransformation(@Nonnull final RealMatrix matrix,
                                                                      final boolean symmetrize,
                                                                      @Nonnull final Logger logger) {
        if (matrix.getRowDimension() != matrix.getColumnDimension()) {
            throw new IllegalArgumentException("The input matrix must be square");
        }
        final RealMatrix finalMatrix;
        if (symmetrize) {
            final double symTol = 10 * matrix.getRowDimension() * matrix.getColumnDimension() * Precision.EPSILON;
            if (!MatrixUtils.isSymmetric(matrix, symTol)) {
                logger.info("The input matrix is not symmetric -- enforcing symmetrization");
                finalMatrix = matrix.add(matrix.transpose()).scalarMultiply(0.5);
            } else {
                finalMatrix = matrix;
            }
        } else {
            finalMatrix = matrix;
        }
        final EigenDecomposition decomposer = new EigenDecomposition(finalMatrix);
        final double[] eigs = decomposer.getRealEigenvalues();
        final RealMatrix VT = decomposer.getVT();
        return ImmutablePair.of(eigs, VT);
    }

    /**
     * Same as {@link #getOrthogonalizerAndSorterTransformation(RealMatrix, boolean, Logger)} but with
     * INDArray input/output.
     *
     * @param matrix a square symmetric real matrix
     * @param symmetrize enforce symmetrization of the matrix
     * @param logger a logger instance
     * @return an {@link INDArray} representation of the rotation operator
     */
    public static ImmutablePair<INDArray, INDArray> getOrthogonalizerAndSorterTransformation(@Nonnull final INDArray matrix,
                                                                                             final boolean symmetrize,
                                                                                             @Nonnull final Logger logger) {
        final ImmutablePair<double[], RealMatrix> out = getOrthogonalizerAndSorterTransformation(
                Nd4jApacheAdapterUtils.convertINDArrayToApacheMatrix(matrix), symmetrize, logger);
        return ImmutablePair.of(Nd4j.create(out.left, new int[] {1, matrix.shape()[0]}),
                Nd4jApacheAdapterUtils.convertApacheMatrixToINDArray(out.right));
    }
}
