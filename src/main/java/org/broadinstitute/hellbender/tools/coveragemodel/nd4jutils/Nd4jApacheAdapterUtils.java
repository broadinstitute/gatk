package org.broadinstitute.hellbender.tools.coveragemodel.nd4jutils;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.broadinstitute.hellbender.utils.Utils;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import javax.annotation.Nonnull;

/**
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class Nd4jApacheAdapterUtils {
    /**
     * INDArray to Apache
     *
     * @param matrix rank-2 INDArray
     * @return Apache matrix
     */
    public static RealMatrix convertINDArrayToApacheMatrix(@Nonnull final INDArray matrix) {
        Utils.validateArg(matrix.rank() == 2, "Input rank is not 2 (not matrix)");
        final int[] shape = matrix.shape();
        final INDArray concreteMatrix = matrix.isView() ? matrix.dup() : matrix;
        final double[] data = concreteMatrix.data().asDouble();
        final char ordering = concreteMatrix.ordering();
        if (ordering == 'c') {
            return new BlockRealMatrix(monoToBiDiArrayRowMajor(data, shape[0], shape[1]));
        } else { /* ordering == 'f' */
            return new BlockRealMatrix(monoToBiDiArrayColumnMajor(data, shape[0], shape[1]));
        }
    }

    /**
     * INDArray to Apache
     *
     * @param vector rank-1 INDArray
     * @return Apache vector
     */
    public static RealVector convertINDArrayToApacheVector(@Nonnull final INDArray vector) {
        Utils.validateArg(vector.isVector(), "Input INDArray is not a vector.");
        return new ArrayRealVector(vector.isView()
                ? vector.dup().data().asDouble()
                : vector.data().asDouble(), false);
    }

    /**
     * Apache to INDArray
     *
     * @param matrix Apache matrix
     * @return rank-2 INDArray
     */
    public static INDArray convertApacheMatrixToINDArray(@Nonnull final RealMatrix matrix) {
        return Nd4j.create(matrix.getData());
    }

    /**
     * Converts a 1D double array to a 2D double array with given shape (in c order = row-major)
     *
     * @param array a 1D double array
     * @param rows number of rows in the output 2D array
     * @param cols number of columns in the output 2D array
     * @return a 2D double array
     */
    private static double[][] monoToBiDiArrayRowMajor(final double[] array, final int rows, final int cols) {
        Utils.validateArg(array.length == rows * cols, "Invalid array length");
        double[][] bidi = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            System.arraycopy(array, (i * cols), bidi[i], 0, cols);
        }
        return bidi;
    }

    /**
     * Converts a 1D double array to a 2D double array with given shape (in Fortran order = column-major)
     *
     * @param array a 1D double array
     * @param rows number of rows in the output 2D array
     * @param cols number of columns in the output 2D array
     * @return a 2D double array
     */
    private static double[][] monoToBiDiArrayColumnMajor(final double[] array, final int rows, final int cols) {
        return transposeBiDiArray(monoToBiDiArrayRowMajor(array, cols, rows));
    }

    /**
     * Transpose a 2D double array
     *
     * @param array a 2D double array
     * @return transposed 2D double array
     */
    private static double[][] transposeBiDiArray(final double[][] array) {
        final int rows = array.length;
        Utils.validateArg(rows > 0, "Array must be non-empty");
        final int cols = array[0].length;
        final double[][] trans = new double[cols][rows];
        for (int i = 0; i < cols; i++) {
            for (int j = 0; j < rows; j++) {
                trans[i][j] = array[j][i];
            }
        }
        return trans;
    }
}
