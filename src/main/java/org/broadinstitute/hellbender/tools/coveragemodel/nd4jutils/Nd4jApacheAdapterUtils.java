package org.broadinstitute.hellbender.tools.coveragemodel.nd4jutils;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
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
     * TODO github/gatk-protected issue #853 -- this method must be optimized
     *
     * @param matrix rank-2 INDArray
     * @return Apache matrix
     */
    public static RealMatrix convertINDArrayToApacheMatrix(@Nonnull final INDArray matrix) {
        if (matrix.rank() != 2) {
            throw new IllegalArgumentException("Input rank is not 2 (not matrix)");
        }
        final int[] shape = matrix.shape();
        final BlockRealMatrix out = new BlockRealMatrix(shape[0], shape[1]);
        for (int i = 0; i < shape[0]; i++) {
            for (int j = 0; j < shape[1]; j++) {
                out.setEntry(i, j, matrix.getDouble(i, j));
            }
        }
        return out;
    }

    /**
     * INDArray to Apache
     *
     * TODO github/gatk-protected issue #853 -- this method must be optimized
     *
     * @param vector rank-1 INDArray
     * @return Apache vector
     */
    public static RealVector convertINDArrayToApacheVector(@Nonnull final INDArray vector) {
        if (!vector.isVector()) {
            throw new IllegalArgumentException("Input INDArray is not a vector.");
        }
        final int length = vector.length();
        final RealVector out = new ArrayRealVector(length);
        for (int i=0; i<length; i++) {
            out.setEntry(i, vector.getDouble(i));
        }
        return out;
    }

    /**
     * Apache to INDArray
     *
     * TODO github/gatk-protected issue #853 -- this method must be optimized
     *
     * @param matrix Apache matrix
     * @return rank-2 INDArray
     */
    public static INDArray convertApacheMatrixToINDArray(@Nonnull final RealMatrix matrix) {
        final int[] shape = new int[] {matrix.getRowDimension(), matrix.getColumnDimension()};
        final INDArray out = Nd4j.create(shape);
        for (int i = 0; i < shape[0]; i++) {
            for (int j = 0; j < shape[1]; j++) {
                out.putScalar(new int[] {i, j}, matrix.getEntry(i, j));
            }
        }
        return out;
    }
}
