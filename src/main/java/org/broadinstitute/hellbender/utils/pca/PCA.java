package org.broadinstitute.hellbender.utils.pca;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.svd.SVD;

import java.util.function.Function;
import java.util.stream.DoubleStream;

/**
 * Principal component analysis.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class PCA {

    /**
     * Array of variable centers subtracted from the data-matrix before PCA.
     */
    private final double[] centers;

    /**
     * Principal component directions or eigenVectors, one per column in the same order
     * as in {@link #variances}.
     */
    private final RealMatrix eigenVectors;

    /**
     * Principal component variances sorted by magnitude (large variance comes first).
     */
    private final double[] variances;

    /**
     * Creates a new PCA result using SVD.
     *
     * <p>
     *     This operation will do all required computation, thus it might take long to complete for
     *     large matrices.
     * </p>
     *
     * <p>
     *     The input matrix rows must represent variables whereas columns represent samples.
     * </p>
     *
     * @param dataMatrix the input data matrix.
     * @param svdFactory the factory to create a SVD from a data-matrix.
     *
     * @throws IllegalArgumentException if either {@code dataMatrix} or {@code svdFactory} is {@code null}.
     */
    public static PCA createPCA(final RealMatrix dataMatrix, final Function<RealMatrix, SVD> svdFactory) {
        Utils.nonNull(dataMatrix, "the input matrix cannot be null");
        Utils.nonNull(svdFactory, "the SVD factory cannot be null");
        final int rowCount = dataMatrix.getRowDimension();
        final int columnCount = dataMatrix.getColumnDimension();
        final double[] centers = new double[rowCount];
        final RealMatrix centered =
                new Array2DRowRealMatrix(rowCount, columnCount);
        for (int i = 0; i < rowCount; i++) {
            final double[] row = dataMatrix.getRow(i);
            final double center = MathUtils.mean(row, 0, row.length);
            for (int j = 0; j < columnCount; j++) {
                row[j] -= center;
            }
            centered.setRow(i, row);
            centers[i] = center;
        }
        final SVD svd = svdFactory.apply(centered);
        final RealMatrix eigenVectors = svd.getU();
        final double inverseDenominator = 1.0 / (columnCount - 1.0);
        final double[] variances = DoubleStream.of(svd.getSingularValues()).map(d -> d * d * inverseDenominator).toArray();
        return new PCA(centers, eigenVectors, variances);
    }

    /**
     * Creates a PCA instance given all its member values.
     *
     * @param centers the new instance centers.
     * @param eigenVectors the new instance eigenVectors.
     * @param variances the new instance variances.
     */
    private PCA(final double[] centers, final RealMatrix eigenVectors, final double[] variances) {
        this.centers = centers;
        this.eigenVectors = eigenVectors;
        this.variances = variances;
    }

    /**
     * Returns the eigen-vectors for the principal components.
     * <p>
     *     The result matrix has one column per principal components.
     * </p>
     * <p>
     *     Each row represent the contributions of the corresponding input variable to that component.
     * </p>
     *
     * @return never {@code null}.
     */
    public RealMatrix getEigenVectors() {
        return eigenVectors;
    }

    /**
     * Returns the variable centers subtracted from the data before performing PCA.
     * @return never {@code null}.
     */
    public RealVector getCenters() {
        return new ArrayRealVector(centers);
    }

    /**
     * Returns the variances of each principal component.
     * @return never {@code null}.
     */
    public RealVector getVariances() {
        return new ArrayRealVector(variances);
    }
}
