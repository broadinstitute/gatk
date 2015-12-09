package org.broadinstitute.hellbender.utils.pca;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hdf5.HDF5File;
import org.broadinstitute.hellbender.utils.svd.SVD;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.function.Function;
import java.util.stream.DoubleStream;

/**
 * Principal component analysis.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class PCA {

    public static final String PCA_GROUP_NAME = "PCA";

    public static final String VARIANCES_FULL_PATH = PCA_GROUP_NAME + "/" + "VARIANCES";

    public static final String SAMPLES_FULL_PATH = PCA_GROUP_NAME + "/" + "SAMPLES";

    public static final String VARIABLES_FULL_PATH = PCA_GROUP_NAME + "/" + "VARIABLES";

    public static final String CENTERS_FULL_PATH = PCA_GROUP_NAME + "/" + "CENTERS";

    public static final String EIGEN_VECTORS_FULL_PATH = PCA_GROUP_NAME + "/" + "EIGEN_VECTORS";

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
     * List of variable names as they appear in the input data-matrix (one per row).
     */
    private final List<String> variables;

    /**
     * List of samples as they appear in the input data-matrix (one per column).
     */
    private final List<String> samples;

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
     * @param dataMatrix the input data matrix. It can be {@code null}.
     * @param svdFactory the factory to create a SVD from a data-matrix. It can be {@code null}.
     *
     * @throws IllegalArgumentException if either {@code dataMatrix} or {@code svdFactory} is {@code null}.
     */
    public static PCA createPCA(final RealMatrix dataMatrix, final Function<RealMatrix, SVD> svdFactory) {
        return createPCA(null, null, dataMatrix, svdFactory);
    }

    /**
     * Creates a new PCA result using SVD.
     * <p>
     * This operation will do all required computation, thus it might take long to complete for
     * large matrices.
     * </p>
     * <p>
     * The input matrix rows must represent variables whereas columns represent samples.
     * </p>
     *
     * @param variables  variable names following the same order as the row in the data-matrix.
     * @param samples    sample names following the same order as the columns in the data-matrix.
     * @param dataMatrix the input data matrix. It can be {@code null}.
     * @param svdFactory the factory to create a SVD from a data-matrix. It can be {@code null}.
     * @throws IllegalArgumentException if either {@code dataMatrix} or {@code svdFactory} is {@code null}. If
     *                                  either {@code variables} or {@code samples} is non-null and inconsistent with the input data matrix
     *                                  (non matching number of rows for variables, columns for samples), contain a {@code null} or contain
     *                                  repeated names.
     */
    public static PCA createPCA(final List<String> variables, final List<String> samples,
                                final RealMatrix dataMatrix, final Function<RealMatrix, SVD> svdFactory) {
        Utils.nonNull(dataMatrix, "the input matrix cannot be null");
        Utils.nonNull(svdFactory, "the SVD factory cannot be null");
        final List<String> variableNames = checkNonNullUniqueNames(variables, "variable");
        final List<String> sampleNames = checkNonNullUniqueNames(samples, "sample");

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
        return new PCA(variableNames, sampleNames, centers, eigenVectors, variances);
    }

    private static List<String> checkNonNullUniqueNames(final List<String> input, final String roleName) {
        if (input == null) {
            return null;
        } else if (input.stream().anyMatch(Objects::isNull)) {
            throw new IllegalArgumentException(String.format("the input %s list must not contain nulls", roleName));
        } else if (input.stream().sorted().distinct().count() != input.size()) {
            throw new IllegalArgumentException(String.format("the input %s list must not contain repeats", roleName));
        } else {
            return Collections.unmodifiableList(new ArrayList<>(input));
        }

    }

    /**
     * Creates a PCA instance given all its member values.
     *
     * @param centers the new instance centers.
     * @param eigenVectors the new instance eigenVectors.
     * @param variances the new instance variances.
     */
    private PCA(final List<String> variables, final List<String> samples, final double[] centers, final RealMatrix eigenVectors, final double[] variances) {
        this.variables = variables;
        this.samples = samples;
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

    /**
     * Writes the content of the PCA into HDF5 file.
     *
     * <p>
     *     The data is stored under the group {@value #PCA_GROUP_NAME}.
     * </p>
     *
     * @param pca the PCA result to store.
     * @param output the output HDF5 file.
     */
    public static void writeHDF5(final PCA pca, final HDF5File output) {
        if (pca.variables != null ) {
            output.makeStringArray(VARIABLES_FULL_PATH, pca.variables.toArray(new String[pca.variables.size()]));
        }
        if (pca.samples != null) {
            output.makeStringArray(SAMPLES_FULL_PATH, pca.samples.toArray(new String[pca.samples.size()]));
        }
        output.makeDoubleArray(VARIANCES_FULL_PATH, pca.variances);
        output.makeDoubleArray(CENTERS_FULL_PATH, pca.centers);
        output.makeDoubleMatrix(EIGEN_VECTORS_FULL_PATH, pca.eigenVectors.getData());
    }
}
