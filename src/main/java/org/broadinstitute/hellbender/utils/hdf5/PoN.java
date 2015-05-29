package org.broadinstitute.hellbender.utils.hdf5;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;

/**
 * Panel of Normal data structure.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public interface PoN {

    /**
     * Target names listed by their numerical index in this PoN.
     *
     * <p>
     * The returned list cannot be modified and its content might change to reflect
     * changes in the underlying PoN storing resource depending on the implementation.
     * </p>
     *
     * @return never {@code null}.
     */
    List<String> targetNames();

    /**
     * Sample names listed by their numerical index in this PoN.
     *
     * <p>
     * The returned list cannot be modified and its content might change to reflect
     * changes in the underlying PoN storing resource depending on the implementation.
     * </p>
     *
     * @return never {@code null}.
     */
    List<String> sampleNames();

    /**
     * Log normal selected samples listed by their numerical index in the log-normal matrices.
     *
     * <p>
     * The return matrix is a modifiable detached copy of the values in the
     * PoN.
     * </p>
     *
     * @return never {@code null}.
     */
    List<String> logNormalSampleNames();

    /**
     * Target factors.
     *
     * <p>
     * The result matrix is a column vector where the ith element is the factor for the ith target
     * as returned by {@link #targetNames}.
     * </p>
     *
     * <p>
     * The return matrix is a modifiable detached copy of the values in the
     * PoN.
     * </p>
     *
     * @return never {@code null}. A matrix with {@code Tx1} dimensions where {@code T}
     * is the number of targets in this PoN.
     */
    RealMatrix targetFactors();

    /**
     * Normalized percent coverage.
     *
     * <p>
     * The result matrix {@code [i,j]} element is the normalized percent coverage of target ith in
     * sample jth.
     * </p>
     *
     * <p>
     * The return matrix is a modifiable detached copy of the values in the
     * PoN.
     * </p>
     *
     * @return never {@code null}. A matrix with {@code TxS} dimensions where {@code T}
     * is the number of targets and {@code S} the number of samples in this PoN.
     */
    RealMatrix normalizedPercentCoverage();

    /**
     * Returns the maximum ratio cut-off.
     *
     * @return a valid ratio cut-off.
     */
    double maximumRatioCutoff();

    /**
     * Returns the PoN version.
     *
     * @return never {@code null}, an string of the form ver.sub-ver
     */
    String version();

    /**
     * Returns the log-normal matrix.
     *
     * <p>
     * The return matrix is a modifiable detached copy of the values in the
     * PoN.
     * </p>
     *
     * <p>
     * The return matrix is a modifiable detached copy of the values in the
     * PoN.
     * </p>
     *
     * @return never {@code null}, a matrix with dimensions {@code SxT} where {@code T} is the number of targets and
     * {@code S} the number of samples considered for the log-normal.
     */
    RealMatrix logNormals();

    /**
     * Returns the log-normal pseudo-inverse matrix.
     *
     * <p>
     * The return matrix is a modifiable detached copy of the values in the
     * PoN.
     * </p>
     *
     * @return never {@code null}, a matrix with dimensions {@code TxS} where {@code T} is the number of targets and
     * {@code S} the number of samples considered for the log-normal.
     */
    RealMatrix logNormalsPseudoInverse();

    /**
     * Returns the reduced PoN matrix.
     *
     * <p>
     * The return matrix is a modifiable detached copy of the values in the
     * PoN.
     * </p>
     *
     * @return never {@code null}, a matrix with dimensions {@code TxE} where {@code T} is the number of targets and
     * {@code E} the number of eigen samples.
     */
    RealMatrix reducedPoN();

    /**
     * Returns the reduced PoN pseudo-inverse matrix.
     *
     * <p>
     * The return matrix is a modifiable detached copy of the values in the
     * PoN.
     * </p>
     *
     * @return never {@code null}, a matrix with dimensions {@code ExT} where {@code T} is the number of targets and
     * {@code E} the number of eigen samples.
     */
    RealMatrix reducedPoNPseudoInverse();

    /**
     * Performs target factor normalization on a matrix.
     *
     * <p>
     *     The input matrix remains unchanged.
     * </p>
     *
     * <p>
     *     The result matrix is a new instance that can be further modified by the calling code.
     * </p>
     */
    default RealMatrix factorNormalization(final RealMatrix input) {
        Utils.nonNull(input, "the input matrix cannot be null");
        final RealMatrix targetFactors = targetFactors();
        if (input.getRowDimension() != targetFactors.getRowDimension()) {
            throw new IllegalArgumentException(String.format("the input row count must match the PoN's target count: %d != %d",
                    input.getRowDimension(),targetFactors.getRowDimension()));
        }
        final RealMatrix result = new Array2DRowRealMatrix(input.getRowDimension(),input.getColumnDimension());
        for (int i = 0; i < result.getRowDimension(); i++) {
            final double[] values = input.getRow(i);
            final double factor = targetFactors.getEntry(i, 0);
            final double inverseFactor = 1.0 / factor;
            for (int j = 0; j < values.length; j++) {
                values[j] *= inverseFactor;
            }
            result.setRow(i,values);
        }
        return result;
    }

}