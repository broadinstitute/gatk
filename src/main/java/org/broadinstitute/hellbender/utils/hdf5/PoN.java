package org.broadinstitute.hellbender.utils.hdf5;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;
import java.util.stream.IntStream;

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
    List<String> getTargetNames();

    /**
     * Reduced PoN Target names listed by their numerical index in the reduced PoN.
     *
     * <p>
     * The returned list cannot be modified and its content might change to reflect
     * changes in the underlying PoN storing resource depending on the implementation.
     * </p>
     *
     * @return never {@code null}.
     */
    List<String> getPanelTargetNames();

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
    List<String> getSampleNames();

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
    List<String> getPanelSampleNames();

    /**
     * Return the target factors.
     *
     * <p>
     * The result matrix is a column vector where the ith element is the factor for the ith target
     * as returned by {@link #getTargetNames}.
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
    RealMatrix getTargetFactors();

    /**
     * Set the target factors in the pon.
     *
     * @param targetFactors the new value for the target factors.
     * @throws IllegalArgumentException if {@code targetFactors} is null or has the wrong
     * dimensions.
     */
    void setTargetFactors(final RealMatrix targetFactors);

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
    RealMatrix getNormalizedCounts();

    /**
     * Returns the PoN version.
     *
     * <p>
     *     The version major is the integer part, whereas the version minor is the
     *     decimal part.
     * </p>
     *
     * @return any valid double value.
     */
    double getVersion();

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
    RealMatrix getLogNormalizedCounts();

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
    RealMatrix getLogNormalizedPInverseCounts();

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
    RealMatrix getReducedPanelCounts();

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
    RealMatrix getReducedPanelPInverseCounts();

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
        final RealMatrix targetFactors = getTargetFactors();
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
            result.setRow(i, values);
        }
        return result;
    }

    /**
     * Calculate the beta-hats that best fit case read counts given the panel of normals.
     *
     * @param input a {@code TxS} matrix where {@code T} is the number of targets and {@code S} the number of count groups (e.g. case samples).
     * @param useReduced whether to use the reduced PoN or the original one.
     * @return never {@code null} a matrix of {@code NxS} dimensions, where N is the number of samples in
     *  the panel and S the original name of count groups.
     */
    default RealMatrix betaHats(final RealMatrix input, final boolean useReduced, final double epsilon) {
        Utils.nonNull(input,"the input counts must not be null");
        if (Double.isNaN(epsilon) || epsilon <= 0) {
            throw new IllegalArgumentException(String.format("invalid epsilon value must be > 0: %f", epsilon));
        }
        final double targetThreshold = (Math.log(epsilon) / Math.log(2)) + 1;
        final RealMatrix normalsInverse = useReduced ? getReducedPanelPInverseCounts() : getLogNormalizedPInverseCounts();

        // copy case samples in order to mask targets in-place and mask (set to zero) targets with coverage below threshold
        final RealMatrix maskedInput = input.copy();
        maskedInput.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int row, final int column, final double value) { return value > targetThreshold ? value : 0; }
        });

        return normalsInverse.multiply(maskedInput);
    }

    /**
     * Applies tangent normalization.
     * <p>
     *     The input row order is suppose to match the panels target order.
     * </p>
     *
     * @param input the input counts to normalize. This matrix dimension is TxS where T is the number of targets
     *              and S the number of count columns.
     * @param betaHats the beta-hats for the projection to use for the normalization. The dimension of this matrix
     *                  is NxS where N is the number of samples in the panel of choice and S is the number of count columns.
     * @param useReduced whether to use the reduced or the full normal panel sample.
     * @return never {@code null}.
     */
    default RealMatrix tangentNormalization(final RealMatrix input, final RealMatrix betaHats, final boolean useReduced) {
        Utils.nonNull(input,"tangent normalization input counts cannot be null");
        Utils.nonNull(betaHats,"tangent beta-hats cannot be null");
        if (input.getColumnDimension() != betaHats.getColumnDimension()) {
            throw new IllegalArgumentException(String.format("the input count column count (%d) does not match the number of columns in the beta-hats (%d)", input.getColumnDimension(), betaHats.getColumnDimension()));
        }
        final RealMatrix normals = useReduced ? getReducedPanelCounts() : getLogNormalizedCounts();
        if (normals.getColumnDimension() != betaHats.getRowDimension()) {
            throw new IllegalArgumentException(String.format("beta-hats component count (%d) does not match the number of samples in the PoN (%d)", normals.getRowDimension(), normals.getColumnDimension()));
        }
        final RealMatrix projection = normals.multiply(betaHats);
        return input.subtract(projection);
    }
}