package org.broadinstitute.hellbender.tools.pon.coverage.pca;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.pon.coverage.CoveragePanelOfNormals;

import java.util.List;

/**
 * Interface for coverage panel of normals data structure.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public interface PCACoveragePoN extends CoveragePanelOfNormals<PCATangentNormalizationResult> {
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
     *  Get a list of the targets that are in this PoN (some filtering will have been applied).
     *
     * @return never {@code null} List of Target instances that are a copy (modifiable; not references) of what was in the PoN.
     */
    List<Target> getTargets();

   /**
     *  Get a list of the targets that were initially submitted to this PoN (no filtering will have been applied).
     *
     * @return never {@code null} List of Target instances that are a copy (modifiable; not references) of what was in the PoN.
     */
    List<Target> getRawTargets();

    /**
     *  Get a list of the targets that are in the final hyperplane after all filtering has been applied..
     *
     * @return never {@code null} List of Target instances that are a copy (modifiable; not references) of what was in the PoN.
     */
    List<Target> getPanelTargets();

    /**
     * Get an array of target factors,
     * where the ith element is the factor for the ith target as returned by {@link #getTargetNames}
     *
     * <p>
     * The returned array cannot be modified and its content might change to reflect
     * changes in the underlying PoN storing resource depending on the implementation.
     * </p>
     *
     * @return never {@code null}.
     */
    double[] getTargetFactors();

    /**
     * Get an array of target variances (post-tangent normalization) that are in the final hyperplane after all filtering has been applied,
     * where the ith element is the variance for the ith target as returned by {@link #getTargetNames}
     *
     * <p>
     * The returned array cannot be modified and its content might change to reflect
     * changes in the underlying PoN storing resource depending on the implementation.
     * </p>
     *
     * @return never {@code null}.
     */
    double[] getTargetVariances();

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
     * @return never {@code null}. A matrix with dimensions {@code TxS} where {@code T}
     * is the number of targets and {@code S} the number of samples in this PoN.
     */
    RealMatrix getNormalizedCounts();

    /**
     * Returns the log-normalized matrix.
     *
     * <p>
     * The return matrix is a modifiable detached copy of the values in the
     * PoN.
     * </p>
     *
     * @return never {@code null}, a matrix with dimensions {@code TxS} where {@code T} is the number of targets and
     * {@code S} the number of samples considered for the log-normalized panel.
     */
    RealMatrix getLogNormalizedCounts();

    /**
     * Returns the log-normalized pseudoinverse matrix.
     *
     * <p>
     * The return matrix is a modifiable detached copy of the values in the
     * PoN.
     * </p>
     *
     * @return never {@code null}, a matrix with dimensions {@code SxT} where {@code T} is the number of targets and
     * {@code S} the number of samples considered for the log-normalized panel.
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
     * {@code E} the number of eigensamples.
     */
    RealMatrix getReducedPanelCounts();

    /**
     * Returns the reduced PoN pseudoinverse matrix.
     *
     * <p>
     * The return matrix is a modifiable detached copy of the values in the
     * PoN.
     * </p>
     *
     * @return never {@code null}, a matrix with dimensions {@code ExT} where {@code T} is the number of targets and
     * {@code E} the number of eigensamples.
     */
    RealMatrix getReducedPanelPInverseCounts();

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
     * Initial target names listed by their numerical index in this PoN.
     *
     * <p>
     * The returned list cannot be modified and its content might change to reflect
     * changes in the underlying PoN storing resource depending on the implementation.
     * </p>
     *
     * @return never {@code null}.
     */
    List<String> getRawTargetNames();

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
     * Log-normalized selected samples listed by their numerical index in the log-normalized matrices.
     *
     * <p>
     * The return matrix is an unmodifiable detached copy of the values in the
     * PoN.
     * </p>
     *
     * @return never {@code null}.
     */
    List<String> getPanelSampleNames();

    @Override
    default PCATangentNormalizationResult normalize(final ReadCountCollection proportionalCoverageProfile, final JavaSparkContext ctx) {
        return PCATangentNormalizationUtils.tangentNormalize(this, proportionalCoverageProfile, true, ctx);     //doFactorNormalization = true
    }

    @Override
    default PCATangentNormalizationResult normalizeNormalsInPoN(final JavaSparkContext ctx) {
        final ReadCountCollection normals = new ReadCountCollection(getTargets(), getSampleNames(), getNormalizedCounts());
        return PCATangentNormalizationUtils.tangentNormalize(this, normals, false, ctx);                       //doFactorNormalization = false, normalizedCounts are already factor normalized
    }
}