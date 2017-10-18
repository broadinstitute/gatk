package org.broadinstitute.hellbender.tools.copynumber.coverage.denoising.svd;

import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.SimpleCountCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.List;

/**
 * Interface for the panel of normals (PoN) for SVD-based coverage denoising.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public interface SVDReadCountPanelOfNormals {
    /**
     * Returns the PoN version.
     */
    double getVersion();

    /**
     * Returns the number of eigensamples.
     */
    int getNumEigensamples();

    /**
     * Returns a modifiable copy of the original matrix of integer read-counts (represented as doubles) used to build the PoN
     * (no filtering will have been applied).  This matrix has has dimensions {@code M_original x N_original},
     * where {@code M_original} is the number of original intervals and {@code N_original} is the number of
     * original samples.
     */
    double[][] getOriginalReadCounts();

    /**
     * Returns a modifiable copy of the list of the original intervals that were used to build this PoN
     * (no filtering will have been applied).  This list has length {@code M_original}.
     */
    List<SimpleInterval> getOriginalIntervals();

    /**
     * Returns a modifiable copy of an array containing the GC content of the original intervals
     * (in the same order as in {@link #getOriginalIntervals()}).  This array has length {@code M_original}.
     */
    double[] getOriginalIntervalGCContent();

    /**
     * Returns a modifiable copy of the list of the intervals contained in this PoN after all filtering has been applied.
     * This list has length {@code M}.
     */
    List<SimpleInterval> getPanelIntervals();

    /**
     * Returns a modifiable copy of an array containing the median (across all samples, before filtering)
     * of the fractional coverage at each panel interval (in the same order as in {@link #getPanelIntervals()}).
     * This is used to standardize samples.  This array has length {@code M}.
     */
    double[] getPanelIntervalFractionalMedians();

    /**
     * Returns a modifiable copy of an array of the singular values of the eigensamples in decreasing order.
     * This array has length {@code K}.
     */
    double[] getSingularValues();

    /**
     * Returns a modifiable copy of an array containing the orthnonormal matrix of eigensample vectors.
     * This matrix has has dimensions {@code M x K},
     * where {@code M} is the number of panel intervals (after filtering)
     * and {@code K} is the number of eigensamples.
     * Columns are sorted by singular value in decreasing order.
     */
    double[][] getEigensampleVectors();

    default SVDDenoisedCopyRatioResult denoise(final SimpleCountCollection readCounts,
                                               final int numEigensamples) {
        return SVDDenoisingUtils.denoise(this, readCounts, numEigensamples);
    }
}