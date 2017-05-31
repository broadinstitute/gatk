package org.broadinstitute.hellbender.tools.pon.coverage;

import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;

/**
 * Coverage panel of normals interface.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public interface CoveragePanelOfNormals<T extends CoveragePoNNormalizationResult> {
    T normalize(final ReadCountCollection proportionalCoverageProfile, final JavaSparkContext ctx);

    default T normalize(final ReadCountCollection proportionalCoverageProfile) {
        return normalize(proportionalCoverageProfile, null);
    }

    T normalizeNormalsInPoN(final JavaSparkContext ctx);

    default T normalizeNormalsInPoN() {
        return normalizeNormalsInPoN(null);
    }
}