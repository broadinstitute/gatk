package org.broadinstitute.hellbender.tools.pon.coverage;

import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;

/**
 * Interface for a coverage profile that has been normalized by a {@link CoveragePanelOfNormals}.
 */
@FunctionalInterface
public interface CoveragePoNNormalizationResult {
    ReadCountCollection getProfile();
}
