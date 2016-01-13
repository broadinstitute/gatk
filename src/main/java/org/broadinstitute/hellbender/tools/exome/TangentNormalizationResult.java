package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.math3.linear.RealMatrix;

/**
 * Stores the results of a tangent normalization.
 */
public class TangentNormalizationResult {
    ReadCountCollection tangentNormalized;
    ReadCountCollection preTangentNormalized;
    RealMatrix tangentBetaHats;
    ReadCountCollection targetFactorNormalizedCounts;

    public TangentNormalizationResult(final ReadCountCollection tangentNormalized, final ReadCountCollection preTangentNormalized, final RealMatrix tangentBetaHats, final ReadCountCollection targetFactorNormalizedCounts) {
        this.tangentNormalized = tangentNormalized;
        this.preTangentNormalized = preTangentNormalized;
        this.tangentBetaHats = tangentBetaHats;
        this.targetFactorNormalizedCounts = targetFactorNormalizedCounts;
    }

    public ReadCountCollection getTangentNormalized() {
        return tangentNormalized;
    }

    public ReadCountCollection getPreTangentNormalized() {
        return preTangentNormalized;
    }

    public RealMatrix getTangentBetaHats() {
        return tangentBetaHats;
    }

    public ReadCountCollection getTargetFactorNormalizedCounts() {
        return targetFactorNormalizedCounts;
    }
}
