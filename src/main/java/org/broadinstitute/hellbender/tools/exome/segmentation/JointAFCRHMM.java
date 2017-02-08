package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.List;

/**
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class JointAFCRHMM extends ClusteringGenomicHMM<JointSegmentationDatum, AFCRHiddenState> {
    private final double logCoverageCauchyWidth;
    private final double log10OutlierProbability;
    private final double log10NonOutlierProbability;

    public JointAFCRHMM(final List<AFCRHiddenState> hiddenStateValues, final double memoryLength, final double logCoverageCauchyWidth, final double outlierProbability) {
        super(hiddenStateValues, memoryLength);
        this.logCoverageCauchyWidth = logCoverageCauchyWidth;
        log10OutlierProbability = Math.log10(outlierProbability);
        log10NonOutlierProbability = Math.log10(1 - outlierProbability);
    }

    @Override
    public double logEmissionProbability(final JointSegmentationDatum datum, final Integer state, final SimpleInterval position) {
        return logEmissionProbability(datum, getHiddenStateValue(state));
    }

    @Override
    public double logEmissionProbability(final JointSegmentationDatum datum, final AFCRHiddenState hiddenState) {
        return datum.isTarget() ?
                CopyRatioHMM.logEmissionProbability(datum.getCopyRatio(), hiddenState.getLog2CopyRatio(), logCoverageCauchyWidth)
                : AlleleFractionHMM.logEmissionProbability(datum.getAllelicCount(), hiddenState.getMinorAlleleFraction(),
                log10OutlierProbability, log10NonOutlierProbability);
    }
}
