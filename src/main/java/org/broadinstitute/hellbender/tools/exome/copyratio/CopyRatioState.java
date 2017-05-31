package org.broadinstitute.hellbender.tools.exome.copyratio;

import org.broadinstitute.hellbender.utils.mcmc.Parameter;
import org.broadinstitute.hellbender.utils.mcmc.ParameterizedState;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * The state of the copy-ratio model, containing: <p>
 *      1.  the global variance <p>
 *      2.  the global outlier probability <p>
 *      3.  means for each segment <p>
 *      4.  outlier indicators for each target <p>
 * <p>
 * See docs/CNVs/CNV-methods.pdf for details.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CopyRatioState extends ParameterizedState<CopyRatioParameter> {
    public static final class SegmentMeans extends ArrayList<Double> {
        private static final long serialVersionUID = 951753L;
        public SegmentMeans(final List<Double> segmentMeans) {
            super(new ArrayList<>(segmentMeans));
        }
    }

    public static final class OutlierIndicators extends ArrayList<Boolean> {
        private static final long serialVersionUID = 357159L;
        public OutlierIndicators(final List<Boolean> outlierIndicators) {
            super(new ArrayList<>(outlierIndicators));
        }
    }

    public CopyRatioState(final double variance, final double outlierProbability,
                          final SegmentMeans segmentMeans, final OutlierIndicators outlierIndicators) {
        super(Arrays.asList(
                new Parameter<>(CopyRatioParameter.VARIANCE, variance),
                new Parameter<>(CopyRatioParameter.OUTLIER_PROBABILITY, outlierProbability),
                new Parameter<>(CopyRatioParameter.SEGMENT_MEANS, segmentMeans),
                new Parameter<>(CopyRatioParameter.OUTLIER_INDICATORS, outlierIndicators)));
    }

    public double variance() {
        return get(CopyRatioParameter.VARIANCE, Double.class);
    }

    public double outlierProbability() {
        return get(CopyRatioParameter.OUTLIER_PROBABILITY, Double.class);
    }

    public double segmentMean(final int segment) {
        return get(CopyRatioParameter.SEGMENT_MEANS, CopyRatioState.SegmentMeans.class).get(segment);
    }

    public boolean targetOutlierIndicator(final int target) {
        return get(CopyRatioParameter.OUTLIER_INDICATORS, CopyRatioState.OutlierIndicators.class).get(target);
    }
}
