package org.broadinstitute.hellbender.tools.exome.copyratio;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.mcmc.ParameterSampler;
import org.broadinstitute.hellbender.utils.mcmc.univariatesamplers.SliceSampler;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class CopyRatioSamplers {
    private CopyRatioSamplers() {}

    //Calculates the exponent for a normal distribution; used in log-likelihood calculation below.
    private static double normalTerm(final double quantity, final double mean, final double variance) {
        return (quantity - mean) * (quantity - mean) / (2. * variance);
    }

    //samples log conditional posterior for the variance parameter, assuming uniform prior; this is given by
    //the product of Gaussian likelihoods for each non-outlier target t:
    //  log[product_{non-outlier t} variance^(-1/2) * exp(-(coverage_t - mean_t)^2 / (2 * variance))] + constant
    //where mean_t is identical for all targets in a segment
    protected static final class VarianceSampler implements ParameterSampler<Double, CopyRatioParameter, CopyRatioState, CopyRatioData> {
        private final double varianceMin;
        private final double varianceMax;
        private final double varianceSliceSamplingWidth;

        public VarianceSampler(final double varianceMin, final double varianceMax, final double varianceSliceSamplingWidth) {
            this.varianceMin = varianceMin;
            this.varianceMax = varianceMax;
            this.varianceSliceSamplingWidth = varianceSliceSamplingWidth;
        }

        @Override
        public Double sample(final RandomGenerator rng, final CopyRatioState state, final CopyRatioData dataCollection) {
            final Function<Double, Double> logConditionalPDF = newVariance -> {
                final double gaussianLogNormalization = 0.5 * Math.log(newVariance);
                double ll = 0.;
                for (int segment = 0; segment < dataCollection.getNumSegments(); segment++) {
                    final List<CopyRatioData.IndexedCoverage> indexedCoveragesInSegment = dataCollection.getIndexedCoveragesInSegment(segment);
                    for (final CopyRatioData.IndexedCoverage c : indexedCoveragesInSegment) {
                        if (!state.targetOutlierIndicator(c.getTargetIndex())) {
                            ll -= normalTerm(c.getCoverage(), state.segmentMean(segment), newVariance) + gaussianLogNormalization;
                        }
                    }
                }
                return ll;
            };
            return new SliceSampler(rng, logConditionalPDF, varianceMin, varianceMax, varianceSliceSamplingWidth).sample(state.variance());
        }
    }


    //samples log conditional posterior for the outlier-probability parameter, assuming Beta(alpha, beta) prior;
    //this is given by:
    //  log Beta(alpha + number of outlier targets, beta + number of non-outlier targets) + constant
    protected static final class OutlierProbabilitySampler implements ParameterSampler<Double, CopyRatioParameter, CopyRatioState, CopyRatioData> {
        private final double outlierProbabilityPriorAlpha;
        private final double outlierProbabilityPriorBeta;

        public OutlierProbabilitySampler(final double outlierProbabilityPriorAlpha, final double outlierProbabilityPriorBeta) {
            this.outlierProbabilityPriorAlpha = outlierProbabilityPriorAlpha;
            this.outlierProbabilityPriorBeta = outlierProbabilityPriorBeta;
        }

        @Override
        public Double sample(final RandomGenerator rng, final CopyRatioState state, final CopyRatioData dataCollection) {
            final int numOutliers = (int) IntStream.range(0, dataCollection.getNumTargets()).filter(state::targetOutlierIndicator).count();
            return new BetaDistribution(rng,
                    outlierProbabilityPriorAlpha + numOutliers,
                    outlierProbabilityPriorBeta + dataCollection.getNumTargets() - numOutliers).sample();
        }
    }

    //samples log conditional posteriors for the segment-mean parameters, assuming uniform priors bounded by minimum and maximum coverage;
    //for each segment s, this is given by the product of Gaussian likelihoods for each non-outlier target t:
    //  log[product_{non-outlier t in s} exp(-(coverage_t - mean_s)^2 / (2 * variance))] + constant
    protected static final class SegmentMeansSampler implements ParameterSampler<CopyRatioState.SegmentMeans, CopyRatioParameter, CopyRatioState, CopyRatioData> {
        private final double coverageMin;
        private final double coverageMax;
        private final double meanSliceSamplingWidth;

        public SegmentMeansSampler(final double coverageMin, final double coverageMax, final double meanSliceSamplingWidth) {
            this.coverageMin = coverageMin;
            this.coverageMax = coverageMax;
            this.meanSliceSamplingWidth = meanSliceSamplingWidth;
        }

        @Override
        public CopyRatioState.SegmentMeans sample(final RandomGenerator rng, final CopyRatioState state, final CopyRatioData dataCollection) {
            final List<Double> means = new ArrayList<>(dataCollection.getNumSegments());
            for (int segment = 0; segment < dataCollection.getNumSegments(); segment++) {
                final List<CopyRatioData.IndexedCoverage> indexedCoveragesInSegment = dataCollection.getIndexedCoveragesInSegment(segment);
                if (indexedCoveragesInSegment.isEmpty()) {
                    means.add(Double.NaN);
                } else {
                    final Function<Double, Double> logConditionalPDF = newMean ->
                            indexedCoveragesInSegment.stream()
                                    .filter(c -> !state.targetOutlierIndicator(c.getTargetIndex()))
                                    .mapToDouble(c -> -normalTerm(c.getCoverage(), newMean, state.variance()))
                                    .sum();
                    //slice sample within range given by minimum and maximum coverages
                    final SliceSampler sampler = new SliceSampler(rng, logConditionalPDF, coverageMin, coverageMax, meanSliceSamplingWidth);
                    means.add(sampler.sample(state.segmentMean(segment)));
                }
            }
            return new CopyRatioState.SegmentMeans(means);
        }
    }

    //samples log conditional posteriors for the outlier-indicator parameters; for each target t, this is given by:
    //          z_t * [log outlier_prob + outlierUniformLogLikelihood]
    //  + (1 - z_t) * [log(1 - outlier_prob) - log(2 * pi * variance)/2 - (coverage_t - mean_t)^2 / (2 * variance)]
    //  + const
    //where z_t is the indicator for target t, and outlier_prob is the outlier probability.
    //note that we compute the normalizing constant, so that we can sample a new indicator value by simply sampling
    //uniformly in [0, 1] and checking whether the resulting value is less than the probability of being an outlier
    //(corresponding to the first line in the unnormalized expression above)
    protected static final class OutlierIndicatorsSampler implements ParameterSampler<CopyRatioState.OutlierIndicators, CopyRatioParameter, CopyRatioState, CopyRatioData> {
        private final double outlierUniformLogLikelihood;

        public OutlierIndicatorsSampler(final double outlierUniformLogLikelihood) {
            this.outlierUniformLogLikelihood = outlierUniformLogLikelihood;
        }

        @Override
        public CopyRatioState.OutlierIndicators sample(final RandomGenerator rng, final CopyRatioState state, final CopyRatioData dataCollection) {
            final double outlierUnnormalizedLogProbability =
                    Math.log(state.outlierProbability()) + outlierUniformLogLikelihood;
            final double notOutlierUnnormalizedLogProbabilityPrefactor =
                    Math.log(1. - state.outlierProbability()) - 0.5 * Math.log(2 * Math.PI * state.variance());
            final List<Boolean> indicators = new ArrayList<>();
            for (int segment = 0; segment < dataCollection.getNumSegments(); segment++) {
                final List<CopyRatioData.IndexedCoverage> indexedCoveragesInSegment = dataCollection.getIndexedCoveragesInSegment(segment);
                for (final CopyRatioData.IndexedCoverage c : indexedCoveragesInSegment) {
                    final double notOutlierUnnormalizedLogProbability =
                            notOutlierUnnormalizedLogProbabilityPrefactor
                                    - normalTerm(c.getCoverage(), state.segmentMean(segment), state.variance());
                    //note: we are working in natural log space, so we divide by ln(10) before using normalizeFromLog10
                    final double conditionalProbability =
                            MathUtils.normalizeFromLog10ToLinearSpace(new double[]{
                                    MathUtils.logToLog10(outlierUnnormalizedLogProbability),
                                    MathUtils.logToLog10(notOutlierUnnormalizedLogProbability)})[0];
                    indicators.add(rng.nextDouble() < conditionalProbability);
                }
            }
            return new CopyRatioState.OutlierIndicators(indicators);
        }
    }
}
