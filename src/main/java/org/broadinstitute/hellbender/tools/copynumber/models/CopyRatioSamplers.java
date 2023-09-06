package org.broadinstitute.hellbender.tools.copynumber.models;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.NaturalLogUtils;
import org.broadinstitute.hellbender.utils.mcmc.MinibatchSliceSampler;
import org.broadinstitute.hellbender.utils.mcmc.ParameterSampler;

import java.util.ArrayList;
import java.util.List;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class CopyRatioSamplers {
    private static final Logger logger = LogManager.getLogger(CopyRatioSamplers.class);

    private static final FunctionCache<Double> logCache = new FunctionCache<>(FastMath::log);

    private static final Function<Double, Double> UNIFORM_LOG_PRIOR = x -> 0.;
    private static final int GLOBAL_MINIBATCH_SIZE = 1000;
    private static final int SEGMENT_MINIBATCH_SIZE = 100;
    private static final double APPROX_THRESHOLD = 0.1;

    private CopyRatioSamplers() {}

    //Calculates the exponent for a normal distribution; used in log-likelihood calculation below.
    private static double normalTerm(final double quantity, 
                                     final double mean, 
                                     final double variance) {
        return (quantity - mean) * (quantity - mean) / (2. * variance);
    }

    //samples log conditional posterior for the variance parameter, assuming uniform prior; this is given by
    //the product of Gaussian likelihoods for each non-outlier point t:
    //  log[product_{non-outlier t} variance^(-1/2) * exp(-(log2cr_t - mean_t)^2 / (2 * variance))] + constant
    //where mean_t is identical for all points in a segment
    static final class VarianceSampler implements ParameterSampler<Double, CopyRatioParameter, CopyRatioState, CopyRatioSegmentedData> {
        private final double varianceMin;
        private final double varianceMax;
        private final double varianceSliceSamplingWidth;

        VarianceSampler(final double varianceMin, 
                        final double varianceMax, 
                        final double varianceSliceSamplingWidth) {
            this.varianceMin = varianceMin;
            this.varianceMax = varianceMax;
            this.varianceSliceSamplingWidth = varianceSliceSamplingWidth;
        }

        @Override
        public Double sample(final RandomGenerator rng, 
                             final CopyRatioState state, 
                             final CopyRatioSegmentedData data) {
            logger.debug("Sampling variance...");
            final List<CopyRatioSegmentedData.IndexedCopyRatio> nonOutlierIndexedCopyRatios =
                    data.getIndexedCopyRatios().stream()
                            .filter(icr -> !state.outlierIndicator(icr.getIndex()))
                            .collect(Collectors.toList());
            final BiFunction<CopyRatioSegmentedData.IndexedCopyRatio, Double, Double> logConditionalPDF = (icr, newVariance) ->
                    -0.5 * logCache.computeIfAbsent(newVariance)
                            - normalTerm(icr.getLog2CopyRatioValue(), state.segmentMean(icr.getSegmentIndex()), newVariance);
            return new MinibatchSliceSampler<>(
                    rng, nonOutlierIndexedCopyRatios, UNIFORM_LOG_PRIOR, logConditionalPDF,
                    varianceMin, varianceMax, varianceSliceSamplingWidth,
                    GLOBAL_MINIBATCH_SIZE, APPROX_THRESHOLD).sample(state.variance());
        }
    }

    //samples log conditional posterior for the outlier-probability parameter, assuming Beta(alpha, beta) prior;
    //this is given by:
    //  log Beta(alpha + number of outlier points, beta + number of non-outlier points) + constant
    static final class OutlierProbabilitySampler implements ParameterSampler<Double, CopyRatioParameter, CopyRatioState, CopyRatioSegmentedData> {
        private final double outlierProbabilityPriorAlpha;
        private final double outlierProbabilityPriorBeta;

        OutlierProbabilitySampler(final double outlierProbabilityPriorAlpha, 
                                  final double outlierProbabilityPriorBeta) {
            this.outlierProbabilityPriorAlpha = outlierProbabilityPriorAlpha;
            this.outlierProbabilityPriorBeta = outlierProbabilityPriorBeta;
        }

        @Override
        public Double sample(final RandomGenerator rng, 
                             final CopyRatioState state, 
                             final CopyRatioSegmentedData data) {
            logger.debug("Sampling outlier probability...");
            final int numOutliers = (int) IntStream.range(0, data.getNumPoints()).filter(state::outlierIndicator).count();
            return new BetaDistribution(rng,
                    outlierProbabilityPriorAlpha + numOutliers,
                    outlierProbabilityPriorBeta + data.getNumPoints() - numOutliers).sample();
        }
    }

    //samples log conditional posteriors for the segment-mean parameters, assuming uniform priors bounded by minimum and maximum log2 copy-ratio values;
    //for each segment s, this is given by the product of Gaussian likelihoods for each non-outlier point t:
    //  log[product_{non-outlier t in s} exp(-(log2cr_t - mean_s)^2 / (2 * variance))] + constant
    static final class SegmentMeansSampler implements ParameterSampler<CopyRatioState.SegmentMeans, CopyRatioParameter, CopyRatioState, CopyRatioSegmentedData> {
        private final double meanMin;
        private final double meanMax;
        private final double meanSliceSamplingWidth;

        SegmentMeansSampler(final double meanMin, 
                            final double meanMax, 
                            final double meanSliceSamplingWidth) {
            this.meanMin = meanMin;
            this.meanMax = meanMax;
            this.meanSliceSamplingWidth = meanSliceSamplingWidth;
        }

        @Override
        public CopyRatioState.SegmentMeans sample(final RandomGenerator rng,
                                                  final CopyRatioState state,
                                                  final CopyRatioSegmentedData data) {
            final List<Double> means = new ArrayList<>(data.getNumSegments());
            final BiFunction<CopyRatioSegmentedData.IndexedCopyRatio, Double, Double> logConditionalPDF = (icr, newMean) ->
                    state.outlierIndicator(icr.getIndex())
                            ? 0.
                            : -normalTerm(icr.getLog2CopyRatioValue(), newMean, state.variance());
            for (int segmentIndex = 0; segmentIndex < data.getNumSegments(); segmentIndex++) {
                final List<CopyRatioSegmentedData.IndexedCopyRatio> indexedCopyRatiosInSegment = data.getIndexedCopyRatiosInSegment(segmentIndex);
                if (indexedCopyRatiosInSegment.isEmpty()) {
                    means.add(Double.NaN);
                } else {
                    logger.debug(String.format("Sampling mean for segment %d...", segmentIndex));
                    final MinibatchSliceSampler<CopyRatioSegmentedData.IndexedCopyRatio> sampler = new MinibatchSliceSampler<>(
                            rng, indexedCopyRatiosInSegment, UNIFORM_LOG_PRIOR, logConditionalPDF,
                            meanMin, meanMax, meanSliceSamplingWidth,
                            SEGMENT_MINIBATCH_SIZE, APPROX_THRESHOLD);
                    means.add(sampler.sample(state.segmentMean(segmentIndex)));
                }
            }
            return new CopyRatioState.SegmentMeans(means);
        }
    }

    //samples log conditional posteriors for the outlier-indicator parameters; for each point t, this is given by:
    //          z_t * [log outlier_prob + outlierUniformLogLikelihood]
    //  + (1 - z_t) * [log((1 - outlier_prob) / (2 * pi * variance)^(1/2)) - (log2cr_t - mean_t)^2 / (2 * variance)]
    //  + const
    //where z_t is the indicator for point t, and outlier_prob is the outlier probability.
    //note that we compute the normalizing constant, so that we can sample a new indicator value by simply sampling
    //uniformly in [0, 1] and checking whether the resulting value is less than the probability of being an outlier
    //(corresponding to the first line in the unnormalized expression above)
    static final class OutlierIndicatorsSampler implements ParameterSampler<CopyRatioState.OutlierIndicators, CopyRatioParameter, CopyRatioState, CopyRatioSegmentedData> {
        private final double outlierUniformLogLikelihood;

        OutlierIndicatorsSampler(final double outlierUniformLogLikelihood) {
            this.outlierUniformLogLikelihood = outlierUniformLogLikelihood;
        }

        @Override
        public CopyRatioState.OutlierIndicators sample(final RandomGenerator rng,
                                                       final CopyRatioState state,
                                                       final CopyRatioSegmentedData data) {
            logger.debug("Sampling outlier indicators...");
            final double outlierUnnormalizedLogProbability =
                    Math.log(state.outlierProbability()) + outlierUniformLogLikelihood;
//            final double notOutlierUnnormalizedLogProbabilityPrefactor =
//                    Math.log(1. - state.outlierProbability()) - 0.5 * Math.log(2 * Math.PI * state.variance());
            final double notOutlierUnnormalizedLogProbabilityPrefactor =
                    Math.log((1. - state.outlierProbability()) / FastMath.sqrt(2 * Math.PI * state.variance()));
            final List<Boolean> indicators = new ArrayList<>(data.getNumPoints());
            for (int segmentIndex = 0; segmentIndex < data.getNumSegments(); segmentIndex++) {
                final List<CopyRatioSegmentedData.IndexedCopyRatio> indexedCopyRatiosInSegment = data.getIndexedCopyRatiosInSegment(segmentIndex);
                for (final CopyRatioSegmentedData.IndexedCopyRatio indexedCopyRatio : indexedCopyRatiosInSegment) {
                    final double notOutlierUnnormalizedLogProbability =
                            notOutlierUnnormalizedLogProbabilityPrefactor
                                    - normalTerm(indexedCopyRatio.getLog2CopyRatioValue(), state.segmentMean(segmentIndex), state.variance());
                    final double conditionalProbability =
                            FastMath.exp(outlierUnnormalizedLogProbability -
                                    NaturalLogUtils.logSumLog(outlierUnnormalizedLogProbability, notOutlierUnnormalizedLogProbability));
                    indicators.add(rng.nextDouble() < conditionalProbability);
                }
            }
            return new CopyRatioState.OutlierIndicators(indicators);
        }
    }
}
