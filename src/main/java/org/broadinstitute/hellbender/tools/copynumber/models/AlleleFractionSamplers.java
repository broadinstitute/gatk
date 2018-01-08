package org.broadinstitute.hellbender.tools.copynumber.models;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.mcmc.ParameterSampler;
import org.broadinstitute.hellbender.utils.mcmc.SliceSampler;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Sampler classes for the allele-fraction model.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class AlleleFractionSamplers {
    private static final Logger logger = LogManager.getLogger(AlleleFractionSamplers.class);

    private static final int NUM_POINTS_GLOBAL_SUBSAMPLE_THRESHOLD = 10000;
    private static final int NUM_POINTS_SEGMENT_SUBSAMPLE_THRESHOLD = 1000;

    private AlleleFractionSamplers() {}

    static final class MeanBiasSampler implements ParameterSampler<Double, AlleleFractionParameter, AlleleFractionState, AlleleFractionSegmentedData> {
        private static final double MIN_MEAN_BIAS = 0.;

        private final double maxMeanBias;
        private final double meanBiasSliceSamplingWidth;

        MeanBiasSampler(final double maxMeanBias,
                        final double meanBiasSliceSamplingWidth) {
            this.maxMeanBias = maxMeanBias;
            this.meanBiasSliceSamplingWidth = meanBiasSliceSamplingWidth;
        }

        @Override
        public Double sample(final RandomGenerator rng,
                             final AlleleFractionState state,
                             final AlleleFractionSegmentedData data) {
            logger.debug("Sampling mean bias...");
            final Function<AlleleFractionGlobalParameters, Double> logLikelihoodEstimate = logLikelihoodFromSubsample(
                    rng, state.minorFractions(), data, NUM_POINTS_GLOBAL_SUBSAMPLE_THRESHOLD);
            return new SliceSampler(rng,
                    x -> logLikelihoodEstimate.apply(state.globalParameters().copyWithNewMeanBias(x)),
                    MIN_MEAN_BIAS, maxMeanBias, meanBiasSliceSamplingWidth)
                    .sample(state.meanBias());
        }
    }

    static final class BiasVarianceSampler implements ParameterSampler<Double, AlleleFractionParameter, AlleleFractionState, AlleleFractionSegmentedData> {
        private static final double MIN_BIAS_VARIANCE = 1E-10;

        private final double maxBiasVariance;
        private final double biasVarianceSliceSamplingWidth;

        BiasVarianceSampler(final double maxBiasVariance,
                            final double biasVarianceSliceSamplingWidth) {
            this.maxBiasVariance = maxBiasVariance;
            this.biasVarianceSliceSamplingWidth = biasVarianceSliceSamplingWidth;
        }

        @Override
        public Double sample(final RandomGenerator rng,
                             final AlleleFractionState state,
                             final AlleleFractionSegmentedData data) {
            logger.debug("Sampling bias variance...");
            final Function<AlleleFractionGlobalParameters, Double> logLikelihoodEstimate = logLikelihoodFromSubsample(
                    rng, state.minorFractions(), data, NUM_POINTS_GLOBAL_SUBSAMPLE_THRESHOLD);
            return new SliceSampler(rng,
                    x -> logLikelihoodEstimate.apply(state.globalParameters().copyWithNewBiasVariance(x)),
                    MIN_BIAS_VARIANCE, maxBiasVariance, biasVarianceSliceSamplingWidth)
                    .sample(state.biasVariance());
        }
    }

    static final class OutlierProbabilitySampler implements ParameterSampler<Double, AlleleFractionParameter, AlleleFractionState, AlleleFractionSegmentedData> {
        private static final double MIN_OUTLIER_PROBABILITY = 0.;

        private final double maxOutlierProbability;
        private final double outlierProbabilitySliceSamplingWidth;

        OutlierProbabilitySampler(final double maxOutlierProbability,
                                  final double outlierProbabilitySliceSamplingWidth) {
            this.maxOutlierProbability = maxOutlierProbability;
            this.outlierProbabilitySliceSamplingWidth = outlierProbabilitySliceSamplingWidth;
        }

        @Override
        public Double sample(final RandomGenerator rng,
                             final AlleleFractionState state,
                             final AlleleFractionSegmentedData data) {
            logger.debug("Sampling outlier probability...");
            final Function<AlleleFractionGlobalParameters, Double> logLikelihoodEstimate = logLikelihoodFromSubsample(
                    rng, state.minorFractions(), data, NUM_POINTS_GLOBAL_SUBSAMPLE_THRESHOLD);
            return new SliceSampler(rng,
                    x -> logLikelihoodEstimate.apply(state.globalParameters().copyWithNewOutlierProbability(x)),
                    MIN_OUTLIER_PROBABILITY, maxOutlierProbability, outlierProbabilitySliceSamplingWidth)
                    .sample(state.outlierProbability());
        }
    }

    // sample minor fractions of all segments
    static final class MinorFractionsSampler implements ParameterSampler<AlleleFractionState.MinorFractions, AlleleFractionParameter, AlleleFractionState, AlleleFractionSegmentedData> {
        private static double MIN_MINOR_FRACTION = 0.;
        private static double MAX_MINOR_FRACTION = 0.5;
        private static final double PRIOR_BETA = 1.;

        private final Function<Double, Double> logPrior;
        private final List<Double> sliceSamplingWidths;

        MinorFractionsSampler(final AlleleFractionPrior prior,
                              final List<Double> sliceSamplingWidths) {
            logPrior = f -> new BetaDistribution(null, prior.getMinorAlleleFractionPriorAlpha(), PRIOR_BETA).logDensity(2 * f);
            this.sliceSamplingWidths = sliceSamplingWidths;
        }

        @Override
        public AlleleFractionState.MinorFractions sample(final RandomGenerator rng, final AlleleFractionState state, final AlleleFractionSegmentedData data) {
            final List<Double> minorFractions = new ArrayList<>(data.getNumSegments());
            for (int segment = 0; segment < data.getNumSegments(); segment++) {
                logger.debug(String.format("Sampling minor fraction for segment %d...", segment));
                final List<AlleleFractionSegmentedData.IndexedAllelicCount> allelicCountsInSegment =
                        data.getIndexedAllelicCountsInSegment(segment);
                if (allelicCountsInSegment.isEmpty()){
                    minorFractions.add(Double.NaN);
                } else {
                    final Function<Double, Double> segmentLogLikelihoodEstimate = segmentLogLikelihoodFromSubsample(
                            rng, state.globalParameters(), allelicCountsInSegment, NUM_POINTS_SEGMENT_SUBSAMPLE_THRESHOLD);
                    final SliceSampler sampler = new SliceSampler(rng,
                            f -> logPrior.apply(f) + segmentLogLikelihoodEstimate.apply(f),
                            MIN_MINOR_FRACTION, MAX_MINOR_FRACTION, sliceSamplingWidths.get(segment));
                    minorFractions.add(sampler.sample(state.segmentMinorFraction(segment)));
                }
            }
            return new AlleleFractionState.MinorFractions(minorFractions);
        }
    }

    private static List<AlleleFractionSegmentedData.IndexedAllelicCount> subsample(final RandomGenerator rng,
                                                                                   final List<AlleleFractionSegmentedData.IndexedAllelicCount> allelicCounts,
                                                                                   final int numPointsSubsampleThreshold) {
        //subsample the data if we are above the threshold
        return allelicCounts.size() > numPointsSubsampleThreshold
                ? IntStream.range(0, numPointsSubsampleThreshold).boxed().map(i -> rng.nextInt(allelicCounts.size())).map(allelicCounts::get).collect(Collectors.toList())
                : allelicCounts;
    }

    private static Function<AlleleFractionGlobalParameters, Double> logLikelihoodFromSubsample(final RandomGenerator rng,
                                                                                               final AlleleFractionState.MinorFractions minorFractions,
                                                                                               final AlleleFractionSegmentedData data,
                                                                                               final int numPointsSubsampleThreshold) {
        final List<AlleleFractionSegmentedData.IndexedAllelicCount> subsampledAllelicCounts =
                subsample(rng, data.getIndexedAllelicCounts(), numPointsSubsampleThreshold);
        final double scalingFactor = (double) data.getNumPoints() / subsampledAllelicCounts.size();
        final Map<Integer, List<AlleleFractionSegmentedData.IndexedAllelicCount>> segmentIndexToSubsampledAllelicCountsInSegmentMap =
                subsampledAllelicCounts.stream()
                        .collect(Collectors.groupingBy(AlleleFractionSegmentedData.IndexedAllelicCount::getSegmentIndex, Collectors.toList()));
        return parameters -> {
            double logLikelihood = 0.;
            for (final int segmentIndex : segmentIndexToSubsampledAllelicCountsInSegmentMap.keySet()) {
                logLikelihood += AlleleFractionLikelihoods.segmentLogLikelihood(
                        parameters, minorFractions.get(segmentIndex), segmentIndexToSubsampledAllelicCountsInSegmentMap.get(segmentIndex));
            }
            return scalingFactor * logLikelihood;
        };
    }

    private static Function<Double, Double> segmentLogLikelihoodFromSubsample(final RandomGenerator rng,
                                                                              final AlleleFractionGlobalParameters parameters,
                                                                              final List<AlleleFractionSegmentedData.IndexedAllelicCount> allelicCountsInSegment,
                                                                              final int numPointsSubsampleThreshold) {
        final List<AlleleFractionSegmentedData.IndexedAllelicCount> subsampledAllelicCountsInSegment =
                subsample(rng, allelicCountsInSegment, numPointsSubsampleThreshold);
        final double scalingFactor = (double) allelicCountsInSegment.size() / subsampledAllelicCountsInSegment.size();
        return minorFraction -> scalingFactor * AlleleFractionLikelihoods.segmentLogLikelihood(parameters, minorFraction, subsampledAllelicCountsInSegment);
    }
}
