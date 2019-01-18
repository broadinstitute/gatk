package org.broadinstitute.hellbender.tools.copynumber.models;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.mcmc.MinibatchSliceSampler;
import org.broadinstitute.hellbender.utils.mcmc.ParameterSampler;

import java.util.ArrayList;
import java.util.List;
import java.util.function.BiFunction;
import java.util.function.Function;

/**
 * Sampler classes for the allele-fraction model.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class AlleleFractionSamplers {
    private static final Logger logger = LogManager.getLogger(AlleleFractionSamplers.class);

    private static final Function<Double, Double> UNIFORM_LOG_PRIOR = x -> 0.;
    private static final int GLOBAL_MINIBATCH_SIZE = 1000;
    private static final int SEGMENT_MINIBATCH_SIZE = 10;
    private static final double APPROX_THRESHOLD = 0.1;

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
            final BiFunction<AlleleFractionSegmentedData.IndexedAllelicCount, Double, Double> logConditionalPDF = (iac, newMeanBias) ->
                    AlleleFractionLikelihoods.hetLogLikelihood(
                            state.globalParameters().copyWithNewMeanBias(newMeanBias),
                            state.segmentMinorFraction(iac.getSegmentIndex()),
                            iac);
            return new MinibatchSliceSampler<>(
                    rng, data.getIndexedAllelicCounts(), UNIFORM_LOG_PRIOR, logConditionalPDF,
                    MIN_MEAN_BIAS, maxMeanBias, meanBiasSliceSamplingWidth,
                    GLOBAL_MINIBATCH_SIZE, APPROX_THRESHOLD).sample(state.globalParameters().getMeanBias());
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
            final BiFunction<AlleleFractionSegmentedData.IndexedAllelicCount, Double, Double> logConditionalPDF = (iac, newBiasVariance) ->
                    AlleleFractionLikelihoods.hetLogLikelihood(
                            state.globalParameters().copyWithNewBiasVariance(newBiasVariance),
                            state.segmentMinorFraction(iac.getSegmentIndex()),
                            iac);
            return new MinibatchSliceSampler<>(
                    rng, data.getIndexedAllelicCounts(), UNIFORM_LOG_PRIOR, logConditionalPDF,
                    MIN_BIAS_VARIANCE, maxBiasVariance, biasVarianceSliceSamplingWidth,
                    GLOBAL_MINIBATCH_SIZE, APPROX_THRESHOLD).sample(state.globalParameters().getBiasVariance());
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
            final BiFunction<AlleleFractionSegmentedData.IndexedAllelicCount, Double, Double> logConditionalPDF = (iac, newOutlierProbability) ->
                    AlleleFractionLikelihoods.hetLogLikelihood(
                            state.globalParameters().copyWithNewOutlierProbability(newOutlierProbability),
                            state.segmentMinorFraction(iac.getSegmentIndex()),
                            iac);
            return new MinibatchSliceSampler<>(
                    rng, data.getIndexedAllelicCounts(), UNIFORM_LOG_PRIOR, logConditionalPDF,
                    MIN_OUTLIER_PROBABILITY, maxOutlierProbability, outlierProbabilitySliceSamplingWidth,
                    GLOBAL_MINIBATCH_SIZE, APPROX_THRESHOLD).sample(state.globalParameters().getOutlierProbability());
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
            final BiFunction<AlleleFractionSegmentedData.IndexedAllelicCount, Double, Double> logConditionalPDF = (iac, newMinorFraction) ->
                    AlleleFractionLikelihoods.hetLogLikelihood(state.globalParameters(), newMinorFraction, iac);
            for (int segmentIndex = 0; segmentIndex < data.getNumSegments(); segmentIndex++) {
                logger.debug(String.format("Sampling minor fraction for segment %d...", segmentIndex));
                final List<AlleleFractionSegmentedData.IndexedAllelicCount> allelicCountsInSegment =
                        data.getIndexedAllelicCountsInSegment(segmentIndex);
                if (allelicCountsInSegment.isEmpty()){
                    minorFractions.add(Double.NaN);
                } else {
                    final MinibatchSliceSampler<AlleleFractionSegmentedData.IndexedAllelicCount> sampler =
                            new MinibatchSliceSampler<>(
                                    rng, allelicCountsInSegment, logPrior, logConditionalPDF,
                                    MIN_MINOR_FRACTION, MAX_MINOR_FRACTION, sliceSamplingWidths.get(segmentIndex),
                                    SEGMENT_MINIBATCH_SIZE, APPROX_THRESHOLD);
                    minorFractions.add(sampler.sample(state.segmentMinorFraction(segmentIndex)));
                }
            }
            return new AlleleFractionState.MinorFractions(minorFractions);
        }
    }
}
