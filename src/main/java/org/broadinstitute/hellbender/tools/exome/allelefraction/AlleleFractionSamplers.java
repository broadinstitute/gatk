package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.apache.commons.math3.random.RandomGenerator;
import org.broadinstitute.hellbender.utils.mcmc.AdaptiveMetropolisSampler;
import org.broadinstitute.hellbender.utils.mcmc.Sampler;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Sampler classes for the allele-fraction model.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AlleleFractionSamplers {
    private AlleleFractionSamplers() {}

    // sample mean bias
    protected static final class MeanBiasSampler implements Sampler<Double, AlleleFractionState, AlleleFractionData> {
        private final AdaptiveMetropolisSampler sampler;

        public MeanBiasSampler(final AlleleFractionState initialState, final double initialStepSize) {
            sampler = new AdaptiveMetropolisSampler(initialState.meanBias(), initialStepSize, 0, Double.POSITIVE_INFINITY);
        }

        public Double sample(final RandomGenerator rng, final AlleleFractionState state, final AlleleFractionData data) {
            final AllelicPanelOfNormals allelicPON = data.getPON();
            if (allelicPON.equals(AllelicPanelOfNormals.EMPTY_PON)) {
                return sampler.sample(x -> AlleleFractionLikelihoods.logLikelihood(state.shallowCopyWithProposedMeanBias(x), data));
            }
            return allelicPON.getMLEMeanBias(); // if PON is available, always return MLE mean bias as "sample"
        }
    }

    // sample bias variance
    protected static final class BiasVarianceSampler implements Sampler<Double, AlleleFractionState, AlleleFractionData> {
        private final AdaptiveMetropolisSampler sampler;

        public BiasVarianceSampler(final AlleleFractionState initialState, final double initialStepSize) {
            sampler = new AdaptiveMetropolisSampler(initialState.biasVariance(), initialStepSize, 0, Double.POSITIVE_INFINITY);
        }

        public Double sample(final RandomGenerator rng, final AlleleFractionState state, final AlleleFractionData data) {
            final AllelicPanelOfNormals allelicPON = data.getPON();
            if (allelicPON.equals(AllelicPanelOfNormals.EMPTY_PON)) {
                return sampler.sample(x -> AlleleFractionLikelihoods.logLikelihood(state.shallowCopyWithProposedBiasVariance(x), data));
            }
            return allelicPON.getMLEBiasVariance(); // if PON is available, always return MLE bias variance as "sample"
        }
    }

    // sample outlier probability
    protected static final class OutlierProbabilitySampler implements Sampler<Double, AlleleFractionState, AlleleFractionData> {
        private final AdaptiveMetropolisSampler sampler;

        public OutlierProbabilitySampler(final AlleleFractionState initialState, final double initialStepSize) {
            sampler = new AdaptiveMetropolisSampler(initialState.outlierProbability(), initialStepSize, 0, Double.POSITIVE_INFINITY);
        }

        public Double sample(final RandomGenerator rng, final AlleleFractionState state, final AlleleFractionData data) {
            return sampler.sample(x -> AlleleFractionLikelihoods.logLikelihood(state.shallowCopyWithProposedOutlierProbability(x), data));
        }
    }

    // sample minor fraction of a single segment
    private static final class PerSegmentMinorFractionSampler implements Sampler<Double, AlleleFractionState, AlleleFractionData> {
        private static final double MAX_MINOR_FRACTION = 0.5;   //by definition!
        private final AdaptiveMetropolisSampler sampler;

        private final int segmentIndex;

        public PerSegmentMinorFractionSampler(final int segmentIndex, final AlleleFractionState initialState, final double initialStepSize) {
            sampler = new AdaptiveMetropolisSampler(initialState.minorFractionInSegment(segmentIndex), initialStepSize, 0, MAX_MINOR_FRACTION);
            this.segmentIndex = segmentIndex;
        }

        public Double sample(final RandomGenerator rng, final AlleleFractionState state, final AlleleFractionData data) {
            if (data.numHetsInSegment(segmentIndex) == 0) {
                return Double.NaN;
            }
            return sampler.sample(AlleleFractionLikelihoods.segmentLogLikelihoodConditionalOnMinorFraction(state, data, segmentIndex));
        }
    }

    // sample minor fractions of all segments
    protected static final class MinorFractionsSampler implements Sampler<AlleleFractionState.MinorFractions, AlleleFractionState, AlleleFractionData> {
        private final List<PerSegmentMinorFractionSampler> perSegmentSamplers = new ArrayList<>();
        private final int numSegments;

        public MinorFractionsSampler(final AlleleFractionState initialState,
                                     final List<Double> initialStepSizes) {
            this.numSegments = initialState.minorFractions().size();
            for (int segment = 0; segment < numSegments; segment++) {
                perSegmentSamplers.add(new PerSegmentMinorFractionSampler(segment, initialState, initialStepSizes.get(segment)));
            }
        }

        public AlleleFractionState.MinorFractions sample(final RandomGenerator rng, final AlleleFractionState state, final AlleleFractionData data) {
            return new AlleleFractionState.MinorFractions(perSegmentSamplers.stream()
                    .map(sampler -> sampler.sample(rng, state, data)).collect(Collectors.toList()));
        }
    }
}
