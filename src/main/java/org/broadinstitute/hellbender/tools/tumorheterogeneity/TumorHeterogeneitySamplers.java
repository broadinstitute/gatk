package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import com.google.common.collect.Sets;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.special.Gamma;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.utils.mcmc.ParameterSampler;
import org.broadinstitute.hellbender.utils.mcmc.SliceSampler;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class TumorHeterogeneitySamplers {
    private static final double EPSILON = TumorHeterogeneityUtils.EPSILON;
    private static final double NUM_ENSEMBLE_SAMPLES_PER_HYPERPARAMETER_SAMPLE = 10;
    private static final double COPY_RATIO_NOISE_FLOOR_MAX = 1E-2;
    private static final double NOISE_FACTOR_MAX = 1E2;

    static final Logger logger = LogManager.getLogger(TumorHeterogeneitySamplers.class);

    private TumorHeterogeneitySamplers() {}

    static final class ConcentrationSampler implements ParameterSampler<Double, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private final double concentrationMin;
        private final double concentrationMax;
        private final double concentrationSliceSamplingWidth;

        private final int numWalkers;
        private int numSamples = 0;

        ConcentrationSampler(final double concentrationMin, final double concentrationMax, final double concentrationSliceSamplingWidth, final int numWalkers) {
            this.concentrationMin = concentrationMin;
            this.concentrationMax = concentrationMax;
            this.concentrationSliceSamplingWidth = concentrationSliceSamplingWidth;
            this.numWalkers = numWalkers;
        }

        @Override
        public Double sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
            if (numSamples % NUM_ENSEMBLE_SAMPLES_PER_HYPERPARAMETER_SAMPLE * numWalkers == 0) {
                final int numPopulations = state.populationMixture().numPopulations();
                final double concentrationPriorAlpha = state.priors().concentrationPriorHyperparameterValues().getAlpha();
                final double concentrationPriorBeta = state.priors().concentrationPriorHyperparameterValues().getBeta();
                final Function<Double, Double> logConditionalPDF = newConcentration -> {
                    final double populationFractionsTerm = IntStream.range(0, numPopulations)
                            .mapToDouble(i -> (newConcentration - 1) * Math.log(state.populationMixture().populationFraction(i) + EPSILON)).sum();
                    return (concentrationPriorAlpha - 1.) * Math.log(newConcentration) - concentrationPriorBeta * newConcentration +
                            Gamma.logGamma(newConcentration * numPopulations) - numPopulations * Gamma.logGamma(newConcentration) + populationFractionsTerm;
                };
                final double concentration = new SliceSampler(rng, logConditionalPDF, concentrationMin, concentrationMax, concentrationSliceSamplingWidth).sample(state.concentration());
                logger.debug("Sampled concentration: " + concentration);
                numSamples++;
                return concentration;
            } else {
                numSamples++;
                return state.concentration();
            }

        }
    }

    static final class CopyRatioNoiseFloorSampler implements ParameterSampler<Double, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private static final double MIN = 0.;
        private static final double MAX = COPY_RATIO_NOISE_FLOOR_MAX;
        private final double copyRatioNoiseFloorSliceSamplingWidth;

        private final int numWalkers;
        private int numSamples = 0;

        CopyRatioNoiseFloorSampler(final double copyRatioNoiseFloorSliceSamplingWidth, final int numWalkers) {
            this.copyRatioNoiseFloorSliceSamplingWidth = copyRatioNoiseFloorSliceSamplingWidth;
            this.numWalkers = numWalkers;
        }

        @Override
        public Double sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
            if (numSamples % NUM_ENSEMBLE_SAMPLES_PER_HYPERPARAMETER_SAMPLE *  numWalkers == 0) {
                final Function<Double, Double> logConditionalPDF = newCopyRatioNoiseFloor -> {
                    final TumorHeterogeneityState newState = new TumorHeterogeneityState(
                            state.concentration(),
                            newCopyRatioNoiseFloor,
                            state.copyRatioNoiseFactor(),
                            state.minorAlleleFractionNoiseFactor(),
                            state.populationMixture(),
                            state.priors());
                    return TumorHeterogeneityUtils.calculateLogPosterior(newState, data);
                };
                final double copyRatioNoiseFloor = new SliceSampler(rng, logConditionalPDF, MIN, MAX, copyRatioNoiseFloorSliceSamplingWidth).sample(state.copyRatioNoiseFloor());
                logger.debug("Sampled CR noise floor: " + copyRatioNoiseFloor);
                numSamples++;
                return copyRatioNoiseFloor;
            } else {
                numSamples++;
                return state.copyRatioNoiseFloor();
            }
//            return 0.;
        }
    }

    static final class CopyRatioNoiseFactorSampler implements ParameterSampler<Double, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private static final double MIN = 1.;
        private static final double MAX = NOISE_FACTOR_MAX;
        private final double copyRatioNoiseFactorSliceSamplingWidth;

        private final int numWalkers;
        private int numSamples = 0;

        CopyRatioNoiseFactorSampler(final double copyRatioNoiseFactorSliceSamplingWidth, final int numWalkers) {
            this.copyRatioNoiseFactorSliceSamplingWidth = copyRatioNoiseFactorSliceSamplingWidth;
            this.numWalkers = numWalkers;
        }

        @Override
        public Double sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
//            if (numSamples % NUM_ENSEMBLE_SAMPLES_PER_HYPERPARAMETER_SAMPLE *  numWalkers == 0) {
//                final Function<Double, Double> logConditionalPDF = newCopyRatioNoiseFactor -> {
//                    final TumorHeterogeneityState newState = new TumorHeterogeneityState(
//                            state.concentration(),
//                            state.copyRatioNoiseFloor(),
//                            newCopyRatioNoiseFactor,
//                            state.minorAlleleFractionNoiseFactor(),
//                            state.populationMixture(),
//                            state.priors());
//                    return TumorHeterogeneityUtils.calculateLogPosterior(newState, data);
//                };
//                final double copyRatioNoiseFactor = new SliceSampler(rng, logConditionalPDF, MIN, MAX, copyRatioNoiseFactorSliceSamplingWidth).sample(state.copyRatioNoiseFactor());
//                logger.debug("Sampled CR noise factor: " + copyRatioNoiseFactor);
//                numSamples++;
//                return copyRatioNoiseFactor;
//            } else {
//                numSamples++;
//                return state.copyRatioNoiseFactor();
//            }
            return 1.;
        }
    }

    static final class MinorAlleleFractionNoiseFactorSampler implements ParameterSampler<Double, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private static final double MIN = 1.;
        private static final double MAX = NOISE_FACTOR_MAX;
        private final double minorAlleleFractionNoiseFactorSliceSamplingWidth;

        private final int numWalkers;
        private int numSamples = 0;

        MinorAlleleFractionNoiseFactorSampler(final double minorAlleleFractionNoiseFactorSliceSamplingWidth, final int numWalkers) {
            this.minorAlleleFractionNoiseFactorSliceSamplingWidth = minorAlleleFractionNoiseFactorSliceSamplingWidth;
            this.numWalkers = numWalkers;
        }

        @Override
        public Double sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
//            if (numSamples % NUM_ENSEMBLE_SAMPLES_PER_HYPERPARAMETER_SAMPLE *  numWalkers == 0) {
//                final Function<Double, Double> logConditionalPDF = newMinorAlleleFractionNoiseFactor -> {
//                    final TumorHeterogeneityState newState = new TumorHeterogeneityState(
//                            state.concentration(),
//                            state.copyRatioNoiseFloor(),
//                            state.copyRatioNoiseFactor(),
//                            newMinorAlleleFractionNoiseFactor,
//                            state.populationMixture(),
//                            state.priors());
//                    return TumorHeterogeneityUtils.calculateLogPosterior(newState, data);
//                };
//                final double minorAlleleFractionNoiseFactor = new SliceSampler(rng, logConditionalPDF, MIN, MAX, minorAlleleFractionNoiseFactorSliceSamplingWidth).sample(state.minorAlleleFractionNoiseFactor());
//                logger.debug("Sampled MAF noise factor: " + minorAlleleFractionNoiseFactor);
//                numSamples++;
//                return minorAlleleFractionNoiseFactor;
//            } else {
//                numSamples++;
//                return state.minorAlleleFractionNoiseFactor();
//            }
            return 1.;
        }
    }

    static final class PopulationMixtureSampler implements ParameterSampler<PopulationMixture, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private static double maxLogPosterior = -1E100;
        private static final double scaleParameter = 2.;

        private final int numWalkers;
        private int numSamples = 0;
        private int numAccepted = 0;

        private final int maxTotalCopyNumber;
        private final List<List<Integer>> totalCopyNumberProductStates;
        private final Map<Integer, Set<PloidyState>> ploidyStateSetsMap = new HashMap<>();

        private final List<List<Double>> walkerPositions;

        PopulationMixtureSampler(final RandomGenerator rng, final int numPopulations, final List<PloidyState> ploidyStates, final PloidyState normalPloidyState, final int numWalkers) {
            final Set<Integer> totalCopyNumberStates = ploidyStates.stream().map(PloidyState::total).collect(Collectors.toSet());
            maxTotalCopyNumber = Collections.max(totalCopyNumberStates);
            totalCopyNumberProductStates = new ArrayList<>(Sets.cartesianProduct(Collections.nCopies(numPopulations, totalCopyNumberStates)));
            for (final int totalCopyNumber : totalCopyNumberStates) {
                final Set<PloidyState> ploidyStateSet = ploidyStates.stream().filter(ps -> ps.total() == totalCopyNumber).collect(Collectors.toSet());
                ploidyStateSetsMap.put(totalCopyNumber, ploidyStateSet);
            }

            this.numWalkers = numWalkers;
            walkerPositions = new ArrayList<>(numWalkers);
            for (int walkerIndex = 0; walkerIndex < numWalkers; walkerIndex++) {
                final List<Double> walkerPosition = IntStream.range(0, numPopulations - 1).boxed()
                        .map(i -> rng.nextGaussian()).collect(Collectors.toList());
                walkerPosition.add(Math.max(EPSILON, Math.min(normalPloidyState.total() + rng.nextGaussian(), maxTotalCopyNumber)));
                walkerPositions.add(walkerPosition);
            }
        }

        @Override
        public PopulationMixture sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
            final int numPopulations = state.populationMixture().numPopulations();
            final PloidyState normalPloidyState = state.priors().normalPloidyState();

            final int currentWalkerIndex = numSamples % numWalkers;
            final int selectedWalkerIndex = IntStream.range(0, numWalkers).boxed().filter(i -> i != currentWalkerIndex).collect(Collectors.toList()).get(rng.nextInt(numWalkers - 1));
            final List<Double> currentWalkerPosition = walkerPositions.get(currentWalkerIndex);
            final List<Double> selectedWalkerPosition = walkerPositions.get(selectedWalkerIndex);

            final double z = Math.pow((scaleParameter - 1.) * rng.nextDouble() + 1, 2.) / scaleParameter;
            final List<Double> proposedWalkerPosition = IntStream.range(0, numPopulations).boxed()
                    .map(i -> selectedWalkerPosition.get(i) + z * (currentWalkerPosition.get(i) - selectedWalkerPosition.get(i)))
                    .collect(Collectors.toList());
            proposedWalkerPosition.set(numPopulations - 1, Math.max(EPSILON, Math.min(proposedWalkerPosition.get(numPopulations - 1), maxTotalCopyNumber)));

            final List<Double> currentTransformedPopulationFractions = currentWalkerPosition.subList(0, numPopulations - 1);
            final PopulationMixture.PopulationFractions currentPopulationFractions =
                    TumorHeterogeneityUtils.calculatePopulationFractionsFromTransformedPopulationFractions(currentTransformedPopulationFractions);
            final double currentInitialPloidy = currentWalkerPosition.get(numPopulations - 1);
            final PopulationMixture.VariantProfileCollection currentVariantProfileCollection =
                    TumorHeterogeneityUtils.proposeVariantProfileCollection(rng, state, data,
                            currentPopulationFractions, currentInitialPloidy, totalCopyNumberProductStates, ploidyStateSetsMap);
            final PopulationMixture currentPopulationMixture = new PopulationMixture(
                    currentPopulationFractions, currentVariantProfileCollection, normalPloidyState);
            final double currentResultingPloidy = currentPopulationMixture.ploidy(data);
            logger.debug("Current population fractions: " + currentPopulationFractions);
            logger.debug("Current ploidy: " + currentResultingPloidy);

            final TumorHeterogeneityState currentState = new TumorHeterogeneityState(
                    state.concentration(),
                    state.copyRatioNoiseFloor(),
                    state.copyRatioNoiseFactor(),
                    state.minorAlleleFractionNoiseFactor(),
                    currentPopulationMixture,
                    state.priors());

            final List<Double> proposedTransformedPopulationFractions = proposedWalkerPosition.subList(0, numPopulations - 1);
            final PopulationMixture.PopulationFractions proposedPopulationFractions =
                    TumorHeterogeneityUtils.calculatePopulationFractionsFromTransformedPopulationFractions(proposedTransformedPopulationFractions);
            final double proposedInitialPloidy = proposedWalkerPosition.get(numPopulations - 1);
            final PopulationMixture.VariantProfileCollection proposedVariantProfileCollection =
                    TumorHeterogeneityUtils.proposeVariantProfileCollection(rng, state, data,
                            proposedPopulationFractions, proposedInitialPloidy, totalCopyNumberProductStates, ploidyStateSetsMap);
            final PopulationMixture proposedPopulationMixture = new PopulationMixture(
                    proposedPopulationFractions, proposedVariantProfileCollection, normalPloidyState);
            final double proposedResultingPloidy = proposedPopulationMixture.ploidy(data);
            logger.debug("Proposed population fractions: " + proposedPopulationFractions);
            logger.debug("Proposed initial ploidy: " + proposedInitialPloidy);
            logger.debug("Proposed resulting ploidy: " + proposedResultingPloidy);

            final TumorHeterogeneityState proposedState = new TumorHeterogeneityState(
                    state.concentration(),
                    state.copyRatioNoiseFloor(),
                    state.copyRatioNoiseFactor(),
                    state.minorAlleleFractionNoiseFactor(),
                    proposedPopulationMixture,
                    state.priors());

            final double proposedLogPosterior =
                    TumorHeterogeneityUtils.calculateLogPosterior(proposedState, data)
                    + TumorHeterogeneityUtils.calculateLogJacobianFactor(proposedPopulationFractions)
                    - 100. * Math.abs(proposedInitialPloidy - proposedResultingPloidy);
            final double currentLogPosterior =
                    TumorHeterogeneityUtils.calculateLogPosterior(currentState, data)
                    + TumorHeterogeneityUtils.calculateLogJacobianFactor(currentPopulationFractions)
                    - 100. * Math.abs(currentInitialPloidy - currentResultingPloidy);
            final double acceptanceProbability = Math.min(1., Math.exp((numPopulations - 1.) * Math.log(z) + proposedLogPosterior - currentLogPosterior));
            logger.debug("Log posterior of current state: " + currentLogPosterior);
            logger.debug("Log posterior of proposed state: " + proposedLogPosterior);
            numSamples += 1;
            if (proposedLogPosterior > maxLogPosterior) {
                maxLogPosterior = currentLogPosterior;
                logger.debug("New maximum: " + maxLogPosterior);
            }
            if (rng.nextDouble() < acceptanceProbability) {
                numAccepted += 1;
                logger.debug("Proposed state accepted.");
                System.out.println(proposedInitialPloidy + " " + proposedResultingPloidy + " " + Math.abs(proposedInitialPloidy - proposedResultingPloidy));
                walkerPositions.set(currentWalkerIndex, proposedWalkerPosition);
                return proposedPopulationMixture;
            }

            logger.debug("Acceptance rate: " + (double) numAccepted / numSamples);
            return currentPopulationMixture;
        }
    }
}
