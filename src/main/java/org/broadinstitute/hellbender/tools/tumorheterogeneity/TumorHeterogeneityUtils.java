package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import com.google.common.collect.Sets;
import org.apache.commons.collections4.ListUtils;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.PopulationMixture.VariantProfileCollection;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyStatePrior;
import org.broadinstitute.hellbender.utils.mcmc.coordinates.CoordinateUtils;
import org.broadinstitute.hellbender.utils.mcmc.coordinates.SimplexPosition;
import org.broadinstitute.hellbender.utils.mcmc.coordinates.WalkerPosition;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Helper class containing package-private constants and methods used in org.broadinstitute.hellbender.tools.tumorheterogeneity.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class TumorHeterogeneityUtils {
    private static final Logger logger = LogManager.getLogger(TumorHeterogeneityUtils.class);

    static final int NUM_POPULATIONS_CLONAL = 2;

    static final double INITIAL_WALKER_BALL_SIZE_CLONAL = 1.;
    static final double INITIAL_WALKER_BALL_SIZE = 1.;

    static final double EPSILON = 1E-10;
    static final double COPY_RATIO_EPSILON = 1E-2; //below this, use mirrored minor-allele fraction posterior

    //fixes concentration to practically unity for clonal-only version
    static final double CONCENTRATION_PRIOR_ALPHA_CLONAL = 1E6;
    static final double CONCENTRATION_PRIOR_BETA_CLONAL = 1E6;

    static final double CONCENTRATION_MIN = EPSILON;
    static final double CONCENTRATION_MAX = 1. + EPSILON;

    static final double COPY_RATIO_NORMALIZATION_MIN = 0.5;
    static final double COPY_RATIO_NORMALIZATION_MAX = 2.;

    static final double COPY_RATIO_NOISE_CONSTANT_MIN = EPSILON;
    static final double COPY_RATIO_NOISE_CONSTANT_MAX = 1E-1;

    static final double PLOIDY_MIN = 1E-1;

    private static final int CONCENTRATION_WALKER_DIMENSION_INDEX = 0;
    private static final int COPY_RATIO_NORMALIZATION_WALKER_DIMENSION_INDEX = 1;
    private static final int COPY_RATIO_NOISE_CONSTANT_WALKER_DIMENSION_INDEX = 2;
    private static final int INITIAL_PLOIDY_WALKER_DIMENSION_INDEX = 3;
    private static final int POPULATION_FRACTIONS_WALKER_DIMENSION_START_INDEX = 4;
    static final int NUM_GLOBAL_PARAMETERS = 4;

    private TumorHeterogeneityUtils() {}

    static double calculateLogPosterior(final TumorHeterogeneityState state,
                                        final TumorHeterogeneityData data) {
        if (isOutsideBounds(state, data)) {
            return Double.NEGATIVE_INFINITY;
        }

        //concentration prior
        final double concentrationPriorAlpha = data.priors().concentrationPriorHyperparameterValues().getAlpha();
        final double concentrationPriorBeta = data.priors().concentrationPriorHyperparameterValues().getBeta();
        final double concentration = state.concentration();
        final double logPriorConcentration =
                concentrationPriorAlpha * FastMath.log(Math.max(EPSILON, concentrationPriorBeta))
                        + (concentrationPriorAlpha - 1.) * FastMath.log(Math.max(EPSILON, concentration))
                        - concentrationPriorBeta * concentration
                        - Gamma.logGamma(concentrationPriorAlpha);

        //copy-ratio normalization prior
        final double copyRatioNormalizationPriorAlpha = data.priors().copyRatioNormalizationPriorHyperparameterValues().getAlpha();
        final double copyRatioNormalizationPriorBeta = data.priors().copyRatioNormalizationPriorHyperparameterValues().getBeta();
        final double copyRatioNormalization = state.copyRatioNormalization();
        final double logPriorCopyRatioNormalization =
                copyRatioNormalizationPriorAlpha * FastMath.log(Math.max(EPSILON, copyRatioNormalizationPriorBeta))
                        + (copyRatioNormalizationPriorAlpha - 1.) * FastMath.log(copyRatioNormalization)
                        - copyRatioNormalizationPriorBeta * copyRatioNormalization
                        - Gamma.logGamma(copyRatioNormalizationPriorAlpha);

        //copy-ratio noise-constant prior
        final double copyRatioNoiseConstantPriorAlpha = data.priors().copyRatioNoiseConstantPriorHyperparameterValues().getAlpha();
        final double copyRatioNoiseConstantPriorBeta = data.priors().copyRatioNoiseConstantPriorHyperparameterValues().getBeta();
        final double copyRatioNoiseConstant = state.copyRatioNoiseConstant();
        final double logPriorCopyRatioNoiseConstant =
                copyRatioNoiseConstantPriorAlpha * FastMath.log(Math.max(EPSILON, copyRatioNoiseConstantPriorBeta))
                        + (copyRatioNoiseConstantPriorAlpha - 1.) * FastMath.log(copyRatioNoiseConstant)
                        - copyRatioNoiseConstantPriorBeta * copyRatioNoiseConstant
                        - Gamma.logGamma(copyRatioNoiseConstantPriorAlpha);

        //population-fractions prior
        final int numPopulations = state.populationMixture().numPopulations();
        final double logPriorPopulationFractionsSum = IntStream.range(0, numPopulations)
                .mapToDouble(i -> (concentration - 1.) * FastMath.log(Math.max(EPSILON, state.populationMixture().populationFraction(i))))
                .sum();
        final double logPriorPopulationFractions =
                Gamma.logGamma(concentration * numPopulations)
                        - numPopulations * Gamma.logGamma(concentration)
                        + logPriorPopulationFractionsSum;

        //variant-profiles prior
        double logPriorVariantProfiles = 0.;
        final VariantProfileCollection variantProfileCollection = state.populationMixture().variantProfileCollection();
        final int numVariantPopulations = variantProfileCollection.numVariantPopulations();
        for (int segmentIndex = 0; segmentIndex < variantProfileCollection.numSegments(); segmentIndex++) {
            final Set<PloidyState> ploidyStatesInSegment = new HashSet<>();
            for (int populationIndex = 0; populationIndex < numVariantPopulations; populationIndex++) {
                final PloidyState ploidyState = variantProfileCollection.ploidyState(populationIndex, segmentIndex);
                logPriorVariantProfiles += data.priors().ploidyStatePrior().logProbability(ploidyState);
                ploidyStatesInSegment.add(ploidyState);
            }
            logPriorVariantProfiles += -data.priors().subcloneVariancePenalty() * ploidyStatesInSegment.size();
        }

        //copy-ratio--minor-allele-fraction likelihood
        double logLikelihoodSegments = 0.;
        final double ploidy = state.ploidy();
        for (int segmentIndex = 0; segmentIndex < data.numSegments(); segmentIndex++) {
            final double mAlleleCopyNumber = state.populationMixture().calculatePopulationAveragedCopyNumberFunction(segmentIndex, PloidyState::m);
            final double nAlleleCopyNumber = state.populationMixture().calculatePopulationAveragedCopyNumberFunction(segmentIndex, PloidyState::n);
            final double totalCopyNumber = mAlleleCopyNumber + nAlleleCopyNumber;
            final double copyRatio = state.copyRatioNormalization() * totalCopyNumber / Math.max(EPSILON, ploidy);
            final double minorAlleleFraction = calculateMinorAlleleFraction(mAlleleCopyNumber, nAlleleCopyNumber);
            logLikelihoodSegments += data.logDensity(segmentIndex, copyRatio, minorAlleleFraction, state.copyRatioNoiseConstant());
        }

        //log penalty for mismatch between initial ploidy and ploidy resulting from proposeVariantProfileCollection
        final double logPloidyMismatchPenalty = -data.priors().ploidyMismatchPenalty() * Math.abs(state.initialPloidy() - state.ploidy());

        logger.debug("Log-posterior components:"
                + " " + logPriorConcentration
                + " " + logPriorCopyRatioNormalization
                + " " + logPriorCopyRatioNoiseConstant
                + " " + logPriorPopulationFractions
                + " " + logPriorVariantProfiles
                + " " + logLikelihoodSegments
                + " " + logPloidyMismatchPenalty);

        return logPriorConcentration + logPriorCopyRatioNormalization + logPriorCopyRatioNoiseConstant
                + logPriorPopulationFractions + logPriorVariantProfiles + logLikelihoodSegments + logPloidyMismatchPenalty;
    }

    /**
     * Calculates the log of the Jacobian factor for the state-to-walker transformation.
     */
    static double calculateLogJacobianFactor(final TumorHeterogeneityState state,
                                             final TumorHeterogeneityData data) {
        final double ploidyMax = data.priors().ploidyStatePrior().maxCopyNumber();
        return CoordinateUtils.calculateLogJacobianFactor(state.concentration(), CONCENTRATION_MIN, CONCENTRATION_MAX)
                + CoordinateUtils.calculateLogJacobianFactor(state.copyRatioNormalization(), COPY_RATIO_NORMALIZATION_MIN, COPY_RATIO_NORMALIZATION_MAX)
                + CoordinateUtils.calculateLogJacobianFactor(state.copyRatioNoiseConstant(), COPY_RATIO_NOISE_CONSTANT_MIN, COPY_RATIO_NOISE_CONSTANT_MAX)
                + CoordinateUtils.calculateLogJacobianFactor(state.populationMixture().ploidy(data), PLOIDY_MIN, ploidyMax)
                + SimplexPosition.calculateLogJacobianFactor(state.populationMixture().populationFractions());
    }

    /**
     * Transforms a walker position (i.e., a point in unbounded N-dimensional space) to
     * a {@link TumorHeterogeneityState} (which is composed of parameters that may be bounded or confined to the
     * unit simplex).  This transformation includes a proposal step to find the maximum a posteriori
     * {@link VariantProfileCollection} conditional on the other model parameters and a proposed "initial" ploidy
     * using {@link TumorHeterogeneityUtils#proposeVariantProfileCollection}.
     */
    static TumorHeterogeneityState transformWalkerPositionToState(final WalkerPosition walkerPosition,
                                                                  final TumorHeterogeneityData data,
                                                                  final List<List<Integer>> copyNumberProductStates,
                                                                  final Map<Integer, Set<PloidyState>> ploidyStateSetsMap) {
        final double concentration = CoordinateUtils.transformWalkerCoordinateToBoundedVariable(
                walkerPosition.get(CONCENTRATION_WALKER_DIMENSION_INDEX), CONCENTRATION_MIN, CONCENTRATION_MAX);
        final double copyRatioNormalization = CoordinateUtils.transformWalkerCoordinateToBoundedVariable(
                walkerPosition.get(COPY_RATIO_NORMALIZATION_WALKER_DIMENSION_INDEX), COPY_RATIO_NORMALIZATION_MIN, COPY_RATIO_NORMALIZATION_MAX);
        final double copyRatioNoiseConstant = CoordinateUtils.transformWalkerCoordinateToBoundedVariable(
                walkerPosition.get(COPY_RATIO_NOISE_CONSTANT_WALKER_DIMENSION_INDEX), COPY_RATIO_NOISE_CONSTANT_MIN, COPY_RATIO_NOISE_CONSTANT_MAX);
        final double ploidyMax = data.priors().ploidyStatePrior().maxCopyNumber();
        final double initialPloidy = CoordinateUtils.transformWalkerCoordinateToBoundedVariable(
                walkerPosition.get(INITIAL_PLOIDY_WALKER_DIMENSION_INDEX), PLOIDY_MIN, ploidyMax);

        final PopulationMixture.PopulationFractions populationFractions = new PopulationMixture.PopulationFractions(
                SimplexPosition.calculateSimplexPositionFromWalkerPosition(
                        new WalkerPosition(walkerPosition.subList(POPULATION_FRACTIONS_WALKER_DIMENSION_START_INDEX, walkerPosition.numDimensions()))));

        final VariantProfileCollection variantProfileCollection = proposeVariantProfileCollection(
                copyRatioNormalization, copyRatioNoiseConstant,  initialPloidy,
                populationFractions, data, copyNumberProductStates, ploidyStateSetsMap);

        final PloidyState normalPloidyState = data.priors().normalPloidyState();
        final PopulationMixture populationMixture = new PopulationMixture(populationFractions, variantProfileCollection, normalPloidyState);
        final double ploidy = populationMixture.ploidy(data);

        return new TumorHeterogeneityState(concentration, copyRatioNormalization, copyRatioNoiseConstant, initialPloidy, ploidy, populationMixture);
    }

    /**
     * Transforms a {@link TumorHeterogeneityState} to a walker position.  This is only used to initialize a ball
     * of walkers around an initial state in {@link TumorHeterogeneityModeller#initializeWalkerBall}.
     */
    static WalkerPosition transformStateToWalkerPosition(final TumorHeterogeneityState state,
                                                         final TumorHeterogeneityData data) {
        final double concentrationWalkerCoordinate = CoordinateUtils.transformBoundedVariableToWalkerCoordinate(
                state.concentration(), CONCENTRATION_MIN, CONCENTRATION_MAX);
        final double copyRatioNormalizationWalkerCoordinate = CoordinateUtils.transformBoundedVariableToWalkerCoordinate(
                state.copyRatioNormalization(), COPY_RATIO_NORMALIZATION_MIN, COPY_RATIO_NORMALIZATION_MAX);
        final double copyRatioNoiseConstantWalkerCoordinate = CoordinateUtils.transformBoundedVariableToWalkerCoordinate(
                state.copyRatioNoiseConstant(), COPY_RATIO_NOISE_CONSTANT_MIN, COPY_RATIO_NOISE_CONSTANT_MAX);
        final double ploidyMax = data.priors().ploidyStatePrior().maxCopyNumber();
        final double initialPloidyWalkerCoordinate = CoordinateUtils.transformBoundedVariableToWalkerCoordinate(
                state.initialPloidy(), PLOIDY_MIN, ploidyMax);
        final WalkerPosition populationFractionsWalkerCoordinates = SimplexPosition.calculateWalkerPositionFromSimplexPosition(state.populationMixture().populationFractions());

        return new WalkerPosition(ListUtils.union(
                Arrays.asList(
                        concentrationWalkerCoordinate,
                        copyRatioNormalizationWalkerCoordinate,
                        copyRatioNoiseConstantWalkerCoordinate,
                        initialPloidyWalkerCoordinate),
                populationFractionsWalkerCoordinates));
    }

    /**
     * Conditional on the global model parameters, population fractions, and a proposed "initial" ploidy, heuristically
     * samples a {@link VariantProfileCollection} from the posterior distribution.
     */
    private static VariantProfileCollection proposeVariantProfileCollection(final double copyRatioNormalization,
                                                                            final double copyRatioNoiseConstant,
                                                                            final double initialPloidy,
                                                                            final PopulationMixture.PopulationFractions populationFractions,
                                                                            final TumorHeterogeneityData data,
                                                                            final List<List<Integer>> copyNumberProductStates,
                                                                            final Map<Integer, Set<PloidyState>> ploidyStateSetsMap) {
        final int numVariantPopulations = populationFractions.numVariantPopulations();
        final int numSegments = data.numSegments();
        final PloidyState normalPloidyState = data.priors().normalPloidyState();
        final PloidyStatePrior ploidyStatePrior = data.priors().ploidyStatePrior();

        //initialize a list of variant profiles to store the result
        final List<PopulationMixture.VariantProfile> variantProfiles =
                TumorHeterogeneityState.initializeNormalProfiles(numVariantPopulations, numSegments, normalPloidyState);

        for (int segmentIndex = 0; segmentIndex < numSegments; segmentIndex++) {
            final int si = segmentIndex;

            //for all possible copy-number product states, calculate the copy-ratio posterior given the proposed ploidy
            final List<Double> logProbabilitiesCopyNumberProductStates = copyNumberProductStates.stream()
                    .map(cnps ->
                            data.copyRatioLogDensity(
                                    si,
                                    copyRatioNormalization * calculateTotalCopyNumber(populationFractions, cnps, normalPloidyState) / Math.max(EPSILON, initialPloidy),
                                    copyRatioNoiseConstant)
                            + cnps.stream().mapToDouble(ploidyStatePrior::logProbability).sum())
                    .collect(Collectors.toList());
            //find maximum-a-posteriori copy-number product state using copy-ratio--only posteriors
            final int maxPosteriorCopyNumberProductStateIndex = IntStream.range(0, copyNumberProductStates.size()).boxed()
                    .max((i, j) -> Double.compare(logProbabilitiesCopyNumberProductStates.get(i), logProbabilitiesCopyNumberProductStates.get(j))).get();
            final List<Integer> copyNumberProductState = copyNumberProductStates.get(maxPosteriorCopyNumberProductStateIndex);

            //for all ploidy-state product states consistent with the sampled copy-number product state, calculate the copy-ratio--minor-allele-fraction posteriors
            final List<List<PloidyState>> ploidyStateProductStates =
                    new ArrayList<>(Sets.cartesianProduct(copyNumberProductState.stream().map(ploidyStateSetsMap::get).collect(Collectors.toList())));
            final List<PloidyState> ploidyStateProductState;
            if (ploidyStateProductStates.size() == 1) {
                //if only one consistent ploidy-state product, not necessary to perform maximum a posteriori calculations
                ploidyStateProductState = ploidyStateProductStates.get(0);
            } else {
                //calculate the copy ratio of the sampled copy-number product state
                final double copyRatio = copyRatioNormalization * calculateTotalCopyNumber(populationFractions, copyNumberProductState, normalPloidyState) / Math.max(EPSILON, initialPloidy);
                //for all possible ploidy-state product states, calculate the copy-ratio--minor-allele-fraction posterior
                final List<Double> logProbabilitiesPloidyStateProductStates = ploidyStateProductStates.stream()
                        .map(psps ->
                                data.logDensity(
                                        si,
                                        copyRatio,
                                        calculateMinorAlleleFraction(populationFractions, psps, normalPloidyState),
                                        copyRatioNoiseConstant)
                                + psps.stream().mapToDouble(ploidyStatePrior::logProbability).sum())
                        .collect(Collectors.toList());
                //find maximum-a-posteriori ploidy-state product state using copy-ratio--minor-allele-fraction posteriors
                final int maxPosteriorPloidyStateProductStateIndex = IntStream.range(0, ploidyStateProductStates.size()).boxed()
                        .max((i, j) -> Double.compare(logProbabilitiesPloidyStateProductStates.get(i), logProbabilitiesPloidyStateProductStates.get(j))).get();
                ploidyStateProductState = ploidyStateProductStates.get(maxPosteriorPloidyStateProductStateIndex);
            }

            //store the sampled ploidy-state product state in the list of variant profiles
            IntStream.range(0, numVariantPopulations).forEach(i -> variantProfiles.get(i).set(si, ploidyStateProductState.get(i)));
        }
        //return a new VariantProfileCollection created from the list of variant profiles
        return new VariantProfileCollection(variantProfiles);
    }

    private static double calculateTotalCopyNumber(final PopulationMixture.PopulationFractions populationFractions,
                                                   final List<Integer> totalCopyNumberProductState,
                                                   final PloidyState normalPloidyState) {
        final int numVariantPopulations = populationFractions.numVariantPopulations();
        //loop performs better than stream
        double totalCopyNumber = 0.;
        for (int populationIndex = 0; populationIndex < numVariantPopulations; populationIndex++) {
            totalCopyNumber += totalCopyNumberProductState.get(populationIndex) * populationFractions.get(populationIndex);
        }
        totalCopyNumber += normalPloidyState.total() * populationFractions.normalFraction();
        return totalCopyNumber;
    }

    private static double calculateMinorAlleleFraction(final PopulationMixture.PopulationFractions populationFractions,
                                                       final List<PloidyState> ploidyStateProductState,
                                                       final PloidyState normalPloidyState) {
        final int numVariantPopulations = populationFractions.numVariantPopulations();
        //performance is not as critical here, so we use streams
        final double mAlleleCopyNumber = IntStream.range(0, numVariantPopulations)
                .mapToDouble(i -> ploidyStateProductState.get(i).m() * populationFractions.get(i))
                .sum() + normalPloidyState.m() * populationFractions.normalFraction();
        final double nAlleleCopyNumber = IntStream.range(0, numVariantPopulations)
                .mapToDouble(i -> ploidyStateProductState.get(i).n() * populationFractions.get(i))
                .sum() + normalPloidyState.n() * populationFractions.normalFraction();
        return calculateMinorAlleleFraction(mAlleleCopyNumber, nAlleleCopyNumber);
    }

    private static double calculateMinorAlleleFraction(final double m, final double n) {
        return Math.min(m, n) / Math.max(EPSILON, m + n);
    }

    private static boolean isOutsideBounds(final TumorHeterogeneityState state,
                                           final TumorHeterogeneityData data) {
        final double ploidyMax = data.priors().ploidyStatePrior().maxCopyNumber();
        return state.concentration() < CONCENTRATION_MIN || state.concentration() > CONCENTRATION_MAX ||
                state.copyRatioNormalization() < COPY_RATIO_NORMALIZATION_MIN || state.copyRatioNoiseConstant() > COPY_RATIO_NORMALIZATION_MAX ||
                state.copyRatioNoiseConstant() < COPY_RATIO_NOISE_CONSTANT_MIN || state.copyRatioNoiseConstant() > COPY_RATIO_NOISE_CONSTANT_MAX ||
                state.ploidy() < PLOIDY_MIN || state.ploidy() > ploidyMax ||
                state.populationMixture().variantProfileCollection().stream().anyMatch(PopulationMixture.VariantProfile::isTotalDeletion);
    }
}
