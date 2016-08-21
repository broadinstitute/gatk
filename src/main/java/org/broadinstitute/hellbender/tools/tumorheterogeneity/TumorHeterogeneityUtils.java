package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import com.google.common.collect.Sets;
import org.apache.commons.math3.analysis.function.Logit;
import org.apache.commons.math3.analysis.function.Sigmoid;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.special.Gamma;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyStatePrior;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class TumorHeterogeneityUtils {
    static final double EPSILON = 1E-10;

    private TumorHeterogeneityUtils() {}

    static double calculateLogPosterior(final TumorHeterogeneityState state,
                                        final TumorHeterogeneityData data) {
        final int numPopulations = state.populationMixture().numPopulations();
        final int numSegments = data.numSegments();

        //concentration prior
        final double concentrationPriorAlpha = state.priors().concentrationPriorHyperparameterValues().getAlpha();
        final double concentrationPriorBeta = state.priors().concentrationPriorHyperparameterValues().getBeta();
        final double concentration = state.concentration();
        final double logPriorConcentration =
                concentrationPriorAlpha * Math.log(concentrationPriorBeta + EPSILON)
                        + (concentrationPriorAlpha - 1.) * Math.log(concentration + EPSILON)
                        - concentrationPriorBeta * concentration
                        - Gamma.logGamma(concentrationPriorAlpha);

        //copy-ratio noise-floor prior
        final double copyRatioNoiseFloorPriorAlpha = state.priors().copyRatioNoiseFloorPriorHyperparameterValues().getAlpha();
        final double copyRatioNoiseFloorPriorBeta = state.priors().copyRatioNoiseFloorPriorHyperparameterValues().getBeta();
        final double copyRatioNoiseFloor = state.copyRatioNoiseFloor();
        final double logPriorCopyRatioNoiseFloor =
                copyRatioNoiseFloorPriorAlpha * Math.log(copyRatioNoiseFloorPriorBeta + EPSILON)
                        + (copyRatioNoiseFloorPriorAlpha - 1.) * Math.log(copyRatioNoiseFloor + EPSILON)
                        - copyRatioNoiseFloorPriorBeta * copyRatioNoiseFloor
                        - Gamma.logGamma(copyRatioNoiseFloorPriorAlpha);

        //copy-ratio noise-factor prior
        final double copyRatioNoiseFactorPriorAlpha = state.priors().copyRatioNoiseFactorPriorHyperparameterValues().getAlpha();
        final double copyRatioNoiseFactorPriorBeta = state.priors().copyRatioNoiseFactorPriorHyperparameterValues().getBeta();
        final double copyRatioNoiseFactor = state.copyRatioNoiseFactor();
        final double logPriorCopyRatioNoiseFactor =
                copyRatioNoiseFactorPriorAlpha * Math.log(copyRatioNoiseFactorPriorBeta + EPSILON)
                        + (copyRatioNoiseFactorPriorAlpha - 1.) * Math.log(copyRatioNoiseFactor - 1. + EPSILON)
                        - copyRatioNoiseFactorPriorBeta * (copyRatioNoiseFactor - 1.)
                        - Gamma.logGamma(copyRatioNoiseFactorPriorAlpha);

        //minor-allele-fraction noise-factor prior
        final double minorAlleleFractionNoiseFactorPriorAlpha = state.priors().minorAlleleFractionNoiseFactorPriorHyperparameterValues().getAlpha();
        final double minorAlleleFractionNoiseFactorPriorBeta = state.priors().minorAlleleFractionNoiseFactorPriorHyperparameterValues().getBeta();
        final double minorAlleleFractionNoiseFactor = state.minorAlleleFractionNoiseFactor();
        final double logPriorMinorAlleleFractionNoiseFactor =
                minorAlleleFractionNoiseFactorPriorAlpha * Math.log(minorAlleleFractionNoiseFactorPriorBeta + EPSILON)
                        + (minorAlleleFractionNoiseFactorPriorAlpha - 1.) * Math.log(minorAlleleFractionNoiseFactor - 1. + EPSILON)
                        - minorAlleleFractionNoiseFactorPriorBeta * (minorAlleleFractionNoiseFactor - 1.)
                        - Gamma.logGamma(minorAlleleFractionNoiseFactorPriorAlpha);

        //population-fractions prior
        final double logPriorPopulationFractionsSum = IntStream.range(0, numPopulations)
                .mapToDouble(i -> (concentration - 1.) * Math.log(state.populationMixture().populationFraction(i) + EPSILON))
                .sum();
        final double logPriorPopulationFractions =
                Gamma.logGamma(concentration * numPopulations)
                        - numPopulations * Gamma.logGamma(concentration)
                        + logPriorPopulationFractionsSum;

        //variant-profiles prior
        final double logPriorVariantProfiles = calculateLogPriorVariantProfiles(
                state.populationMixture().variantProfileCollection(),
                state.priors().ploidyStatePrior(),
                data.segmentLengths());

        //copy-ratio--minor-allele-fraction likelihood
        final double logLikelihoodSegments = calculateLogLikelihoodSegments(state, data);

//        System.out.println(logPriorConcentration + " " + logPriorCopyRatioNoiseFloor + " " + logPriorCopyRatioNoiseFactor + " " + logPriorMinorAlleleFractionNoiseFactor + " " +
//                logPriorPopulationFractions + " " + logPriorVariantProfiles + " " + logLikelihoodSegments);

        return logPriorConcentration + logPriorCopyRatioNoiseFloor + logPriorCopyRatioNoiseFactor + logPriorMinorAlleleFractionNoiseFactor +
                logPriorPopulationFractions + logPriorVariantProfiles + logLikelihoodSegments;
    }

    static double calculateLogPriorVariantProfiles(final PopulationMixture.VariantProfileCollection variantProfileCollection,
                                                   final PloidyStatePrior ploidyStatePrior,
                                                   final List<Integer> segmentLengths) {
        double logPriorVariantProfiles = 0.;
        for (int populationIndex = 0; populationIndex < variantProfileCollection.numVariantPopulations(); populationIndex++) {
            for (int segmentIndex = 0; segmentIndex < variantProfileCollection.numSegments(); segmentIndex++) {
                final PloidyState ploidyState = variantProfileCollection.ploidyState(populationIndex, segmentIndex);
                logPriorVariantProfiles += segmentLengths.get(segmentIndex) * ploidyStatePrior.logProbability(ploidyState);
            }
        }
        return logPriorVariantProfiles;
    }

    static double calculateLogLikelihoodSegments(final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
        double logLikelihoodSegments = 0.;
        final double ploidy = state.populationMixture().ploidy(data);
        for (int segmentIndex = 0; segmentIndex < data.numSegments(); segmentIndex++) {
            final double totalCopyNumber = state.populationMixture().calculatePopulationAveragedCopyNumberFunction(segmentIndex, PloidyState::total);
            final double mAlleleCopyNumber = state.populationMixture().calculatePopulationAveragedCopyNumberFunction(segmentIndex, PloidyState::m);
            final double nAlleleCopyNumber = state.populationMixture().calculatePopulationAveragedCopyNumberFunction(segmentIndex, PloidyState::n);
            final double copyRatio = totalCopyNumber / (ploidy + EPSILON);
            final double minorAlleleFraction = calculateMinorAlleleFraction(mAlleleCopyNumber, nAlleleCopyNumber);
            logLikelihoodSegments += data.logDensity(segmentIndex, copyRatio, minorAlleleFraction, state.copyRatioNoiseFloor(), state.copyRatioNoiseFactor(), state.minorAlleleFractionNoiseFactor());
        }
        return logLikelihoodSegments;
    }

    private static double calculateMinorAlleleFraction(final double m, final double n) {
        return Math.min(m, n) / (m + n + EPSILON);
    }

    static PopulationMixture.VariantProfileCollection proposeVariantProfileCollection(final RandomGenerator rng,
                                                                                      final TumorHeterogeneityState currentState,
                                                                                      final TumorHeterogeneityData data,
                                                                                      final PopulationMixture.PopulationFractions proposedPopulationFractions,
                                                                                      final double proposedInitialPloidy,
                                                                                      final List<List<Integer>> totalCopyNumberProductStates,
                                                                                      final Map<Integer, Set<PloidyState>> ploidyStateSetsMap) {
        final int numPopulations = currentState.populationMixture().numPopulations();
        final int numSegments = data.numSegments();
        final PloidyState normalPloidyState = currentState.priors().normalPloidyState();
        final PloidyStatePrior ploidyStatePrior = currentState.priors().ploidyStatePrior();
        final List<PopulationMixture.VariantProfile> variantProfiles = new ArrayList<>(Collections.nCopies(numPopulations - 1,
                new PopulationMixture.VariantProfile(Collections.nCopies(numSegments, normalPloidyState))));

        for (int segmentIndex = 0; segmentIndex < numSegments; segmentIndex++) {
            final int si = segmentIndex;
            final int segmentLength = data.length(si);
            final double[] log10ProbabilitiesCopyRatio = totalCopyNumberProductStates.stream()
                    .mapToDouble(tcnps -> calculateTotalCopyNumber(proposedPopulationFractions, tcnps, normalPloidyState) / (proposedInitialPloidy + EPSILON))
                    .map(cr -> data.copyRatioLogDensity(si, cr, currentState.copyRatioNoiseFloor(), currentState.copyRatioNoiseFactor()))
                    .map(MathUtils::logToLog10)
                    .toArray();
            final double[] probabilitiesCopyRatio = MathUtils.normalizeFromLog10ToLinearSpace(log10ProbabilitiesCopyRatio);
            final Function<List<Integer>, Double> probabilityFunctionCopyRatio = totalCopyNumberProductState ->
                    probabilitiesCopyRatio[totalCopyNumberProductStates.indexOf(totalCopyNumberProductState)];
            final List<Integer> totalCopyNumberProductState = GATKProtectedMathUtils.randomSelect(totalCopyNumberProductStates, probabilityFunctionCopyRatio, rng);
            final double totalCopyRatio = calculateTotalCopyNumber(proposedPopulationFractions, totalCopyNumberProductState, normalPloidyState) / (proposedInitialPloidy + EPSILON);

            final List<List<PloidyState>> ploidyStateProductStates =
                    new ArrayList<>(Sets.cartesianProduct(totalCopyNumberProductState.stream().map(ploidyStateSetsMap::get).collect(Collectors.toList())));
            final double[] log10ProbabilitiesPloidyStateProductStates = ploidyStateProductStates.stream()
                    .mapToDouble(psps ->
                            data.logDensity(si, totalCopyRatio, calculateMinorAlleleFraction(proposedPopulationFractions, psps, normalPloidyState),
                                    currentState.copyRatioNoiseFloor(), currentState.copyRatioNoiseFactor(), currentState.minorAlleleFractionNoiseFactor())
                                    + psps.stream().mapToDouble(ps -> segmentLength * ploidyStatePrior.logProbability(ps)).sum())
                    .map(MathUtils::logToLog10)
                    .toArray();
            final double[] probabilitiesPloidyStateProductStates = MathUtils.normalizeFromLog10ToLinearSpace(log10ProbabilitiesPloidyStateProductStates);
            final Function<List<PloidyState>, Double> probabilityFunctionProductStates = ploidyStateProductState ->
                    probabilitiesPloidyStateProductStates[ploidyStateProductStates.indexOf(ploidyStateProductState)];
            final List<PloidyState> ploidyStateProductState = GATKProtectedMathUtils.randomSelect(ploidyStateProductStates, probabilityFunctionProductStates, rng);

            IntStream.range(0, numPopulations - 1).forEach(i -> variantProfiles.get(i).set(si, ploidyStateProductState.get(i)));
        }
        return new PopulationMixture.VariantProfileCollection(variantProfiles);
    }

    private static double calculateTotalCopyNumber(final PopulationMixture.PopulationFractions populationFractions,
                                               final List<Integer> totalCopyNumberProductState,
                                               final PloidyState normalPloidyState) {
        final int numPopulations = populationFractions.size();
        return IntStream.range(0, numPopulations - 1).boxed()
                .mapToDouble(i -> totalCopyNumberProductState.get(i) * populationFractions.get(i))
                .sum() + normalPloidyState.total() * populationFractions.normalFraction();
    }

    private static double calculateMinorAlleleFraction(final PopulationMixture.PopulationFractions populationFractions,
                                                   final List<PloidyState> ploidyStateProductState,
                                                   final PloidyState normalPloidyState) {
        final int numPopulations = populationFractions.size();
        final double mAlleleCopyNumber = IntStream.range(0, numPopulations - 1).boxed()
                .mapToDouble(i -> ploidyStateProductState.get(i).m() * populationFractions.get(i))
                .sum() + normalPloidyState.m() * populationFractions.normalFraction();
        final double nAlleleCopyNumber = IntStream.range(0, numPopulations - 1).boxed()
                .mapToDouble(i -> ploidyStateProductState.get(i).n() * populationFractions.get(i))
                .sum() + normalPloidyState.n() * populationFractions.normalFraction();
        return Math.min(mAlleleCopyNumber, nAlleleCopyNumber) / (mAlleleCopyNumber + nAlleleCopyNumber + EPSILON);
    }

    static double calculateLogJacobianFactor(final PopulationMixture.PopulationFractions populationFractions) {
        final List<Double> breakProportions = calculateBreakProportionsFromPopulationFractions(populationFractions);
        return IntStream.range(0, populationFractions.size() - 1).boxed()
                .mapToDouble(i -> Math.log(populationFractions.get(i)) + Math.log(1. - breakProportions.get(i))).sum();

    }

    static List<Double> calculateTransformedPopulationFractionsFromPopulationFractions(final PopulationMixture.PopulationFractions populationFractions) {
        final List<Double> breakProportions = calculateBreakProportionsFromPopulationFractions(populationFractions);
        return calculateTransformedPopulationFractionsFromBreakProportions(breakProportions);
    }

    static PopulationMixture.PopulationFractions calculatePopulationFractionsFromTransformedPopulationFractions(final List<Double> transformedPopulationFractions) {
        final List<Double> breakProportions = calculateBreakProportionsFromTransformedPopulationFractions(transformedPopulationFractions);
        final List<Double> populationFractions = calculatePopulationFractionsFromBreakProportions(breakProportions);
        return new PopulationMixture.PopulationFractions(populationFractions);
    }

    private static List<Double> calculatePopulationFractionsFromBreakProportions(final List<Double> breakProportions) {
        final int numPopulations = breakProportions.size() + 1;
        final List<Double> populationFractions = new ArrayList<>();
        double cumulativeSum = 0.;
        for (int populationIndex = 0; populationIndex < numPopulations - 1; populationIndex++) {
            final double populationFraction = (1. - cumulativeSum) * breakProportions.get(populationIndex);
            populationFractions.add(populationFraction);
            cumulativeSum += populationFraction;
        }
        populationFractions.add(1. - cumulativeSum);
        return new PopulationMixture.PopulationFractions(populationFractions);
    }

    private static List<Double> calculateBreakProportionsFromPopulationFractions(final PopulationMixture.PopulationFractions populationFractions) {
        final int numPopulations = populationFractions.size();
        final List<Double> breakProportions = new ArrayList<>();
        double cumulativeSum = 0.;
        for (int populationIndex = 0; populationIndex < numPopulations - 1; populationIndex++) {
            final double breakProportion = populationFractions.get(populationIndex) / (1. - cumulativeSum);
            breakProportions.add(breakProportion);
            cumulativeSum += populationFractions.get(populationIndex);
        }
        return breakProportions;
    }

    private static List<Double> calculateBreakProportionsFromTransformedPopulationFractions(final List<Double> transformedPopulationFractions) {
        final int numPopulations = transformedPopulationFractions.size() + 1;
        return IntStream.range(0, numPopulations - 1).boxed()
                .map(i -> new Sigmoid().value(transformedPopulationFractions.get(i) + Math.log(1. / (numPopulations - (i + 1)))))
                .collect(Collectors.toList());
    }

    private static List<Double> calculateTransformedPopulationFractionsFromBreakProportions(final List<Double> breakProportions) {
        final int numPopulations = breakProportions.size() + 1;
        return IntStream.range(0, numPopulations - 1).boxed()
                .map(i -> new Logit().value(breakProportions.get(i)) - Math.log(1. / (numPopulations - (i + 1))))
                .collect(Collectors.toList());
    }
}
