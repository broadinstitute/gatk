package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class PopulationMixture {
    private static final double POPULATION_FRACTION_NORMALIZATION_EPSILON = 1E-4;

    private final int numPopulations;   //variant populations + normal population
    private final int numSegments;

    private final PopulationFractions populationFractions;
    private final VariantProfileCollection variantProfileCollection;
    private final PloidyState normalPloidyState;

    public PopulationMixture(final List<Double> populationFractions,
                             final List<VariantProfile> variantProfileCollection,
                             final PloidyState normalPloidyState) {
        Utils.nonNull(populationFractions);
        Utils.nonNull(variantProfileCollection);
        Utils.nonNull(normalPloidyState);
        this.populationFractions = new PopulationFractions(populationFractions);
        this.variantProfileCollection = new VariantProfileCollection(variantProfileCollection);
        Utils.validateArg(this.populationFractions.numPopulations == this.variantProfileCollection.numVariantPopulations + 1,
                "Number of populations must be equal to number of variant populations + 1.");
        numPopulations = this.populationFractions.numPopulations;
        numSegments = this.variantProfileCollection.numSegments;
        this.normalPloidyState = normalPloidyState;
    }

    public int numPopulations() {
        return numPopulations;
    }

    public int numSegments() {
        return numSegments;
    }

    public double populationFraction(final int populationIndex) {
        validatePopulationIndex(populationIndex, numPopulations);
        return populationFractions.get(populationIndex);
    }

    public PopulationFractions populationFractions() {
        return new PopulationFractions(populationFractions);
    }

    public VariantProfileCollection variantProfileCollection() {
        return new VariantProfileCollection(variantProfileCollection);
    }

    public boolean isNormalPopulation(final int populationIndex) {
        validatePopulationIndex(populationIndex, numPopulations);
        return populationIndex == numPopulations - 1;
    }

    public PloidyState ploidyState(final int populationIndex, final int segmentIndex) {
        validatePopulationIndex(populationIndex, numPopulations);
        validateSegmentIndex(segmentIndex, numSegments);
        return isNormalPopulation(populationIndex) ?
                normalPloidyState :
                variantProfileCollection.ploidyState(populationIndex, segmentIndex);
    }

    /*===============================================================================================================*
     * METHODS                                                                                                       *
     *===============================================================================================================*/

    public int calculateCopyNumberFunction(final int segmentIndex,
                                           final Integer populationIndex,
                                           final Function<PloidyState, Integer> copyNumberFunction) {
        validateSegmentIndex(segmentIndex, numSegments);
        validatePopulationIndex(populationIndex, numPopulations);
        return isNormalPopulation(populationIndex) ?
                copyNumberFunction.apply(normalPloidyState) :
                copyNumberFunction.apply(ploidyState(populationIndex, segmentIndex));
    }

    public double calculatePopulationAveragedCopyNumberFunction(final int segmentIndex,
                                                                final Function<PloidyState, Integer> copyNumberFunction) {
        return calculatePopulationAveragedCopyNumberFunctionExcludingPopulation(null, segmentIndex, copyNumberFunction);
    }

    public double ploidy(final TumorHeterogeneityData data) {
        validateData(data, numSegments);
        return IntStream.range(0, numSegments).mapToDouble(i -> calculatePopulationAndGenomicAveragedCopyNumberFunction(null, i, PloidyState::total, data)).sum();
    }

    public double populationPloidy(final int populationIndex, final TumorHeterogeneityData data) {
        validatePopulationIndex(populationIndex, numPopulations);
        validateData(data, numSegments);
        return isNormalPopulation(populationIndex) ?
                normalPloidyState.total() :
                variantProfileCollection.get(populationIndex).ploidy(data);
    }

    double calculatePopulationAveragedCopyNumberFunctionExcludingPopulation(final Integer populationIndexToExclude,
                                                                            final int segmentIndex,
                                                                            final Function<PloidyState, Integer> copyNumberFunction) {
        validateSegmentIndex(segmentIndex, numSegments);
        final int populationIndexToExcludeValue;
        if (populationIndexToExclude == null) {
            populationIndexToExcludeValue = -1; //this will not exclude any populations
        } else {
            validatePopulationIndex(populationIndexToExclude, numPopulations);
            populationIndexToExcludeValue = populationIndexToExclude;
        }

        return IntStream.range(0, numPopulations)
                .filter(populationIndex -> populationIndex != populationIndexToExcludeValue)
                .mapToDouble(populationIndex -> populationFraction(populationIndex) * calculateCopyNumberFunction(segmentIndex, populationIndex, copyNumberFunction))
                .sum();
    }

    private double calculatePopulationAndGenomicAveragedCopyNumberFunction(final Integer populationIndexToExclude,
                                                                           final int segmentIndex,
                                                                           final Function<PloidyState, Integer> copyNumberFunction,
                                                                           final TumorHeterogeneityData data) {
        validateData(data, numSegments);
        validateSegmentIndex(segmentIndex, numSegments);

        return data.fractionalLength(segmentIndex) * calculatePopulationAveragedCopyNumberFunctionExcludingPopulation(populationIndexToExclude, segmentIndex, copyNumberFunction);
    }

    /*===============================================================================================================*
     * INNER CLASSES                                                                                                 *
     *===============================================================================================================*/

    public static final class PopulationFractions extends ArrayList<Double> {
        //list of doubles, size = number of populations (including normal), i-th element = fraction of population i by cell number
        //normal population is last element
        private static final long serialVersionUID = 79454L;
        private final int numPopulations;
        public PopulationFractions(final List<Double> populationFractions) {
            super(Utils.nonNull(new ArrayList<>(populationFractions)));
            Utils.validateArg(populationFractions.size() > 1, "Number of populations must be strictly greater than 1.");
            final double populationFractionNormalization = populationFractions.stream().mapToDouble(Double::doubleValue).sum();
            Utils.validateArg(Math.abs(1. - populationFractionNormalization) <= POPULATION_FRACTION_NORMALIZATION_EPSILON,
                    "Population fractions must sum to unity.");
            numPopulations = populationFractions.size();
        }

        public double normalFraction() {
            return get(numPopulations - 1);
        }

        public double tumorFraction() {
            return 1. - normalFraction();
        }
    }

    public static final class VariantProfileCollection extends ArrayList<VariantProfile> {
        //list of VariantProfiles, size = number of populations (excluding normal), i-th element = VariantProfile for population i
        private static final long serialVersionUID = 76498L;
        private final int numSegments;
        private final int numVariantPopulations;
        public VariantProfileCollection(final List<VariantProfile> variantProfiles) {
            super(Utils.nonNull(new ArrayList<>()));
            Utils.validateArg(variantProfiles.size() > 0, "Number of variants must be positive.");
            final int numSegmentsForFirstVariant = variantProfiles.get(0).numSegments;
            Utils.validateArg(variantProfiles.stream().map(s -> s.numSegments).allMatch(n -> n == numSegmentsForFirstVariant),
                    "Number of segments must be same for all variants.");
            Utils.validateArg(numSegmentsForFirstVariant > 0, "Number of segments must be positive.");
            numSegments = numSegmentsForFirstVariant;
            numVariantPopulations = variantProfiles.size();
            variantProfiles.forEach(vp -> add(new VariantProfile(vp)));
        }

        public int numSegments() {
            return numSegments;
        }

        public int numVariantPopulations() {
            return numVariantPopulations;
        }

        public PloidyState ploidyState(final int populationIndex, final int segmentIndex) {
            validatePopulationIndex(populationIndex, numVariantPopulations);
            validateSegmentIndex(segmentIndex, numSegments);
            return get(populationIndex).ploidyState(segmentIndex);
        }
    }

    /**
     * For each variant population, represents ploidy states in each segment.
     */
    static class VariantProfile extends ArrayList<PloidyState> {
        //list of integers, size = number of segments, i-th element = ploidy state of segment i
        private static final long serialVersionUID = 78476L;
        private final int numSegments;

        VariantProfile(final List<PloidyState> ploidyStates) {
            super(Utils.nonNull(new ArrayList<>(ploidyStates)));
            Utils.validateArg(ploidyStates.size() > 0, "Number of segments must be positive.");
            numSegments = ploidyStates.size();
        }

        public PloidyState ploidyState(final int segmentIndex) {
            validateSegmentIndex(segmentIndex, numSegments);
            return get(segmentIndex);
        }

        public int numSegments() {
            return numSegments;
        }

        public double ploidy(final TumorHeterogeneityData data) {
            validateData(data, numSegments);
            return IntStream.range(0, numSegments).mapToDouble(i -> data.fractionalLength(i) * get(i).total()).sum();
        }
    }

    /*===============================================================================================================*
     * OTHER PRIVATE METHODS                                                                                         *
     *===============================================================================================================*/

    private static void validatePopulationIndex(final int populationIndex, final int numPopulations) {
        Utils.validateArg(0 <= populationIndex && populationIndex < numPopulations, "Population index out of range.");
    }

    private static void validateSegmentIndex(int segmentIndex, final int numSegments) {
        Utils.validateArg(0 <= segmentIndex && segmentIndex < numSegments, "Segment index out of range.");
    }

    private static void validateData(final TumorHeterogeneityData data, final int numSegments) {
        Utils.nonNull(data);
        Utils.validateArg(data.numSegments() == numSegments,
                "Tumor-heterogeneity state and data collection must have same number of segments.");
    }
}
