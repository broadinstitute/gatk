package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.coordinates.SimplexPosition;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Represents a mixture of {@link VariantProfile} populations with corresponding population fractions.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class PopulationMixture {
    private final int numPopulations;   //variant populations + normal population
    private final int numVariantPopulations;   //variant populations + normal population
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
        numPopulations = this.populationFractions.numPopulations();
        numVariantPopulations = this.populationFractions.numVariantPopulations();
        Utils.validateArg(numPopulations == this.variantProfileCollection.numVariantPopulations + 1,
                "Number of populations must be equal to number of variant populations + 1.");
        numSegments = this.variantProfileCollection.numSegments;
        this.normalPloidyState = normalPloidyState;
    }

    public int numPopulations() {
        return numPopulations;
    }

    public int numVariantPopulations() {
        return numVariantPopulations;
    }

    public int numSegments() {
        return numSegments;
    }

    public double populationFraction(final int populationIndex) {
        validatePopulationIndex(populationIndex, numPopulations);
        return populationFractions.get(populationIndex);
    }

    public PopulationFractions populationFractions() {
        return populationFractions;
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

    public PloidyState normalPloidyState() {
        return normalPloidyState;
    }

    /**
     * Return a new {@link PopulationMixture}, with the population fractions
     * for variant profiles that are completely normal reassigned to the normal population.
     */
    public PopulationMixture collapseNormalPopulations(final PloidyState normalPloidyState) {
        final List<Double> adjustedPopulationFractions = new ArrayList<>(numPopulations);
        double normalFraction = populationFractions.get(numPopulations - 1);
        for (int populationIndex = 0; populationIndex < numPopulations - 1; populationIndex++) {
            final double populationFraction = populationFractions.get(populationIndex);
            if (variantProfileCollection.get(populationIndex).isNormal(normalPloidyState)) {
                normalFraction += populationFraction;
                adjustedPopulationFractions.add(0.);
            } else {
                adjustedPopulationFractions.add(populationFraction);
            }
        }
        adjustedPopulationFractions.add(normalFraction);
        return new PopulationMixture(new PopulationFractions(adjustedPopulationFractions), variantProfileCollection, normalPloidyState);
    }

    @Override
    public String toString() {
        return populationFractions.toString();
    }

    /*===============================================================================================================*
     * METHODS                                                                                                       *
     *===============================================================================================================*/

    public double ploidy(final TumorHeterogeneityData data) {
        validateData(data, numSegments);
        return IntStream.range(0, numSegments).mapToDouble(i -> calculateLengthAndPopulationAveragedCopyNumberFunction(i, PloidyState::total, data)).sum();

    }

    double calculatePopulationAveragedCopyNumberFunction(final int segmentIndex,
                                                         final Function<PloidyState, Integer> copyNumberFunction) {
        validateSegmentIndex(segmentIndex, numSegments);

        return IntStream.range(0, numPopulations)
                .mapToDouble(populationIndex -> populationFraction(populationIndex) * calculateCopyNumberFunction(segmentIndex, populationIndex, copyNumberFunction))
                .sum();
    }

    private int calculateCopyNumberFunction(final int segmentIndex,
                                            final Integer populationIndex,
                                            final Function<PloidyState, Integer> copyNumberFunction) {
        return isNormalPopulation(populationIndex) ?
                copyNumberFunction.apply(normalPloidyState) :
                copyNumberFunction.apply(ploidyState(populationIndex, segmentIndex));
    }

    private double calculateLengthAndPopulationAveragedCopyNumberFunction(final int segmentIndex,
                                                                          final Function<PloidyState, Integer> copyNumberFunction,
                                                                          final TumorHeterogeneityData data) {
        return data.fractionalLength(segmentIndex) * calculatePopulationAveragedCopyNumberFunction(segmentIndex, copyNumberFunction);
    }

    /*===============================================================================================================*
     * INNER CLASSES                                                                                                 *
     *===============================================================================================================*/

    public static final class PopulationFractions extends SimplexPosition {
        //list of doubles normalized to unity, size = number of populations (including normal), i-th element = fraction of population i by cell number
        //normal population is last element
        private static final long serialVersionUID = 723562368L;
        private final int numPopulations;
        private final int numVariantPopulations;

        public PopulationFractions(final List<Double> populationFractions) {
            super(populationFractions);
            numPopulations = populationFractions.size();
            numVariantPopulations = numPopulations - 1;
        }

        public int numPopulations() {
            return numPopulations;
        }

        public int numVariantPopulations() {
            return numVariantPopulations;
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
            super(Collections.unmodifiableList(Utils.nonNull(variantProfiles).stream().map(VariantProfile::new).collect(Collectors.toList())));
            Utils.validateArg(variantProfiles.size() > 0, "Number of variants must be positive.");
            final int numSegmentsForFirstVariant = variantProfiles.get(0).numSegments;
            Utils.validateArg(numSegmentsForFirstVariant > 0, "Number of segments must be positive.");
            Utils.validateArg(variantProfiles.stream().allMatch(vp -> vp.numSegments == numSegmentsForFirstVariant),
                    "Number of segments must be same for all variants.");
            numSegments = numSegmentsForFirstVariant;
            numVariantPopulations = variantProfiles.size();
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

        public boolean equals(final VariantProfileCollection other) {
            for (int populationIndex = 0; populationIndex < numVariantPopulations; populationIndex++) {
                for (int segmentIndex = 0; segmentIndex < numSegments; segmentIndex++) {
                    if (!get(populationIndex).get(segmentIndex).equals(other.get(populationIndex).get(segmentIndex))) {
                        return false;
                    }
                }
            }
            return true;
        }
    }

    public static class VariantProfile extends ArrayList<PloidyState> {
        //list of integers, size = number of segments, i-th element = ploidy state of segment i
        private static final long serialVersionUID = 78476L;
        private static final PloidyState HOMOZYGOUS_DELETION_PLOIDY_STATE = new PloidyState(0, 0);
        private final int numSegments;

        public VariantProfile(final List<PloidyState> ploidyStates) {
            super(Collections.unmodifiableList(new ArrayList<>(Utils.nonNull(ploidyStates))));
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

        public boolean isNormal(final PloidyState normalPloidyState) {
            return stream().allMatch(ps -> ps.equals(normalPloidyState));
        }

        public boolean isTotalDeletion() {
            return stream().allMatch(ps -> ps.equals(HOMOZYGOUS_DELETION_PLOIDY_STATE));
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
