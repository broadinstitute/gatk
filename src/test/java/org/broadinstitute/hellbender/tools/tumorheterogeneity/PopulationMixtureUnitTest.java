package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyStatePrior;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.mcmc.posteriorsummary.DecileCollection;
import org.broadinstitute.hellbender.utils.mcmc.posteriorsummary.PosteriorSummary;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.IntStream;

/**
 * Unit tests for {@link PopulationMixture}.  Checks that tests of state validity are correctly performed.
 */
public class PopulationMixtureUnitTest {
    private static final double EPSILON = 1E-10;
    private static final PosteriorSummary DUMMY_POSTERIOR_SUMMARY = new PosteriorSummary(Double.NaN, Double.NaN, Double.NaN);
    private static final DecileCollection DUMMY_DECILE_COLLECTION =
            new DecileCollection(Collections.singletonList(Double.NaN), DecileCollection.ConstructionMode.SAMPLES);
    private static final PloidyState NORMAL_PLOIDY_STATE = new PloidyState(1, 1);
    private static final PopulationMixture.VariantProfile DUMMY_VARIANT_PROFILE = new PopulationMixture.VariantProfile(Collections.singletonList(new PloidyState(0, 0)));
    private static final TumorHeterogeneityData DATA;
    private static final PloidyStatePrior DUMMY_PLOIDY_STATE_PRIOR;
    private static final double DUMMY_PARAMETER = 1.;
    private static final TumorHeterogeneityPriorCollection DUMMY_PRIORS;

    static {
        DUMMY_POSTERIOR_SUMMARY.setDeciles(DUMMY_DECILE_COLLECTION);
        //need valid segment-mean posterior summary to construct TumorHeterogeneityData, but it is not used in tests
        final PosteriorSummary segmentMeanPosteriorSummary = new PosteriorSummary(0., -0.1, 0.1);
        segmentMeanPosteriorSummary.setDeciles(new DecileCollection(Arrays.asList(0., -0.1, 0.1), DecileCollection.ConstructionMode.SAMPLES));
        final ACNVModeledSegment segment1 = new ACNVModeledSegment(new SimpleInterval("1", 1, 25), segmentMeanPosteriorSummary, DUMMY_POSTERIOR_SUMMARY);
        final ACNVModeledSegment segment2 = new ACNVModeledSegment(new SimpleInterval("1", 26, 100), segmentMeanPosteriorSummary, DUMMY_POSTERIOR_SUMMARY);
        final Map<PloidyState, Double> unnormalizedLogProbabilityMassFunctionMap = new HashMap<>();
        unnormalizedLogProbabilityMassFunctionMap.put(NORMAL_PLOIDY_STATE, 0.);
        DUMMY_PLOIDY_STATE_PRIOR = new PloidyStatePrior(unnormalizedLogProbabilityMassFunctionMap);
        DUMMY_PRIORS = new TumorHeterogeneityPriorCollection(
                NORMAL_PLOIDY_STATE, DUMMY_PLOIDY_STATE_PRIOR,
                DUMMY_PARAMETER, DUMMY_PARAMETER, DUMMY_PARAMETER, DUMMY_PARAMETER, DUMMY_PARAMETER, DUMMY_PARAMETER, DUMMY_PARAMETER, DUMMY_PARAMETER);
        DATA = new TumorHeterogeneityData(Arrays.asList(segment1, segment2), DUMMY_PRIORS);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testSinglePopulation() {
        //fail if only one population (must have at least one variant and one normal)
        final PopulationMixture.PopulationFractions populationFractions = new PopulationMixture.PopulationFractions(Collections.singletonList(1.));
        final PopulationMixture.VariantProfileCollection variantProfileCollection = new PopulationMixture.VariantProfileCollection(Collections.emptyList());
        new PopulationMixture(populationFractions, variantProfileCollection, NORMAL_PLOIDY_STATE);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testUnnormalizedPopulationFractions() {
        //fail if population fractions not normalized to unity
        final PopulationMixture.PopulationFractions populationFractions = new PopulationMixture.PopulationFractions(Arrays.asList(0.1, 0.1));
        final PopulationMixture.VariantProfileCollection variantProfileCollection = new PopulationMixture.VariantProfileCollection(Collections.singletonList(DUMMY_VARIANT_PROFILE));
        new PopulationMixture(populationFractions, variantProfileCollection, NORMAL_PLOIDY_STATE);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNoVariants() {
        //fail if number of variant populations is not positive
        final PopulationMixture.PopulationFractions populationFractions = new PopulationMixture.PopulationFractions(Arrays.asList(0.1, 0.9));
        final PopulationMixture.VariantProfileCollection variantProfileCollection = new PopulationMixture.VariantProfileCollection(Collections.emptyList());
        new PopulationMixture(populationFractions, variantProfileCollection, NORMAL_PLOIDY_STATE);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInconsistentNumbersOfPopulationsAndVariantPopulations() {
        //fail if number of variant populations + 1 is not equal to number of population fractions
        final PopulationMixture.PopulationFractions populationFractions = new PopulationMixture.PopulationFractions(Arrays.asList(0.2, 0.1, 0.7));
        final PopulationMixture.VariantProfileCollection variantProfileCollection = new PopulationMixture.VariantProfileCollection(Collections.singletonList(DUMMY_VARIANT_PROFILE));
        new PopulationMixture(populationFractions, variantProfileCollection, NORMAL_PLOIDY_STATE);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDifferentNumberOfSegmentsAcrossVariants() {
        //fail if number of segments is not the same for all variant populations
        final PopulationMixture.PopulationFractions populationFractions = new PopulationMixture.PopulationFractions(Arrays.asList(0.2, 0.1, 0.7));
        final PopulationMixture.VariantProfile variantProfile1 = DUMMY_VARIANT_PROFILE;
        final PopulationMixture.VariantProfile variantProfile2 = new PopulationMixture.VariantProfile(Arrays.asList(new PloidyState(0, 1), new PloidyState(0, 0)));
        final PopulationMixture.VariantProfileCollection variantProfileCollection = new PopulationMixture.VariantProfileCollection(Arrays.asList(variantProfile1, variantProfile2));
        new PopulationMixture(populationFractions, variantProfileCollection, NORMAL_PLOIDY_STATE);
    }

    @Test
    public void testPloidyCalculation() {
        final PopulationMixture.PopulationFractions populationFractions = new PopulationMixture.PopulationFractions(Arrays.asList(0.2, 0.1, 0.7));
        final PopulationMixture.VariantProfile variantProfile1 = new PopulationMixture.VariantProfile(Arrays.asList(new PloidyState(1, 2), new PloidyState(0, 0)));
        final PopulationMixture.VariantProfile variantProfile2 = new PopulationMixture.VariantProfile(Arrays.asList(new PloidyState(0, 0), new PloidyState(0, 1)));
        final PopulationMixture.VariantProfileCollection variantProfileCollection = new PopulationMixture.VariantProfileCollection(Arrays.asList(variantProfile1, variantProfile2));
        Assert.assertEquals(new PopulationMixture(populationFractions, variantProfileCollection, NORMAL_PLOIDY_STATE).ploidy(DATA), 1.625, EPSILON);
    }

    @Test
    public void testCollapseNormalPopulations() {
        final PopulationMixture.PopulationFractions populationFractions = new PopulationMixture.PopulationFractions(Arrays.asList(0.2, 0.1, 0.7));
        final PopulationMixture.VariantProfile variantProfile1 = new PopulationMixture.VariantProfile(Arrays.asList(new PloidyState(1, 2), new PloidyState(0, 0)));
        final PopulationMixture.VariantProfile variantProfile2 = new PopulationMixture.VariantProfile(Arrays.asList(NORMAL_PLOIDY_STATE, NORMAL_PLOIDY_STATE));
        final PopulationMixture.VariantProfileCollection variantProfileCollection = new PopulationMixture.VariantProfileCollection(Arrays.asList(variantProfile1, variantProfile2));
        final PopulationMixture collapsedPopulationMixture = new PopulationMixture(populationFractions, variantProfileCollection, NORMAL_PLOIDY_STATE).collapseNormalPopulations(NORMAL_PLOIDY_STATE);
        Assert.assertTrue(collapsedPopulationMixture.variantProfileCollection().get(1).isNormal(NORMAL_PLOIDY_STATE));
        IntStream.range(0, populationFractions.size()).forEach(i ->
                Assert.assertEquals(collapsedPopulationMixture.populationFraction(i), Arrays.asList(0.2, 0., 0.8).get(i), EPSILON));
    }

    @Test
    public void testNormalAndTumorFractions() {
        final PopulationMixture.PopulationFractions populationFractions = new PopulationMixture.PopulationFractions(Arrays.asList(0.2, 0.1, 0.7));
        final PopulationMixture.VariantProfile variantProfile1 = new PopulationMixture.VariantProfile(Arrays.asList(new PloidyState(1, 2), new PloidyState(0, 0)));
        final PopulationMixture.VariantProfile variantProfile2 = new PopulationMixture.VariantProfile(Arrays.asList(new PloidyState(0, 0), new PloidyState(0, 1)));
        final PopulationMixture.VariantProfileCollection variantProfileCollection = new PopulationMixture.VariantProfileCollection(Arrays.asList(variantProfile1, variantProfile2));
        final PopulationMixture populationMixture = new PopulationMixture(populationFractions, variantProfileCollection, NORMAL_PLOIDY_STATE);
        Assert.assertEquals(populationMixture.populationFractions().normalFraction(), 0.7, EPSILON);
        Assert.assertEquals(populationMixture.populationFractions().tumorFraction(), 0.3, EPSILON);
    }
}
