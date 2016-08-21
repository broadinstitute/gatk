package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.mcmc.DecileCollection;
import org.broadinstitute.hellbender.utils.mcmc.PosteriorSummary;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;

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

    static {
        DUMMY_POSTERIOR_SUMMARY.setDeciles(DUMMY_DECILE_COLLECTION);
        //need valid segment-mean posterior summary to construct TumorHeterogeneityData, but it is not used in tests
        final PosteriorSummary segmentMeanPosteriorSummary = new PosteriorSummary(0., -0.1, 0.1);
        segmentMeanPosteriorSummary.setDeciles(new DecileCollection(Arrays.asList(0., -0.1, 0.1), DecileCollection.ConstructionMode.SAMPLES));
        final ACNVModeledSegment segment1 = new ACNVModeledSegment(new SimpleInterval("1", 1, 25), segmentMeanPosteriorSummary, DUMMY_POSTERIOR_SUMMARY);
        final ACNVModeledSegment segment2 = new ACNVModeledSegment(new SimpleInterval("1", 26, 100), segmentMeanPosteriorSummary, DUMMY_POSTERIOR_SUMMARY);
        DATA = new TumorHeterogeneityData(Arrays.asList(segment1, segment2));
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
    public void testCalculatePopulationAndGenomicAveragedPloidy() {
        final PopulationMixture.PopulationFractions populationFractions = new PopulationMixture.PopulationFractions(Arrays.asList(0.2, 0.1, 0.7));
        final PopulationMixture.VariantProfile variantProfile1 = new PopulationMixture.VariantProfile(Arrays.asList(new PloidyState(1, 2), new PloidyState(0, 0)));
        final PopulationMixture.VariantProfile variantProfile2 = new PopulationMixture.VariantProfile(Arrays.asList(new PloidyState(0, 0), new PloidyState(0, 1)));
        final PopulationMixture.VariantProfileCollection variantProfileCollection = new PopulationMixture.VariantProfileCollection(Arrays.asList(variantProfile1, variantProfile2));
        Assert.assertEquals(new PopulationMixture(populationFractions, variantProfileCollection, NORMAL_PLOIDY_STATE).ploidy(DATA), 1.625, EPSILON);
    }
}