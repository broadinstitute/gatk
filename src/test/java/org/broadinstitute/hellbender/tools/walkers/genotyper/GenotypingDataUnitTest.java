package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

/**
 * Test {@link org.broadinstitute.gatk.tools.walkers.genotyper.InfiniteRandomMatingPopulationModel}
 */
public final class GenotypingDataUnitTest {

    @Test(dataProvider="ploidyAndMaximumAlleleAndReadCountsData")
    public void testInstantiation(final int[] ploidies, final int[] readCounts) {
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = ReadLikelihoodsUnitTester.readLikelihoods(2, readCounts);
        final SampleList sampleList = likelihoods;
        final PloidyModel ploidyModel = new HeterogeneousPloidyModel(sampleList, ploidies);
        final GenotypingData<Allele> data = new GenotypingData<>(ploidyModel, likelihoods);
        Assert.assertEquals(data.asListOfAlleles(), likelihoods.asListOfAlleles());
        Assert.assertEquals(data.asListOfSamples(), likelihoods.asListOfSamples());
        Assert.assertEquals(data.readLikelihoods(), likelihoods);
        Assert.assertEquals(data.ploidyModel(), ploidyModel);
        for (int i = 0; i < data.numberOfSamples(); i++) {
            Assert.assertEquals(data.indexOfSample(data.getSample(i)), i);
        }

        for (int i = 0; i < data.numberOfAlleles(); i++) {
            Assert.assertEquals(data.indexOfAllele(data.getAllele(i)), i);
        }
    }

    private static final int[][] PLOIDIES = {
            {1, 1, 1, 1},
            {1, 2, 3, 4},
            {2, 2, 2, 2},
            {2, 1, 2, 1},
            {1},
            {2},
            {},
    };


    private static final int[][] READ_COUNTS = {
            { 10 , 100, 50, 20 },
            { 0, 100, 10, 1 },
            { 1, 2, 3, 4 },
            { 10, 20, 50, 40 },
            { 10  },
            { 20 },
            { }
    };

    @DataProvider(name="ploidyAndMaximumAlleleAndReadCountsData")
    public Object[][] ploidyAndMaximumAlleleAndReadCountsData() {
        final List<Object[]> result = new ArrayList<>(PLOIDIES.length * 2);
        for (int i = 0; i < PLOIDIES.length; i++)
            result.add(new Object[] {PLOIDIES[i], READ_COUNTS[i]});
        return result.toArray(new Object[0][]);
    }

}