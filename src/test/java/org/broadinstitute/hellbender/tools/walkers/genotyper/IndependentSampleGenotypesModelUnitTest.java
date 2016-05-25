package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Test {@link IndependentSampleGenotypesModel}
 */
public final class IndependentSampleGenotypesModelUnitTest {

    @Test(dataProvider="ploidyAndMaximumAlleleAndReadCountsData")
    public void testCalculateLikelihoods(final int[] ploidies, final int alleleCount, final int discardAlleleCount, final int[] readCounts) {
        final ReadLikelihoods<Allele> likelihoods = ReadLikelihoodsUnitTester.readLikelihoods(alleleCount, readCounts);
        final AlleleList<Allele> genotypingAlleleList = discardAlleleCount == 0 ? likelihoods : discardAllelesAtRandom(likelihoods,discardAlleleCount);
        final SampleList sampleList = SampleListUnitTester.sampleList(ploidies.length);
        final PloidyModel ploidyModel = new HeterogeneousPloidyModel(sampleList,ploidies);
        final GenotypingData<Allele> data = new GenotypingData<>(ploidyModel,likelihoods);
        final IndependentSampleGenotypesModel model = new IndependentSampleGenotypesModel();
        final GenotypingLikelihoods<Allele> gLikelihoods = model.calculateLikelihoods(genotypingAlleleList,data);
        Assert.assertNotNull(gLikelihoods);
        AlleleListUnitTester.assertAlleleList(gLikelihoods, genotypingAlleleList.asListOfAlleles());
        SampleListUnitTester.assertSampleList(gLikelihoods,sampleList.asListOfSamples());
        final int sampleCount = gLikelihoods.numberOfSamples();
        for (int i = 0; i < sampleCount; i++) {
            final GenotypeLikelihoods sampleLikelihoods = gLikelihoods.sampleLikelihoods(i);
            Assert.assertNotNull(sampleLikelihoods);
            final double[] values = sampleLikelihoods.getAsVector();
            Assert.assertNotNull(values);
            Assert.assertEquals(values.length, new GenotypeLikelihoodCalculators().getInstance(ploidies[i], genotypingAlleleList.numberOfAlleles()).genotypeCount());
            for (int j = 0; j < values.length; j++)
                Assert.assertTrue(values[j] <= 0);
        }
    }

    private AlleleList<Allele> discardAllelesAtRandom(final AlleleList<Allele> likelihoods, final int discardAlleleCount) {
        final Random rnd = Utils.getRandomGenerator();
        final List<Allele> subset = new ArrayList<>(likelihoods.asListOfAlleles());
        for (int i = 0; i < discardAlleleCount; i++) {
            subset.remove(rnd.nextInt(subset.size()));
        }
        return new IndexedAlleleList<>(subset);
    }

    /**
     * Each entry contains to value, where the first is the total number of alleles and the second
     * The number to discard some arbitrary number of alleles for genotyping for the {@link #testCalculateLikelihoods}.
     */
    private static final int[][] ALLELE_COUNTS = {
            {1, 0},
            {2, 1},
            {5, 2},
            {10, 4},
            {1, 0},
            {2, 1},
            {10, 7}
    };

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
        for (int i = 0; i < PLOIDIES.length; i++) {
            result.add(new Object[] {PLOIDIES[i], ALLELE_COUNTS[i][0], 0, READ_COUNTS[i]});
            final int discardAlleleCount = ALLELE_COUNTS[i][1];
            if (discardAlleleCount == 0) continue;
            result.add(new Object[] { PLOIDIES[i], ALLELE_COUNTS[i][0], ALLELE_COUNTS[i][1], READ_COUNTS[i]});
        }
        return result.toArray(new Object[0][]);
    }

}