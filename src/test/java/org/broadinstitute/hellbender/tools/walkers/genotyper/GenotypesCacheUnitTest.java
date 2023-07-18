package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

public final class GenotypesCacheUnitTest extends GATKBaseTest {


    // ploidy, genotypeIndex, expected allele counts array
    @DataProvider(name = "randomAccessData")
    public Object[][] randomAccessData() {
        return new Object[][] {
                {1, 0, new int[] {1}},
                {1, 7, new int[] {0,0,0,0,0,0,0,1}},
                {2, 4, new int[] {0,1,1}},
                {2, 4, new int[] {0,1,1}},
                {3, 13, new int[] {1,0,1,1}},
                {4, (int) GenotypeIndexCalculator.indexOfFirstGenotypeWithAllele(4,5), new int[] {3,0,0,0,0,1}},
                {4, (int) GenotypeIndexCalculator.indexOfFirstGenotypeWithAllele(4,4), new int[] {3,0,0,0,1}},
                {4, (int) GenotypeIndexCalculator.indexOfFirstGenotypeWithAllele(4,3), new int[] {3,0,0,1}}

        };
    }

    @Test(dataProvider = "randomAccessData")
    public void testCache(final int ploidy, final int genotypeIndex, final int[] expectedAlleleCounts) {
        final GenotypeAlleleCounts gac = GenotypesCache.get(ploidy,genotypeIndex);
        final int distinctCount = (int) Arrays.stream(expectedAlleleCounts).filter(n -> n > 0).count();
        Assert.assertEquals(gac.distinctAlleleCount(), distinctCount);

        for (int n = 0; n < expectedAlleleCounts.length; n++) {
            Assert.assertEquals(gac.alleleCountFor(n), expectedAlleleCounts[n]);
        }

        for (int n = expectedAlleleCounts.length; n < expectedAlleleCounts.length + 3; n++) {
            Assert.assertEquals(gac.alleleCountFor(n), 0);
        }

        final GenotypeAlleleCounts next = gac.next();
        Assert.assertTrue(next.equals(GenotypesCache.get(ploidy, genotypeIndex+1)));
        final GenotypeAlleleCounts nextNext = next.next();

        GenotypeAlleleCounts plusEleven = next.copy();
        for (int i = 0; i < 10; i++) {
            plusEleven = plusEleven.next();
        }
        Assert.assertTrue(plusEleven.equals(GenotypesCache.get(ploidy, genotypeIndex+11)));
        Assert.assertTrue(nextNext.equals(GenotypesCache.get(ploidy, genotypeIndex+2)));
    }

}
