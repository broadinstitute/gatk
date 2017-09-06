package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public final class GenotypeLikelihoodCalculatorsUnitTest extends BaseTest {

    GenotypeLikelihoodCalculators calcs;

    @BeforeClass
    public void init(){
        calcs = new GenotypeLikelihoodCalculators();
    }


    @DataProvider(name="genotypeCount")
    public Iterator<Object[]> ploidyAndMaximumAlleleAndReadCountsData() {
        final int[][] expected =
       {{0, 1, 1,  1,  1,  1,   1 },
        {0, 1, 2,  3,  4,  5,   6 },
        {0, 1, 3,  6, 10, 15,  21 },
        {0, 1, 4, 10, 20, 35,  56 },
        {0, 1, 5, 15, 35, 70, 126 }};

        final List<Object[]> result = new ArrayList<>();
        for (int i = 0; i < expected.length; i++){
            for (int j = 0; j < expected[i].length; j++){
                result.add(new Object[]{i, j, expected[i][j]});
            }
        }

        //now add it backwards too just to check ordering issues
        for (int i = expected.length-1; i >= 0; i--){
            for (int j = expected[i].length-1; j >= 0; j--){
                result.add(new Object[]{i, j, expected[i][j]});
            }
        }

        return result.iterator();
    }
    @Test(dataProvider = "genotypeCount")
    public void testGenotypeCountSharedInstance(int ploidy, int alleleCount, int expected) throws Exception {
        Assert.assertEquals(calcs.genotypeCount(ploidy, alleleCount), expected);
    }

    @Test(dataProvider = "genotypeCount")
    public void testInstanceSharedInstance(int ploidy, int alleleCount, int expected) throws Exception {
        if (ploidy > 0) {
            final GenotypeLikelihoodCalculator inst = calcs.getInstance(ploidy, alleleCount);
            Assert.assertEquals(inst.genotypeCount(), expected);
            Assert.assertEquals(inst.ploidy(), ploidy);
            Assert.assertEquals(inst.alleleCount(), alleleCount);
        }
    }

    @Test(dataProvider = "genotypeCount")
    public void testGenotypeCountNewInstance(int ploidy, int alleleCount, int expected) throws Exception {
        Assert.assertEquals(new GenotypeLikelihoodCalculators().genotypeCount(ploidy, alleleCount), expected);
    }

    @Test(dataProvider = "genotypeCount")
    public void testInstanceNewInstance(int ploidy, int alleleCount, int expected) throws Exception {
        if (ploidy > 0) {
            final GenotypeLikelihoodCalculator inst = new GenotypeLikelihoodCalculators().getInstance(ploidy, alleleCount);
            Assert.assertEquals(inst.genotypeCount(), expected);
            Assert.assertEquals(inst.ploidy(), ploidy);
            Assert.assertEquals(inst.alleleCount(), alleleCount);
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testGenotypeCountOverflow() throws Exception {
        final int genotypeCount = new GenotypeLikelihoodCalculators().genotypeCount(10_000, 10_000);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testGenotypeCountNegativePloidy() throws Exception {
        new GenotypeLikelihoodCalculators().genotypeCount(-1, 1);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testGenotypeCountNegativeAlleleCount() throws Exception {
        new GenotypeLikelihoodCalculators().genotypeCount(1, -1);
    }

    @Test
    public void testComputeMaxAcceptableAlleleCount(){
        Assert.assertEquals(1024, GenotypeLikelihoodCalculators.computeMaxAcceptableAlleleCount(1, 1024));
        Assert.assertEquals(44, GenotypeLikelihoodCalculators.computeMaxAcceptableAlleleCount(2, 1024));
        Assert.assertEquals(17, GenotypeLikelihoodCalculators.computeMaxAcceptableAlleleCount(3, 1024));
        Assert.assertEquals(5, GenotypeLikelihoodCalculators.computeMaxAcceptableAlleleCount(10, 1024));
        Assert.assertEquals(3, GenotypeLikelihoodCalculators.computeMaxAcceptableAlleleCount(20, 1024));
        Assert.assertEquals(2, GenotypeLikelihoodCalculators.computeMaxAcceptableAlleleCount(100, 1024));
    }
}
