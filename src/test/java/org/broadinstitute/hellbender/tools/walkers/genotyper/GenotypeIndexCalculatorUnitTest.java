package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.apache.commons.math3.exception.MathArithmeticException;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class GenotypeIndexCalculatorUnitTest {

    // ploidy, allele, index of first genotype with that allele
    @DataProvider(name = "firstGenotypeWithAlleleData")
    public Object[][] firstGenotypeWithAlleleData() {
        return new Object[][] {
                {1, 0, 0 },
                {1, 5, 5},
                {2, 0, 0 },
                {2, 1, 1 },
                {2, 2, 3},
                {4, 0, 0},
                {4, 1, 1},
                {4, 2, 5},
                {4, 3, 15}
        };
    }

    @Test(dataProvider = "firstGenotypeWithAlleleData")
    public void testIndexOfFirstGenotypeWithAllele(final int ploidy, final int allele, final int expected) {
        Assert.assertEquals(GenotypeIndexCalculator.indexOfFirstGenotypeWithAllele(ploidy, allele), expected);
    }

    // ploidy, allele count
    @DataProvider(name = "genotypeCountData")
    public Object[][] genotypeCountData() {
        return new Object[][] {
                {1, 1},
                {1, 5},
                {2, 1},
                {2, 5},
                {3, 5}
        };
    }

    // a ploidy-P genotype with A alleles can be decomposed as 0 <= N <= P copies of the Ath allele and a ploidy P-N genotype
    // using only A-1 alleles.  If this recursion is satisfied, then the genotype method is correct.
    @Test(dataProvider = "genotypeCountData")
    public void testGenotypeCount(final int ploidy, final int alleleCount) {
        final int direct = GenotypeIndexCalculator.genotypeCount(ploidy, alleleCount);
        if (ploidy == 1) {
            Assert.assertEquals(direct, alleleCount);
        } else {
            // the '1' below is the N = P term
            final int recursive = 1 + new IndexRange(0, ploidy).sumInt(n -> GenotypeIndexCalculator.genotypeCount(ploidy - n, alleleCount - 1));
            Assert.assertEquals(direct, recursive);
        }
    }

    @Test(expectedExceptions = MathArithmeticException.class)
    public void testGenotypeCountOverflow() throws Exception {
        final int genotypeCount = GenotypeIndexCalculator.genotypeCount(10_000, 10_000);
    }

    // alleles list, expected index
    @DataProvider(name = "allelesToIndexData")
    public Object[][] allelesToIndexData() {
        return new Object[][] {
                {new int[] {0}, 0},
                {new int[] {0,0}, 0},
                {new int[] {0,0,0}, 0},
                {new int[] {0,1}, 1},
                {new int[] {1,0}, 1},
                {new int[] {1,1}, 2},
                {new int[] {100}, 100},
                {new int[] {0,100}, 5050},
                {new int[] {100,0}, 5050},
                {new int[] {0,1,2}, 5},
                {new int[] {2,0,0}, 4},
                {new int[] {1,2,1}, 6},
                {new int[] {2,1,2}, 8},
                {new int[] {2,2,2}, 9},
        };
    }

    @Test(dataProvider = "allelesToIndexData")
    public void testAllelesToIndex(final int[] alleles, final int index) {
        Assert.assertEquals(GenotypeIndexCalculator.allelesToIndex(alleles), index);
    }

    // allele counts array, expected index
    @DataProvider(name = "alleleCountsToIndexData")
    public Object[][] alleleCountsToIndexData() {
        return new Object[][] {
                {new int[] {0,1}, 0},
                {new int[] {0,2}, 0},
                {new int[] {0,3}, 0},
                {new int[] {0,1,1,1}, 1},
                {new int[] {1,1,0,1}, 1},
                {new int[] {1,2}, 2},
                {new int[] {100,1}, 100},
                {new int[] {0,1,100,1}, 5050},
                {new int[] {100,1,0,1}, 5050},
                {new int[] {0,1,1,1,2,1}, 5},
                {new int[] {2,1,0,2}, 4},
                {new int[] {1,2,2,1}, 6},
                {new int[] {2,2,1,1}, 8},
                {new int[] {2,3}, 9},
        };
    }

    @Test(dataProvider = "alleleCountsToIndexData")
    public void testAlleleCountsToIndex(final int[] counts, final int index) {
        Assert.assertEquals(GenotypeIndexCalculator.alleleCountsToIndex(counts), index);
    }

    @Test
    public void testComputeMaxAcceptableAlleleCount(){
        Assert.assertEquals(1024, GenotypeIndexCalculator.computeMaxAcceptableAlleleCount(1, 1024));
        Assert.assertEquals(44, GenotypeIndexCalculator.computeMaxAcceptableAlleleCount(2, 1024));
        Assert.assertEquals(17, GenotypeIndexCalculator.computeMaxAcceptableAlleleCount(3, 1024));
        Assert.assertEquals(5, GenotypeIndexCalculator.computeMaxAcceptableAlleleCount(10, 1024));
        Assert.assertEquals(3, GenotypeIndexCalculator.computeMaxAcceptableAlleleCount(20, 1024));
        Assert.assertEquals(2, GenotypeIndexCalculator.computeMaxAcceptableAlleleCount(100, 1024));
    }

    // ploidy, new to old allele reordering, expected result
    @DataProvider(name = "newToOldMapData")
    public Object[][] newToOldMapData() {
        return new Object[][] {
                {1, new int[] {0}, new int[] {0}},
                {1, new int[] {1}, new int[] {1}},
                {2, new int[] {0}, new int[] {0}},
                {2, new int[] {1}, new int[] {2}},
                {2, new int[] {2}, new int[] {5}},
                {2, new int[] {0,1}, new int[] {0,1,2}},
                {2, new int[] {1,0}, new int[] {2,1,0}},
                {2, new int[] {0,2}, new int[] {0,3,5}}
        };
    }

    @Test(dataProvider = "newToOldMapData")
    public void testNewToOldIndexMap(final int ploidy, final int[] newToOldAlleleMap, final int[] expected) {
        final int[] result = GenotypeIndexCalculator.newToOldGenotypeMap(ploidy, newToOldAlleleMap);
        Assert.assertEquals(result, expected);
    }

}