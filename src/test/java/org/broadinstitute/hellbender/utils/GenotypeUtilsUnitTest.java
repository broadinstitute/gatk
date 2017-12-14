package org.broadinstitute.hellbender.utils;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GenotypeUtilsUnitTest extends GATKBaseTest {
    private static final Allele Aref = Allele.create("A", true);
    private static final Allele T = Allele.create("T");
    private static final double DELTA_PRECISION = .001;


    @DataProvider
    public Object[][] getGTsWithAndWithoutLikelihoods(){
        final List<Allele> DIPLOID_ALLELES = Arrays.asList(Aref, T);
        final List<Allele> HAPLOID_ALLELES = Arrays.asList(Aref);
        final int[] SOME_PLS = {0, 4, 1};
        final int[] HAPLOID_PL = {0, 4};
        return new Object[][]{
                {getGenotypeBuilder().noPL().alleles(DIPLOID_ALLELES).make(), false},
                {getGenotypeBuilder().noPL().alleles(HAPLOID_ALLELES).make(), false},
                {getGenotypeBuilder().PL(HAPLOID_PL).alleles(HAPLOID_ALLELES).make(), false},
                {getGenotypeBuilder().PL(SOME_PLS).alleles(DIPLOID_ALLELES).make(), true},
                {GenotypeBuilder.createMissing("sample", 2), false}
        };
    }

    private static GenotypeBuilder getGenotypeBuilder() {
        return new GenotypeBuilder("sample", Arrays.asList(Aref, T));
    }

    @Test(dataProvider = "getGTsWithAndWithoutLikelihoods")
    public void testIsDiploidWithLikelihoods(Genotype g, boolean expected) throws Exception {
        Assert.assertEquals(GenotypeUtils.isDiploidWithLikelihoods(g), expected);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIsDiploidWithLikelihoodsWithNull() {
        GenotypeUtils.isDiploidWithLikelihoods(null);
    }



    @DataProvider
    public Object[][] getGenotypeCountsParameters(){
        final VariantContext vc = new VariantContextBuilder("in memory", "1", 100, 100,
                                                                        Arrays.asList(Aref, T)).make();
        final ArrayList<Genotype> genotypesArray = new ArrayList<>();
        for(int i = 0; i<100; i++){
            final Genotype g = new GenotypeBuilder("sample" + i, Arrays.asList(T, T)).PL(new int[]{100,10,0}).make();
            genotypesArray.add(g);
        }
        final GenotypesContext genotypes = GenotypesContext.create(genotypesArray);

        return new Object[][]{
                {vc, genotypes, false, new GenotypeCounts(0.000, 9.091, 90.909)},
                {vc, genotypes, true, new GenotypeCounts(0.000, 0.000, 100.000)},
        };
    }

    @Test(dataProvider = "getGenotypeCountsParameters")
    public void testRounding(VariantContext vc, GenotypesContext gt, boolean round, GenotypeCounts expected) {
        final GenotypeCounts actual = GenotypeUtils.computeDiploidGenotypeCounts(vc, gt, round);
        Assert.assertEquals(actual.getRefs(), expected.getRefs(), DELTA_PRECISION);
        Assert.assertEquals(actual.getHets(), expected.getHets(), DELTA_PRECISION);
        Assert.assertEquals(actual.getHoms(), expected.getHoms(), DELTA_PRECISION);
    }


}