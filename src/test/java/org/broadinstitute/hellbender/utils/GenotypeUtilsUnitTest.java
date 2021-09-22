package org.broadinstitute.hellbender.utils;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GenotypeUtilsUnitTest extends GATKBaseTest {
    private static final Allele Aref = Allele.create("A", true);
    private static final Allele T = Allele.create("T");
    private static final Allele TT = Allele.create("TT");
    private static final Allele TA = Allele.create("TA");
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

    @DataProvider
    public Object[][] getEdgeCases() {
        //not an edge case, but the other test only has homVar genotypes
        int[] biallelicRef_PLs = {0, 12, 139};
        int[] biallelicHet_PLs = {102, 0, 373};
        final Genotype g1 = new GenotypeBuilder("sample1", Arrays.asList(Aref,Aref)).PL(biallelicRef_PLs).make();
        final Genotype g2 = new GenotypeBuilder("sample1", Arrays.asList(Aref, T)).PL(biallelicHet_PLs).make();
        GenotypeCounts c1 = new GenotypeCounts(1.0, 0, 0);
        GenotypeCounts c2 = new GenotypeCounts(0, 1.0, 0);
        final VariantContext vc_biallelic = new VariantContextBuilder("in memory", "1", 100, 100,
                Arrays.asList(Aref, T)).make();


        int[] GQ0Ref_PLs = {0, 0, 437};
        int[] GQ0Var_PLs = {437, 0, 0};
        int[] noCall_PLs = {0, 0, 0};
        final Genotype g3 = new GenotypeBuilder("sample1", Arrays.asList(Aref,Aref)).PL(GQ0Ref_PLs).make();
        final Genotype g4 = new GenotypeBuilder("sample1", Arrays.asList(Aref, T)).PL(GQ0Var_PLs).make();
        final Genotype g5 = new GenotypeBuilder("sample1", Arrays.asList(Allele.NO_CALL, Allele.NO_CALL)).PL(noCall_PLs).make();
        GenotypeCounts c3 = new GenotypeCounts(1.0, 0, 0);
        GenotypeCounts c4 = new GenotypeCounts(0, 1.0, 0);
        GenotypeCounts c5 = new GenotypeCounts(0, 0, 0);

        int[] multiAlleleicRef_PLs = {0, 3, 705, 3, 705, 705};
        int[] multiAllelicHet1_PLs = {184, 0, 179, 202, 200, 403};
        int[] multiAllelicHet2_PLs = {168, 210, 524, 0, 314, 293};
        int[] hetNonRef_PLs = {952, 462, 426, 489, 0, 456};
        final Genotype g6 = new GenotypeBuilder("sample1", Arrays.asList(Aref,Aref)).PL(multiAlleleicRef_PLs).make();
        final Genotype g7 = new GenotypeBuilder("sample1", Arrays.asList(Aref, T)).PL(multiAllelicHet1_PLs).make();
        final Genotype g8 = new GenotypeBuilder("sample1", Arrays.asList(Allele.NO_CALL, Allele.NO_CALL)).PL(multiAllelicHet2_PLs).make();
        final Genotype g9 = new GenotypeBuilder("sample1", Arrays.asList(T, TT)).PL(hetNonRef_PLs).make();
        GenotypeCounts c6 = new GenotypeCounts(1.0, 0, 0);
        GenotypeCounts c7 = new GenotypeCounts(0, 1.0, 0);
        GenotypeCounts c8 = new GenotypeCounts(0, 1.0, 0);
        GenotypeCounts c9 = new GenotypeCounts(0, 0, 1.0);
        final VariantContext vc_twoAlts = new VariantContextBuilder("in memory", "1", 100, 100,
                Arrays.asList(Aref, T, TT)).make();

        int[] threeAltsRef_PLs = {0, 99, 1485, 99, 1485, 1485, 99, 1485, 1485, 1485};
        int[] threeAltsHet3_PLs = {676, 729, 2236, 729, 2236, 2236, 0, 1508, 1508, 1453};
        int[] threeAltsHetNonRef_PLs = {631, 339, 319, 254, 0, 227, 591, 337, 252, 589};
        final Genotype g10 = new GenotypeBuilder("sample1", Arrays.asList(Aref,Aref)).PL(threeAltsRef_PLs).make();
        final Genotype g11 = new GenotypeBuilder("sample1", Arrays.asList(Aref, TA)).PL(threeAltsHet3_PLs).make();
        final Genotype g12 = new GenotypeBuilder("sample1", Arrays.asList(T, TT)).PL(threeAltsHetNonRef_PLs).make();
        GenotypeCounts c10 = new GenotypeCounts(1.0, 0, 0);
        GenotypeCounts c11 = new GenotypeCounts(0, 1.0, 0);
        GenotypeCounts c12 = new GenotypeCounts(0, 0, 1.0);
        final VariantContext vc_threeAlts = new VariantContextBuilder("in memory", "1", 100, 100,
                Arrays.asList(Aref, T, TT, TA)).make();

        return new Object[][]{
                {vc_biallelic, g1, c1},
                {vc_biallelic, g2, c2},
                {vc_biallelic, g3, c3},
                {vc_biallelic, g4, c4},
                {vc_biallelic, g5, c5},
                {vc_twoAlts, g6, c6},
                {vc_twoAlts, g7, c7},
                {vc_twoAlts, g8, c8},
                {vc_twoAlts, g9, c9},
                {vc_threeAlts, g10, c10},
                {vc_threeAlts, g11, c11},
                {vc_threeAlts, g12, c12}
        };
    }

    //Make sure GQ0 samples get counted only once
    @Test (dataProvider = "getEdgeCases")
    public void testEdgeCasesRounded(VariantContext vc, Genotype g, GenotypeCounts expected) {
        boolean rounded = true;
        final ArrayList<Genotype> genotypesArray = new ArrayList<>();
        genotypesArray.add(g);
        final GenotypesContext genotypes = GenotypesContext.create(genotypesArray);

        final GenotypeCounts actual = GenotypeUtils.computeDiploidGenotypeCounts(vc, genotypes, rounded);
        Assert.assertEquals(actual.getRefs(), expected.getRefs(), DELTA_PRECISION);
        Assert.assertEquals(actual.getHets(), expected.getHets(), DELTA_PRECISION);
        Assert.assertEquals(actual.getHoms(), expected.getHoms(), DELTA_PRECISION);
    }

    @Test
    public void testMultiallelicFakePLs() {
        final Genotype g  = VariantContextTestUtils.makeG("sample1", Aref, Aref, 17);
        final double[] fakePLs = GenotypeUtils.makeApproximateDiploidLog10LikelihoodsFromGQ(g, 3);
        Assert.assertTrue(fakePLs[0] == 0);
        Assert.assertTrue(fakePLs.length == GenotypeLikelihoods.numLikelihoods(3, 2));
        Assert.assertTrue(fakePLs[1] > fakePLs[fakePLs.length-1]);  //het likelihood should be less than hom var likelihood (but great in log10-sapce)
    }

}