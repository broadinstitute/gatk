package org.broadinstitute.hellbender.utils.variant;

import com.google.common.primitives.Ints;
import htsjdk.variant.variantcontext.*;
import net.sf.cglib.core.Local;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import shaded.cloud_nio.com.google.errorprone.annotations.Var;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class LocalAllelerTest extends BaseTest {

    private static final Allele ALT_AAT = Allele.create("AAT");
    private  static final VariantContextBuilder vcbNoNonRef = getVcBuilder().alleles(Arrays.asList(Allele.REF_A, Allele.ALT_C, Allele.ALT_T, ALT_AAT));
    private static final VariantContextBuilder vcbWithNonRef = getVcBuilder().alleles(Arrays.asList(Allele.REF_A, Allele.ALT_C, Allele.ALT_T, ALT_AAT, Allele.NON_REF_ALLELE));

    // GT 0/2 LAA [2]   LGT 0/1
    private static final Genotype het0_2 = makeGenotypeWithAlleles(Arrays.asList(Allele.REF_A, Allele.ALT_T));

    // GT 0/0 LAA []    LGT 0/0
    private static final Genotype hom0_0 = makeGenotypeWithAlleles(Arrays.asList(Allele.REF_A, Allele.REF_A));

    // GT 1/2 LAA [1,2] LGT 1/2
    private static final Genotype het1_2 = makeGenotypeWithAlleles(Arrays.asList(Allele.ALT_C, Allele.ALT_T));

    // GT 2/2 LAA [2]   LGT 1/1
    private static final Genotype hom2_2 = makeGenotypeWithAlleles(Arrays.asList(Allele.ALT_T, Allele.ALT_T));

    // GT 1/3 LAA [1,3] LGT 1/2
    private static final Genotype het1_3 = makeGenotypeWithAlleles(Arrays.asList(Allele.ALT_C, ALT_AAT));

    // GT 0 LAA [0] LGT 0
    private static final Genotype hap_ref = makeGenotypeWithAlleles(Arrays.asList(Allele.REF_A));

    // GT 3 LAA [1] LGT 1
    private static final Genotype hap_alt = makeGenotypeWithAlleles(Arrays.asList(ALT_AAT));

    @DataProvider
    public Object[][] getTestCasesLAAandLGT(){
        return new Object[][]{
                makeGenotypes(vcbNoNonRef, het0_2, Arrays.asList(2), "0/1"),
                makeGenotypes(vcbNoNonRef, hom0_0, Arrays.asList(), "0/0"),
                makeGenotypes(vcbNoNonRef, het1_2, Arrays.asList(1,2), "1/2"),
                makeGenotypes(vcbNoNonRef, hom2_2, Arrays.asList(2), "1/1"),
                makeGenotypes(vcbNoNonRef, het1_3, Arrays.asList(1,3), "1/2"),

                makeGenotypes(vcbWithNonRef, het0_2, Arrays.asList(2,4), "0/1"),
                makeGenotypes(vcbWithNonRef, hom0_0, Arrays.asList(4), "0/0"),
                makeGenotypes(vcbWithNonRef, het1_2, Arrays.asList(1,2,4), "1/2"),
                makeGenotypes(vcbWithNonRef, hom2_2, Arrays.asList(2,4), "1/1"),
                makeGenotypes(vcbWithNonRef, het1_3, Arrays.asList(1,3,4), "1/2")
        };
    }

    private static VariantContextBuilder getVcBuilder() {
        return new VariantContextBuilder("handmade", "chr1", 100, 100, Collections.emptyList());
    }

    private static Genotype makeGenotypeWithAlleles(List<Allele> alleles) {
        GenotypeBuilder gb = new GenotypeBuilder("sample");
        return gb.alleles(alleles).make();
    }

    private static Object[] makeGenotypes(VariantContextBuilder vcb, Genotype rootGenotype, List<Integer> LAA, String LGT){
        return makeGenotypes(vcb, rootGenotype, LAA, LGT, Collections.emptyList(), Collections.emptyList(), Collections.emptyList(), Collections.emptyList());
    }

    private static Object[] makeGenotypesWithAD(VariantContextBuilder vcb, Genotype rootGenotype,
                                                List<Integer> LAA, String LGT,
                                                List<Integer> AD, List<Integer> LAD) {
        return makeGenotypes(vcb, rootGenotype, LAA, LGT, AD, LAD, Collections.emptyList(), Collections.emptyList());
    }

    private static Object[] makeGenotypesWithPL(VariantContextBuilder vcb, Genotype rootGenotype,
                                                List<Integer> LAA, String LGT,
                                                List<Integer> PL, List<Integer> LPL) {
        return makeGenotypes(vcb, rootGenotype, LAA, LGT, Collections.emptyList(), Collections.emptyList(), PL, LPL);
    }


    private static Object[] makeGenotypes(VariantContextBuilder vcb, Genotype rootGenotype,
                                                List<Integer> LAA, String LGT,
                                                List<Integer> AD, List<Integer> LAD,
                                                List<Integer> PL, List<Integer> LPL) {
        GenotypeBuilder gb = new GenotypeBuilder(rootGenotype);
        if( ! AD.isEmpty()){
            gb.AD(Ints.toArray(AD));
        }

        if(!PL.isEmpty()){
            gb.PL(Ints.toArray(PL));
        }

        Genotype newRoot = gb.make();

        gb.attribute(LocalAlleler.LAA, LAA)
                .attribute(LocalAlleler.LGT, LGT);

        if(!LAD.isEmpty()){
            gb.attribute(LocalAlleler.LAD, LAD);
        }

        if(!LPL.isEmpty()){
            gb.attribute(LocalAlleler.LPL, Ints.toArray(LPL));
        }

        return new Object[]{vcb.genotypes(newRoot).make(), newRoot, gb.make()};
    }

    @Test(dataProvider = "getTestCasesLAAandLGT")
    public void testAddLocalFields(VariantContext vc, Genotype original, Genotype expected) {
        Genotype actual = LocalAlleler.addLocalFields(original, vc);
        VariantContextTestUtils.assertGenotypesAreEqual(actual, expected);
    }

    @DataProvider
    public Object[][] getTestCasesLAAandLGTandAD(){
        return new Object[][]{
                makeGenotypesWithAD(vcbNoNonRef, het0_2, Arrays.asList(2), "0/1", Arrays.asList(0,1,2,3), Arrays.asList(0,2)),
                makeGenotypesWithAD(vcbNoNonRef, hom0_0, Arrays.asList(), "0/0", Arrays.asList(0,1,2,3), Arrays.asList(0)),
                makeGenotypesWithAD(vcbNoNonRef, het1_2, Arrays.asList(1,2), "1/2", Arrays.asList(0,1,2,3), Arrays.asList(0,1,2)),
                makeGenotypesWithAD(vcbNoNonRef, hom2_2, Arrays.asList(2), "1/1", Arrays.asList(0,1,2,3), Arrays.asList(0,2)),
                makeGenotypesWithAD(vcbNoNonRef, het1_3, Arrays.asList(1,3), "1/2", Arrays.asList(0,1,2,3), Arrays.asList(0,1,3)),
                makeGenotypesWithAD(vcbNoNonRef, hap_ref, Arrays.asList(), "0", Arrays.asList(0,1,2,3), Arrays.asList(0)),
                makeGenotypesWithAD(vcbNoNonRef, hap_alt, Arrays.asList(3), "1", Arrays.asList(0,1,2,3), Arrays.asList(0,3)),

                makeGenotypesWithAD(vcbWithNonRef, het0_2, Arrays.asList(2,4), "0/1", Arrays.asList(0,1,2,3,4), Arrays.asList(0,2,4)),
                makeGenotypesWithAD(vcbWithNonRef, hom0_0, Arrays.asList(4), "0/0", Arrays.asList(0,1,2,3,4), Arrays.asList(0,4)),
                makeGenotypesWithAD(vcbWithNonRef, het1_2, Arrays.asList(1,2,4), "1/2", Arrays.asList(0,1,2,3,4), Arrays.asList(0,1,2,4)),
                makeGenotypesWithAD(vcbWithNonRef, hom2_2, Arrays.asList(2,4), "1/1", Arrays.asList(0,1,2,3,4), Arrays.asList(0,2,4)),
                makeGenotypesWithAD(vcbWithNonRef, het1_3, Arrays.asList(1,3,4), "1/2", Arrays.asList(0,1,2,3,4), Arrays.asList(0,1,3,4)),
                makeGenotypesWithAD(vcbWithNonRef, hap_ref, Arrays.asList(4), "0", Arrays.asList(0,1,2,3,4), Arrays.asList(0,4)),
                makeGenotypesWithAD(vcbWithNonRef, hap_alt, Arrays.asList(3,4), "1", Arrays.asList(0,1,2,3,4), Arrays.asList(0,3,4)),
        };
    }

    @Test(dataProvider = "getTestCasesLAAandLGTandAD")
    public void testADToLAD(VariantContext vc, Genotype original, Genotype expected){
        Genotype actual = LocalAlleler.addLocalFields(original, vc);
        VariantContextTestUtils.assertGenotypesAreEqual(actual, expected);
    }

    @DataProvider
    public Object[][] getTestCasesLAAandLGTandPL(){
        return new Object[][]{
                makeGenotypesWithPL(vcbNoNonRef, het0_2, Arrays.asList(2), "0/1", Arrays.asList(0,1,11,2,12,22,3,13,23,33), Arrays.asList(0,2,22)),
                makeGenotypesWithPL(vcbNoNonRef, hom0_0, Arrays.asList(), "0/0", Arrays.asList(0,1,11,2,12,22,3,13,23,33), Arrays.asList(0)),
                makeGenotypesWithPL(vcbNoNonRef, het1_2, Arrays.asList(1,2), "1/2", Arrays.asList(0,1,11,2,12,22,3,13,23,33), Arrays.asList(0,1,11,2,12,22)),
                makeGenotypesWithPL(vcbNoNonRef, hom2_2, Arrays.asList(2), "1/1", Arrays.asList(0,1,11,2,12,22,3,13,23,33), Arrays.asList(0,2,22)),
                makeGenotypesWithPL(vcbNoNonRef, het1_3, Arrays.asList(1,3), "1/2", Arrays.asList(0,1,11,2,12,22,3,13,23,33), Arrays.asList(0,1,11,3,13,33)),
                makeGenotypesWithPL(vcbNoNonRef, hap_ref, Arrays.asList(), "0", Arrays.asList(0,1,2,3), Arrays.asList(0)),
                makeGenotypesWithPL(vcbNoNonRef, hap_alt, Arrays.asList(3), "1", Arrays.asList(0,1,2,3), Arrays.asList(0,3)),

                makeGenotypesWithPL(vcbWithNonRef, het0_2, Arrays.asList(2,4), "0/1", Arrays.asList(0,1,11,2,12,22,3,13,23,33,4,14,24,34,44), Arrays.asList(0,2,22,4,24,44)),
                makeGenotypesWithPL(vcbWithNonRef, hom0_0, Arrays.asList(4), "0/0", Arrays.asList(0,1,11,2,12,22,3,13,23,33,4,14,24,34,44), Arrays.asList(0,4,44)),
                makeGenotypesWithPL(vcbWithNonRef, het1_2, Arrays.asList(1,2,4), "1/2", Arrays.asList(0,1,11,2,12,22,3,13,23,33,4,14,24,34,44), Arrays.asList(0,1,11,2,12,22,4,14,24,44)),
                makeGenotypesWithPL(vcbWithNonRef, hom2_2, Arrays.asList(2,4), "1/1", Arrays.asList(0,1,11,2,12,22,3,13,23,33,4,14,24,34,44), Arrays.asList(0,2,22,4,24,44)),
                makeGenotypesWithPL(vcbWithNonRef, het1_3, Arrays.asList(1,3,4), "1/2", Arrays.asList(0,1,11,2,12,22,3,13,23,33,4,14,24,34,44), Arrays.asList(0,1,11,3,13,33,4,14,34,44)),
                makeGenotypesWithPL(vcbWithNonRef, hap_ref, Arrays.asList(4), "0", Arrays.asList(0,1,2,3,4), Arrays.asList(0,4)),
                makeGenotypesWithPL(vcbWithNonRef, hap_alt, Arrays.asList(3,4), "1", Arrays.asList(0,1,2,3,4), Arrays.asList(0,3,4)),
        };
    }

    @Test(dataProvider = "getTestCasesLAAandLGTandPL")
    public void testPLToLPL(VariantContext vc, Genotype original, Genotype expected){
        Genotype actual = LocalAlleler.addLocalFields(original, vc);
        VariantContextTestUtils.assertGenotypesAreEqual(actual, expected);
    }

    @Test
    public void testRemoveNonLocal(){
        Object[] values = makeGenotypes(vcbWithNonRef, het1_3, Arrays.asList(1, 3, 4), "1/2",
                Arrays.asList(0, 1, 2, 3, 4), Arrays.asList(0, 1, 3, 4),
                Arrays.asList(0, 1, 11, 2, 12, 22, 3, 13, 23, 33, 4, 14, 24, 34, 44), Arrays.asList(0, 1, 11, 3, 13, 33, 4, 14, 34, 44));
        VariantContext vc = (VariantContext)values[0];
        Genotype originalGenotype = (Genotype)values[1];
        Genotype expectedGenotype = (Genotype)values[2];
        Genotype localizedGenotype = LocalAlleler.addLocalFields(originalGenotype, vc, true, false);
        Assert.assertTrue(localizedGenotype.getAlleles().isEmpty());
        Assert.assertFalse(localizedGenotype.hasAD());
        Assert.assertFalse(localizedGenotype.hasPL());

        Assert.assertEquals(localizedGenotype.getExtendedAttribute(LocalAlleler.LGT), expectedGenotype.getExtendedAttribute(LocalAlleler.LGT));
        Assert.assertEquals(localizedGenotype.getExtendedAttribute(LocalAlleler.LAA), expectedGenotype.getExtendedAttribute(LocalAlleler.LAA));
        Assert.assertEquals(localizedGenotype.getExtendedAttribute(LocalAlleler.LAD), expectedGenotype.getExtendedAttribute(LocalAlleler.LAD));
        Assert.assertEquals(localizedGenotype.getExtendedAttribute(LocalAlleler.LPL), expectedGenotype.getExtendedAttribute(LocalAlleler.LPL));
    }
}