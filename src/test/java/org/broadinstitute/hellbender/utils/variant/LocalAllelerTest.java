package org.broadinstitute.hellbender.utils.variant;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class LocalAllelerTest extends BaseTest {

    @DataProvider
    public Object[][] getTestCasesLAAandLGT(){
        final Allele ALT_AAT = Allele.create("AAT");
        final VariantContextBuilder vcbNoNonRef = getVcBuilder().alleles(Arrays.asList(Allele.REF_A, Allele.ALT_C, Allele.ALT_T, ALT_AAT));
        final VariantContextBuilder vcbWithNonRef = getVcBuilder().alleles(Arrays.asList(Allele.REF_A, Allele.ALT_C, Allele.ALT_T, ALT_AAT, Allele.NON_REF_ALLELE));

        // GT 0/2 LAA [2]   LGT 0/1
        Genotype het0_2 = makeGenotypeWithAlleles(Arrays.asList(Allele.REF_A, Allele.ALT_T));

        // GT 0/0 LAA []    LGT 0/0
        Genotype hom0_0 = makeGenotypeWithAlleles(Arrays.asList(Allele.REF_A, Allele.REF_A));

        // GT 1/2 LAA [1,2] LGT 1/2
        Genotype het1_2 = makeGenotypeWithAlleles(Arrays.asList(Allele.ALT_C, Allele.ALT_T));

        // GT 2/2 LAA [2]   LGT 1/1
        Genotype hom2_2 = makeGenotypeWithAlleles(Arrays.asList(Allele.ALT_T, Allele.ALT_T));

        // GT 1/3 LAA [1,3] LGT 1/2
        Genotype het1_3 = makeGenotypeWithAlleles(Arrays.asList(Allele.ALT_C, ALT_AAT));


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

    private VariantContextBuilder getVcBuilder() {
        return new VariantContextBuilder("handmade", "chr1", 100, 100, Collections.emptyList());
    }

    private Genotype makeGenotypeWithAlleles(List<Allele> alleles) {
        GenotypeBuilder gb = new GenotypeBuilder("sample");
        return gb.alleles(alleles).make();
    }

    private static Object[] makeGenotypes(VariantContextBuilder vcb, Genotype het0_2, List<Integer> LAA, String LGT) {
        return new Object[]{vcb.genotypes(het0_2).make(), het0_2, new GenotypeBuilder(het0_2)
                .attribute("LAA", LAA).attribute("LGT", LGT).make()};
    }

    @Test(dataProvider = "getTestCasesLAAandLGT")
    public void testAddLocalFields(VariantContext vc, Genotype original, Genotype expected) {
        Genotype actual = LocalAlleler.addLocalFields(original, vc);
        VariantContextTestUtils.assertGenotypesAreEqual(actual, expected);
    }

    @Test
    public void testADToLAD(){

    }
}