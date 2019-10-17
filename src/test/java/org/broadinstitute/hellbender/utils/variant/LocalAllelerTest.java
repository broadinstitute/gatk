package org.broadinstitute.hellbender.utils.variant;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class LocalAllelerTest extends BaseTest {

    @DataProvider
    public Object[][] getTestCases(){
        VariantContextBuilder vcBuilder = new VariantContextBuilder("handmade", "chr1", 100, 100, Collections.emptyList());
        final Allele ALT_AAT = Allele.create("AAT");
        List<Allele> alleles = Arrays.asList(Allele.REF_A, Allele.ALT_C, Allele.ALT_T, ALT_AAT);
        vcBuilder.alleles(alleles);
        GenotypeBuilder genotypeBuilder = new GenotypeBuilder("s1");
        // GT 0/2 LAA [2]   LGT 0/1
        Genotype het0_2 = genotypeBuilder.alleles( Arrays.asList(Allele.REF_A, Allele.ALT_T)).make();

        // GT 0/0 LAA []    LGT 0/0
        Genotype hom0_0 = genotypeBuilder.alleles( Arrays.asList(Allele.REF_A, Allele.REF_A)).make();

        // GT 1/2 LAA [1,2] LGT 1/2
        Genotype het1_2 = genotypeBuilder.alleles( Arrays.asList(Allele.ALT_C, Allele.ALT_T)).make();

        // GT 2/2 LAA [2]   LGT 1/1
        Genotype hom2_2 = genotypeBuilder.alleles( Arrays.asList(Allele.ALT_T, Allele.ALT_T)).make();

        // GT 1/3 LAA [1,3] LGT 1/2
        Genotype het1_3 = genotypeBuilder.alleles( Arrays.asList(Allele.ALT_C, ALT_AAT)).make();
        return new Object[][]{
                makeGenotypes(vcBuilder, het0_2, Arrays.asList(2), "0/1"),
                makeGenotypes(vcBuilder, hom0_0, Arrays.asList(), "0/0"),
                makeGenotypes(vcBuilder, het1_2, Arrays.asList(1,2), "1/2" ),
                makeGenotypes(vcBuilder, hom2_2, Arrays.asList(2), "1/1"),
                makeGenotypes(vcBuilder, het1_3, Arrays.asList(1,3), "1/2")
        };
    }

    private static Object[] makeGenotypes(VariantContextBuilder vcb, Genotype het0_2, List<Integer> LAA, String LGT) {
        return new Object[]{vcb.genotypes(het0_2).make(), het0_2, new GenotypeBuilder(het0_2).attribute("LAA", LAA).attribute("LGT", LGT).make()};
    }

    @Test(dataProvider = "getTestCases")
    public void testAddLocalFields(VariantContext vc, Genotype original, Genotype expected) {
        Genotype actual = LocalAlleler.addLocalFields(original, vc);
        VariantContextTestUtils.assertGenotypesAreEqual(actual, expected);
    }
}