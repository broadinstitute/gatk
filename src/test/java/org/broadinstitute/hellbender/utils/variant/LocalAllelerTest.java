package org.broadinstitute.hellbender.utils.variant;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.testng.annotations.Test;

import java.util.Arrays;

import static org.testng.Assert.*;

public class LocalAllelerTest extends BaseTest {

    public Object[][] getTestCases(){
        VariantContextBuilder vcBuilder = new VariantContextBuilder();
        vcBuilder.alleles(Arrays.asList(Allele.REF_A, Allele.ALT_C, Allele.ALT_T, Allele.create("AAT")));
        GenotypeBuilder genotypeBuilder = new GenotypeBuilder("s1", Arrays.asList(Allele.REF_A, Allele.ALT_T));
        new Object[][]{
                genotypeBuilder.make(), genotypeBuilder.attributes("LAA", "")
        }
    }

    @Test
    public void testAddLocalFields() {
    }
}