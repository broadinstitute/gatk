package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

public class MergeMutect2CallsWithMC3UnitTest extends GATKBaseTest {

    @DataProvider(name = "shouldVariantsBeMatchedDataProvider")
    Object[][] shouldVariantsBeMatchedDataProvider() {
        return new Object[][] {
                //eval contains the truth alt allele
                {Arrays.asList(Allele.REF_A, Allele.ALT_C), Arrays.asList(Allele.REF_A, Allele.ALT_C), true},
                {Arrays.asList(Allele.REF_A, Allele.ALT_C), Arrays.asList(Allele.REF_A, Allele.ALT_C, Allele.ALT_G), true},
                {Arrays.asList(Allele.REF_A, Allele.ALT_C), Arrays.asList(Allele.REF_A, Allele.ALT_G, Allele.ALT_C), true},
                //only the first truth alt allele needs to be contained in eval
                {Arrays.asList(Allele.REF_A, Allele.ALT_C, Allele.ALT_G), Arrays.asList(Allele.REF_A, Allele.ALT_C), true},
                {Arrays.asList(Allele.REF_A, Allele.ALT_C, Allele.ALT_G), Arrays.asList(Allele.REF_A, Allele.ALT_G), false},
                //eval does not contain the truth alt allele
                {Arrays.asList(Allele.REF_A, Allele.ALT_C), Arrays.asList(Allele.REF_A, Allele.ALT_G), false},
                //different reference alleles
                {Arrays.asList(Allele.REF_A, Allele.ALT_C), Arrays.asList(Allele.REF_C, Allele.ALT_A), false}
        };
    }

    @Test(dataProvider = "shouldVariantsBeMatchedDataProvider")
    public void testShouldVariantsBeMatched(final List<Allele> truthAlleles, final List<Allele> evalAlleles, final boolean expected) {
        final VariantContext truth = new VariantContextBuilder("truth", "chr1", 100, 100, truthAlleles).make();
        final VariantContext eval = new VariantContextBuilder("eval", "chr1", 100, 100, evalAlleles).make();

        Assert.assertEquals(new MergeMutect2CallsWithMC3().shouldVariantsBeMatched(truth, eval), expected);
    }
}
