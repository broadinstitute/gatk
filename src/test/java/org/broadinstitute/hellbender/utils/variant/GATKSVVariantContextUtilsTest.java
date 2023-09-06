package org.broadinstitute.hellbender.utils.variant;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class GATKSVVariantContextUtilsTest {

    @DataProvider(name = "getSymbolicAlleleSymbolsTestData")
    public Object[][] getSymbolicAlleleSymbolsTestData() {
        return new Object[][]{
                {Allele.create("<DUP>"), new String[] {"DUP"}},
                {Allele.create("<DUP:TANDEM>"), new String[] {"DUP", "TANDEM"}},
                {Allele.create("<INS:MEI:LINE>"), new String[] {"INS", "MEI", "LINE"}}
        };
    }

    @Test(dataProvider= "getSymbolicAlleleSymbolsTestData")
    public void getSymbolicAlleleSymbolsTest(final Allele allele, final String[] result) {
        final String[] test = GATKSVVariantContextUtils.getSymbolicAlleleSymbols(allele);
        Assert.assertEquals(test, result);
    }

    @Test
    public void isCnvTypeTest() {
        Assert.assertTrue(GATKSVVariantContextUtils.isCnvType(StructuralVariantType.CNV));
        Assert.assertTrue(GATKSVVariantContextUtils.isCnvType(StructuralVariantType.DEL));
        Assert.assertTrue(GATKSVVariantContextUtils.isCnvType(StructuralVariantType.DUP));
        Assert.assertFalse(GATKSVVariantContextUtils.isCnvType(StructuralVariantType.INS));
        Assert.assertFalse(GATKSVVariantContextUtils.isCnvType(StructuralVariantType.INV));
        Assert.assertFalse(GATKSVVariantContextUtils.isCnvType(StructuralVariantType.BND));
    }
}