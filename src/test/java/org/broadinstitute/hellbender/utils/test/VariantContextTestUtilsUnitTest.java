package org.broadinstitute.hellbender.utils.test;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

public class VariantContextTestUtilsUnitTest {

    private static final Allele ARef = Allele.create("A", true);
    private static final Allele T = Allele.create("T");
    private static final Allele C = Allele.create("C");
    private static final Allele G = Allele.create("G");

    @DataProvider(name="scientificNotationValuesToNormalize")
    public Object[][] getScientificNotationValuesToNormalize(){
        final Object aSpecificObject = new Object();
        return new Object[][] {
                {"-2.172e+00", -2.172},
                {"-2.172", -2.172},
                {"-2.172e+01", -21.72},
                {-21.72, -21.72},
                {10, 10},
                {"SomeValue", "SomeValue"},
                {aSpecificObject,  aSpecificObject}

        };
    }

    @Test(dataProvider = "scientificNotationValuesToNormalize")
    public void testNormalizeScientificNotation(Object toNormalize, Object expected){
        Assert.assertEquals(VariantContextTestUtils.normalizeScientificNotation(toNormalize), expected);
    }

    @DataProvider(name="integerValuesToNormalize")
    public Object[][] getIntegerValuesToNormalize(){
        final Object aSpecificObject = new Object();
        return new Object[][] {
                {"27", new Integer(27)},
                {"-27", new Integer(-27)},
                {"0", new Integer(0)},
                {"-27014", new Integer(-27014)},
                {1, 1},
                {1, 1},
                {-1, -1},
                {1, 1},
                {"-2.172e+00", "-2.172e+00"},
                {"-2.172", "-2.172"},
                {-21.72, -21.72},
                {10, 10},
                {"SomeValue", "SomeValue"},
                {aSpecificObject,  aSpecificObject}
        };
    }

    @Test(dataProvider = "integerValuesToNormalize")
    public void testNormalizeInteger(Object toNormalize, Object expected){
        Assert.assertEquals(VariantContextTestUtils.normalizeToInteger(toNormalize), expected);
    }

    @DataProvider
    public Object[][] getEqualVariantContexts(){

        return new Object[][] {
                {GenotypeBuilder.create("sample", Arrays.asList(T, G)), GenotypeBuilder.create("sample", Arrays.asList(G,T))}
        };
    }


    @Test(dataProvider = "getEqualVariantContexts")
    public void testGenotypeEquality(Genotype left, Genotype right){
        VariantContextTestUtils.assertGenotypesAreEqual(left,right);
    }


    @DataProvider
    public Object[][] getTListsToRemap(){
        return new Object[][]{
                {Arrays.asList(1, 2, 3, 4), Arrays.asList(1, 3, 4, 2), Arrays.asList(0, 2, 3, 1)},
                {Arrays.asList(1, 2, 3, 4), Arrays.asList(1, 2, 3, 4), Arrays.asList(0, 1, 2, 3)},
                {Arrays.asList("a", "b", "c","d"), Arrays.asList("a", "d", "c", "b"), Arrays.asList(0, 3,2,1)}
        };
    }


    @Test(dataProvider = "getTListsToRemap")
    public void testTRemapping(List<Integer> original, List<Integer> expectedRemappedValues, List<Integer> mapping){
        Assert.assertEquals(VariantContextTestUtils.remapRTypeValues(original, mapping), expectedRemappedValues);
    }

    @DataProvider
    public Object[][] getAListsToRemap(){
        return new Object[][]{
                {Arrays.asList(1, 2, 3), Arrays.asList(2, 3, 1), Arrays.asList(0, 2, 3, 1)},
                {Arrays.asList(1, 2, 3), Arrays.asList(1, 2, 3), Arrays.asList(0, 1, 2, 3)},
                {Arrays.asList("a", "b", "c"), Arrays.asList("c", "b" , "a"), Arrays.asList(0, 3, 2, 1)}
        };
    }

    @Test(dataProvider = "getAListsToRemap")
    public void testARemapping(List<Integer> original, List<Integer> expectedRemappedValues, List<Integer> mapping){
        Assert.assertEquals(VariantContextTestUtils.remapATypeValues(original, mapping), expectedRemappedValues);
    }


    @DataProvider
    public Object[][] getGListsToRemap() {
        return new Object[][]{
            {Arrays.asList(0,1,2), Arrays.asList(2,1,0), Arrays.asList(ARef, T), Arrays.asList(T, ARef), 2},
            {Arrays.asList(0,1,2,3,4,5), Arrays.asList(2,1,0,4,3,5), Arrays.asList(ARef, T, C), Arrays.asList(T, ARef, C), 2}
        };
    }

    @Test(dataProvider = "getGListsToRemap")
    public void testGListsToRemap(List<Object> original, List<Object> expected, List<Allele> originalAlleles, List<Allele> remappedAlleles, int ploidy){
        Assert.assertEquals(VariantContextTestUtils.remapGTypeValues(original, ploidy,
                                                                     VariantContextTestUtils.createAlleleIndexMap(
                                                                             originalAlleles, remappedAlleles)), expected);
    }

    @DataProvider
    public Object[][] getGenotypes(){
        final HashSet<VCFHeaderLine> vcfHeaderLines = new HashSet<>();
        VCFStandardHeaderLines.addStandardFormatLines(vcfHeaderLines, true, Arrays.asList(VCFConstants.DEPTH_KEY,
                                                                                          VCFConstants.GENOTYPE_ALLELE_DEPTHS,
                                                                                          VCFConstants.GENOTYPE_PL_KEY,
                                                                                          VCFConstants.GENOTYPE_QUALITY_KEY));
        vcfHeaderLines.add(new VCFFormatHeaderLine("ZZ", VCFHeaderLineCount.A, VCFHeaderLineType.Float, "A test tag"));
        final VCFHeader header = new VCFHeader(vcfHeaderLines);
        final String sample = "sample";
        final Genotype unorderedAlleles = new GenotypeBuilder(sample, Arrays.asList(T, C)).make();
        final Genotype orderedAlleles = new GenotypeBuilder(sample, Arrays.asList(C, T)).make();
        final List<Integer> reordering = Arrays.asList(0, 2, 1);
        return new Object[][]{
                {unorderedAlleles, reordering, header, orderedAlleles},
                {new GenotypeBuilder(unorderedAlleles).DP(10).make(), reordering, header,  new GenotypeBuilder(orderedAlleles).DP(10).make()},
                {new GenotypeBuilder(unorderedAlleles).AD(new int[]{5, 3, 1}).make(), reordering, header, new GenotypeBuilder(orderedAlleles).AD(new int[]{5, 1, 3}).make()},
                {new GenotypeBuilder(unorderedAlleles).attribute("ZZ", new double[]{1.0, 2.0}).make(), reordering, header, new GenotypeBuilder(orderedAlleles).attribute("ZZ", new double[]{2.0, 1.0}).make()}
        };
    }

    @Test(dataProvider = "getGenotypes")
    public void testRearrangeGenotypes(Genotype gt, List<Integer> alleleIndexMap, VCFHeader header, Genotype expected){
        VariantContextTestUtils.assertGenotypesAreEqual(VariantContextTestUtils.reorderGenotypeAlleles(header, alleleIndexMap, gt), expected);
    }


    @DataProvider
    public static Object[][] getVariantContexts() {
        final HashSet<VCFHeaderLine> vcfHeaderLines = new HashSet<>();
        VCFStandardHeaderLines.addStandardFormatLines(vcfHeaderLines, true, Arrays.asList(VCFConstants.DEPTH_KEY,
                                                                                          VCFConstants.GENOTYPE_ALLELE_DEPTHS,
                                                                                          VCFConstants.GENOTYPE_PL_KEY,
                                                                                          VCFConstants.GENOTYPE_QUALITY_KEY));
        vcfHeaderLines.add(new VCFFormatHeaderLine("ZZ", VCFHeaderLineCount.A, VCFHeaderLineType.Float, "A test tag"));
        final VCFHeader header = new VCFHeader(vcfHeaderLines);
        final String sample = "sample";
        final Genotype unorderedAlleles = new GenotypeBuilder(sample, Arrays.asList(T, C)).make();
        final Genotype orderedAlleles = new GenotypeBuilder(sample, Arrays.asList(C, T)).make();
        final List<Integer> reordering = Arrays.asList(0, 2, 1);
        final VariantContext baseVC = new VariantContextBuilder("test", "1", 100, 100, Arrays.asList(ARef, T, C)).make();
        return new Object[][]{
                {baseVC, header, new VariantContextBuilder(baseVC).alleles(Arrays.asList(ARef, T, C)).make()},
        };
    }


    @Test(dataProvider = "getVariantContexts")
    public void testRearrangeVariantContexts(VariantContext vc, VCFHeader header, VariantContext expected){
        VariantContextTestUtils.assertVariantContextsAreEqual(VariantContextTestUtils.sortAlleles(vc, header), expected);
    }
}


