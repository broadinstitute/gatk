package org.broadinstitute.hellbender.testutils;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class VariantContextTestUtilsUnitTest extends GATKBaseTest {

    private static final Allele ARef = Allele.create("A", true);
    private static final Allele T = Allele.create("T");
    private static final Allele C = Allele.create("C");
    private static final Allele CC = Allele.create("CC");
    private static final Allele G = Allele.create("G");
    private static final Allele GG = Allele.create("GG");
    private static final String SAMPLE_1 = "NA1";
    private static final String SAMPLE_2 = "NA2";

    @DataProvider
    public Object[][] testReadEntireVCFIntoMemoryData() {
        return new Object[][] {
                { dbsnp_138_b37_20_21_vcf, 160, 9594 }
        };
    }

    @Test(dataProvider = "testReadEntireVCFIntoMemoryData")
    public void testReadEntireVCFIntoMemory(final String vcf, final int expectedNumHeaderLines, final int expectedNumRecords) {
        testReadEntireVCFIntoMemoryHelper(vcf, expectedNumHeaderLines, expectedNumRecords);
    }

    @DataProvider
    public Object[][] testReadEntireVCFIntoMemoryFromGCSData() {
        return new Object[][] {
                { "large/dbsnp_138.b37.20.21.vcf", 160, 9594 }
        };
    }

    @Test(groups = {"bucket"}, dataProvider = "testReadEntireVCFIntoMemoryFromGCSData")
    public void testReadEntireVCFIntoMemoryFromGCS(final String vcf, final int expectedNumHeaderLines, final int expectedNumRecords) {
        testReadEntireVCFIntoMemoryHelper(getGCPTestInputPath() + vcf, expectedNumHeaderLines, expectedNumRecords);
    }

    private void testReadEntireVCFIntoMemoryHelper(final String vcf, final int expectedNumHeaderLines, final int expectedNumRecords) {
        final Pair<VCFHeader, List<VariantContext>> result = VariantContextTestUtils.readEntireVCFIntoMemory(vcf);

        Assert.assertEquals(result.getLeft().getMetaDataInInputOrder().size(), expectedNumHeaderLines, "Wrong number of header lines in VCFHeader returned from VariantContextTestUtils.readEntireVCFIntoMemory()");
        Assert.assertEquals(result.getRight().size(), expectedNumRecords, "Wrong number of VariantContext records returned from VariantContextTestUtils.readEntireVCFIntoMemory()");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testReadEntireVCFIntoMemoryNonVCFInput() {
        VariantContextTestUtils.readEntireVCFIntoMemory(publicTestDir + "org/broadinstitute/hellbender/engine/example_features.bed");
    }

    @Test(expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testReadEntireVCFIntoMemoryNonExistentFile() {
        VariantContextTestUtils.readEntireVCFIntoMemory(getSafeNonExistentFile("testReadEntireVCFIntoMemoryNonExistentFile.vcf").getAbsolutePath());
    }

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
    public Object[][] getTListsToRemap(){
        return new Object[][]{
                {Arrays.asList(1, 2, 3, 4), Arrays.asList(1, 3, 4, 2), Arrays.asList(0, 2, 3, 1)},
                {Arrays.asList(1, 2, 3, 4), Arrays.asList(1, 2, 3, 4), Arrays.asList(0, 1, 2, 3)},
                {Arrays.asList("a", "b", "c","d"), Arrays.asList("a", "d", "c", "b"), Arrays.asList(0, 3,2,1)}
        };
    }


    @Test(dataProvider = "getTListsToRemap")
    public void testRemapping(List<Integer> original, List<Integer> expectedRemappedValues, List<Integer> mapping){
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
        Assert.assertEquals(VariantContextTestUtils.remapGTypeValues(original, originalAlleles, ploidy, remappedAlleles), expected);
    }

    @DataProvider
    public Object[][] getASTListsToRemap(){
        return new Object[][]{
                {"1|2|3|4", "1|3|4|2", Arrays.asList(0, 2, 3, 1)},
                {"|1, 2|3, 4", "|3, 4|1, 2", Arrays.asList(0, 2, 1)},
        };
    }

    @Test(dataProvider = "getTListsToRemap")
    public void testASRemapping(List<Integer> original, List<Integer> expectedRemappedValues, List<Integer> mapping){
        Assert.assertEquals(VariantContextTestUtils.remapRTypeValues(original, mapping), expectedRemappedValues);
    }

    @DataProvider
    public Object[][] alleleRemapExamples(){
        final double[] genotypeLikelihoods1 = {30,0,190};
        VariantContextBuilder builderA = new VariantContextBuilder("a","20",10433328,10433328,  Arrays.asList(ARef,G,C));
        VariantContextBuilder builderB = new VariantContextBuilder("a","20",10433328,10433328,  Arrays.asList(ARef,C,G));
        VariantContextBuilder builderC = new VariantContextBuilder("a","20",10433328,10433328,  Arrays.asList(ARef,C,T));
        VariantContextBuilder builderD = new VariantContextBuilder("a","20",10433328,10433328,  Arrays.asList(ARef,C,CC,G,GG));
        VariantContextBuilder builderE = new VariantContextBuilder("a","20",10433328,10433328,  Arrays.asList(ARef,GG,C,G,CC));
        VariantContext A = builderA.make();
        VariantContext B = builderB.make();
        VariantContext C = builderC.make();

        // TESTING ATTRIBUTES
        builderA.attribute("AS_RAW_MQ","40|20|10").attribute("RAW_MQ","20").attribute("PG", "10,20,30,40,50,60");
        builderB.attribute("AS_RAW_MQ","40|10|20").attribute("RAW_MQ","20").attribute("PG", "10,40,60,20,50,30");
        builderD.attribute("AS_RAW_MQ","0|1|2|3|4");
        builderE.attribute("AS_RAW_MQ","0|4|1|3|2");
        VariantContext A_RAWMQ = builderA.make();
        VariantContext B_RAWMQ = builderB.make();
        VariantContext D_RAWMQ = builderD.make();
        VariantContext E_RAWMQ = builderE.make();

        // TESTING GENOTYPES
        builderA.genotypes(Arrays.asList(new GenotypeBuilder(SAMPLE_1).alleles(Arrays.asList(ARef, G)).PL(new double[]{30,0,190,0,0,0}).GQ(30).AD(new int[]{1,2,3}).make(),
                                         new GenotypeBuilder(SAMPLE_2).alleles(Arrays.asList(VariantContextTestUtilsUnitTest.G, VariantContextTestUtilsUnitTest.C)).PL(new double[]{10,20,30,40,50,60}).GQ(30).AD(new int[]{3,4,5}).make()));
        builderB.genotypes(Arrays.asList(new GenotypeBuilder(SAMPLE_1).alleles(Arrays.asList(ARef, G)).PL(new double[]{30,0,0,0,0,190}).GQ(30).AD(new int[]{1,3,2}).make(),
                                         new GenotypeBuilder(SAMPLE_2).alleles(Arrays.asList(VariantContextTestUtilsUnitTest.G, VariantContextTestUtilsUnitTest.C)).PL(new double[]{10,40,60,20,50,30}).GQ(30).AD(new int[]{3,5,4}).make()));
        VariantContext A_withGenotypes = builderA.make();
        VariantContext B_withGenotypes = builderB.make();

        // NEGATIVE TESTS
        VariantContext B_RAWMQ_BAD = new VariantContextBuilder(B_RAWMQ).attribute("AS_RAW_MQ","40|10|10").make();
        VariantContext B_PG_BAD = new VariantContextBuilder(B_RAWMQ).attribute("PG","10,20,30,40,50,60").make();
        VariantContext B_withGenotypes_BADGQ = builderB.copy().genotypes(Arrays.asList(new GenotypeBuilder(SAMPLE_1).alleles(Arrays.asList(ARef, G)).PL(new double[]{30,0,0,0,0,190}).GQ(30).AD(new int[]{1,3,2}).make(),
                                                                                        new GenotypeBuilder(SAMPLE_2).alleles(Arrays.asList(VariantContextTestUtilsUnitTest.G, VariantContextTestUtilsUnitTest.C)).PL(new double[]{10,40,60,20,50,30}).AD(new int[]{3,5,4}).GQ(40).make())).make();
        VariantContext B_withGenotypes_BADPL = builderB.copy().genotypes(Arrays.asList(new GenotypeBuilder(SAMPLE_1).alleles(Arrays.asList(ARef, G)).PL(new double[]{30,0,190,0,0,0}).GQ(30).AD(new int[]{1,3,2}).make(),
                                                                                        new GenotypeBuilder(SAMPLE_2).alleles(Arrays.asList(VariantContextTestUtilsUnitTest.G, VariantContextTestUtilsUnitTest.C)).PL(new double[]{10,20,30,40,50,60}).GQ(30).AD(new int[]{3,5,4}).make())).make();
        VariantContext B_withGenotypes_BADAD = builderB.copy().genotypes(Arrays.asList(new GenotypeBuilder(SAMPLE_1).alleles(Arrays.asList(ARef, G)).PL(new double[]{30,0,0,0,0,190}).GQ(30).AD(new int[]{1,3,2}).make(),
                                                                                        new GenotypeBuilder(SAMPLE_2).alleles(Arrays.asList(VariantContextTestUtilsUnitTest.G, VariantContextTestUtilsUnitTest.C)).PL(new double[]{10,40,60,20,50,30}).GQ(30).AD(new int[]{3,4,4}).make())).make();


        return new Object[][]{{A,B,true},
                {A,C,false},
                {A_RAWMQ,B_RAWMQ,true},
                {D_RAWMQ,E_RAWMQ,true},
                {A_RAWMQ,B_RAWMQ_BAD,false},
                {A_RAWMQ,B_PG_BAD,false},
                {A_withGenotypes,B_withGenotypes,true},
                {A_withGenotypes,B_withGenotypes_BADGQ,false},
                {A_withGenotypes,B_withGenotypes_BADPL,false},
                {A_withGenotypes,B_withGenotypes_BADAD,false},};
    }

    @Test(dataProvider = "alleleRemapExamples")
    public void testOrderSortAlleles(VariantContext actual, VariantContext expected, boolean shouldSucceed) {
        VCFHeader header = VariantContextTestUtils.getCompleteHeader();

        if (shouldSucceed) {
            VariantContextTestUtils.assertVariantContextsAreEqualAlleleOrderIndependent(actual, expected, Collections.emptyList(), header);
        } else {
            Assert.assertThrows(AssertionError.class, () -> VariantContextTestUtils.assertVariantContextsAreEqualAlleleOrderIndependent(actual, expected, Collections.emptyList(), header));
        }
    }
}


