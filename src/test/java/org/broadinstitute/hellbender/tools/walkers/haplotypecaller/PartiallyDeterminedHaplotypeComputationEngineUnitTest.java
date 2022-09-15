package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.haplotype.PartiallyDeterminedHaplotype;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static org.testng.Assert.*;

public class PartiallyDeterminedHaplotypeComputationEngineUnitTest extends GATKBaseTest {

    VariantContext SNP_C_100 = new VariantContextBuilder("a","20",100, 100, Arrays.asList(Allele.REF_A,Allele.ALT_C)).make();
    VariantContext SNP_G_101 = new VariantContextBuilder("a","20",101, 101, Arrays.asList(Allele.REF_A,Allele.ALT_G)).make();
    VariantContext SNP_G_102 = new VariantContextBuilder("a","20",102, 102, Arrays.asList(Allele.REF_A,Allele.ALT_G)).make();
    VariantContext SNP_C_105 = new VariantContextBuilder("a","20",105, 105, Arrays.asList(Allele.REF_A,Allele.ALT_C)).make();
    VariantContext SNP_C_106 = new VariantContextBuilder("a","20",106, 106, Arrays.asList(Allele.REF_A,Allele.ALT_C)).make();
    VariantContext SNP_C_109 = new VariantContextBuilder("a","20",109, 109, Arrays.asList(Allele.REF_A,Allele.ALT_C)).make();
    VariantContext SNP_C_107 = new VariantContextBuilder("a","20",107, 107, Arrays.asList(Allele.REF_A,Allele.ALT_C)).make();

    VariantContext DEL_AA_105 = new VariantContextBuilder("a","20",105, 106, Arrays.asList(Allele.create("AA", true),Allele.ALT_A)).make();
    VariantContext DEL_AA_100 = new VariantContextBuilder("a","20",100, 101, Arrays.asList(Allele.create("AA", true),Allele.ALT_A)).make();
    VariantContext DEL_AAA_102 = new VariantContextBuilder("a","20",102, 104, Arrays.asList(Allele.create("AAA", true),Allele.ALT_A)).make();

    VariantContext INS_TT_105 = new VariantContextBuilder("a","20",105, 105, Arrays.asList(Allele.REF_A, Allele.create("AT"))).make();
    VariantContext INS_TT_103 = new VariantContextBuilder("a","20",103, 103, Arrays.asList(Allele.REF_A, Allele.create("AT"))).make();
    VariantContext INS_TT_100 = new VariantContextBuilder("a","20",100, 100, Arrays.asList(Allele.REF_A, Allele.create("AT"))).make();
    VariantContext INS_GGG_106 = new VariantContextBuilder("a","20",106, 106, Arrays.asList(Allele.REF_A, Allele.create("AGG"))).make();

    // TODO THESE ARE FOR INVALID TEST CASES
    VariantContext SNP_C_99 = new VariantContextBuilder("a","20",99, 99, Arrays.asList(Allele.REF_A,Allele.ALT_C)).make();
    VariantContext SNP_C_120 = new VariantContextBuilder("a","20",120, 120, Arrays.asList(Allele.REF_A,Allele.ALT_C)).make();



    @DataProvider
    public Object[][] testConstructHaplotypeFromVariantsDataProvider() {
        return new Object[][] {
                { Collections.emptyList(), "AAAAAAAAAA", "10M" },
                { Arrays.asList(SNP_C_105), "AAAAACAAAA", "5M1X4M" },
                { Arrays.asList(SNP_C_100), "CAAAAAAAAA", "1X9M" },
                { Arrays.asList(SNP_C_109), "AAAAAAAAAC", "9M1X" },
                { Arrays.asList(SNP_C_105, SNP_C_106), "AAAAACCAAA", "5M2X3M" },

                { Arrays.asList(DEL_AA_105), "AAAAAAAAA", "6M1D3M" },
                { Arrays.asList(DEL_AA_100), "AAAAAAAAA", "1M1D8M" },
                { Arrays.asList(DEL_AA_105, SNP_C_109), "AAAAAAAAC", "6M1D2M1X" },
                { Arrays.asList(DEL_AA_105, SNP_C_107, SNP_C_109), "AAAAAACAC", "6M1D1X1M1X" },

                { Arrays.asList(INS_TT_105),  "AAAAAATAAAA", "6M1I4M" },
                { Arrays.asList(INS_GGG_106), "AAAAAAAGGAAA", "7M2I3M" },
                { Arrays.asList(DEL_AA_100, INS_GGG_106, SNP_C_109), "AAAAAAGGAAC", "1M1D5M2I2M1X" },

                //this tests that SNPS can be inserted immediately prior to (and following) indels
                { Arrays.asList( SNP_C_105, DEL_AA_105 ), "AAAAACAAA", "5M1X1D3M" },
                { Arrays.asList( SNP_C_100, DEL_AA_100 ), "CAAAAAAAA", "1X1D8M" },
                { Arrays.asList( SNP_C_100, DEL_AA_100, SNP_G_102 ), "CGAAAAAAA", "1X1D1X7M" },
                { Arrays.asList( SNP_C_105, INS_TT_105 ), "AAAAACTAAAA", "5M1X1I4M" },
                { Arrays.asList( SNP_C_100, INS_TT_100, SNP_G_101 ), "CTGAAAAAAAA", "1X1I1X8M" },
                { Arrays.asList( SNP_C_100, INS_TT_100, SNP_G_101, SNP_C_105, DEL_AA_105 ), "CTGAAACAAA", "1X1I1X3M1X1D3M" },

        };
    }
    @Test(dataProvider = "testConstructHaplotypeFromVariantsDataProvider")
    public void basicConstructHaplotypeFromVariants(List<VariantContext> variants, String expectedBases, String expectedCigar) {
        Haplotype ref = new Haplotype("AAAAAAAAAA".getBytes(), true, 500, TextCigarCodec.decode("10M"));
        ref.setGenomeLocation(new SimpleInterval("20", 100, 110));

        Haplotype result = PartiallyDeterminedHaplotypeComputationEngine.constructHaplotypeFromVariants(ref, variants);
        Assert.assertEquals(result.getBases(), expectedBases.getBytes());
        Assert.assertEquals(result.getCigar(), TextCigarCodec.decode(expectedCigar));
    }
    //TODO TESTS TO MAKE:
    // ASSSERT IT FAILS IF STARTS BEFORE OR AFTER
    // ASSERT it gives a good exception if the order is incorrect

    //This test is here for the sake of them being related operations

    @DataProvider
    public Object[][] testGeneratePDHaplotypeDataProvider() {
        return new Object[][] {
                {Arrays.asList(SNP_C_105, SNP_C_106), SNP_C_106, false, "AAAAAACAAA", new byte[]{0,0,0,0,0,17,0,0,0,0}, "10M"},
                {Arrays.asList(SNP_C_105, SNP_C_106), SNP_C_106, true , "AAAAAAAAAA", new byte[]{0,0,0,0,0,17,0,0,0,0}, "10M"},

                {Arrays.asList(INS_TT_103, SNP_C_105, SNP_C_106), INS_TT_103, false, "AAAATAAAAAA", new byte[]{0,0,0,0,0,0,17,17,0,0,0}, "4M1I6M"},
                {Arrays.asList(INS_TT_103, SNP_C_105, SNP_C_106), INS_TT_103, true , "AAAAAAAAAA",  new byte[]{0,0,0,0,0,17,17,0,0,0}, "10M"},
                {Arrays.asList(INS_TT_103, SNP_C_105, SNP_C_106), SNP_C_105,  false, "AAAATACAAAA", new byte[]{0,0,0,0,6,0,0,17,0,0,0}, "4M1I6M"},
                {Arrays.asList(INS_TT_103, SNP_C_105, SNP_C_106), SNP_C_105,  true , "AAAATAAAAAA", new byte[]{0,0,0,0,6,0,0,17,0,0,0}, "4M1I6M"},

                {Arrays.asList(DEL_AAA_102, SNP_C_105, SNP_C_106), DEL_AAA_102, false, "AAAAAAAA"  , new byte[]{0,0,0,17,17,0,0,0}, "8M"},
                {Arrays.asList(DEL_AAA_102, SNP_C_105, SNP_C_106), DEL_AAA_102, true , "AAAAAAAAAA", new byte[]{0,0,0,0,0,17,17,0,0,0}, "10M"},
                {Arrays.asList(DEL_AAA_102, SNP_C_105, SNP_C_106), SNP_C_105,  false,  "AAAAACAAAA", new byte[]{0,0,0,2,4,0,17,0,0,0}, "10M"},
                {Arrays.asList(DEL_AAA_102, SNP_C_105, SNP_C_106), SNP_C_105,  true ,  "AAAAAAAAAA", new byte[]{0,0,0,2,4,0,17,0,0,0}, "10M"},
                {Arrays.asList(DEL_AAA_102, SNP_C_105, SNP_C_106), SNP_C_106,  false,  "AAAAAACAAA", new byte[]{0,0,0,2,4,17,0,0,0,0}, "10M"},
                {Arrays.asList(DEL_AAA_102, SNP_C_105, SNP_C_106), SNP_C_106,  true ,  "AAAAAAAAAA", new byte[]{0,0,0,2,4,17,0,0,0,0}, "10M"},

                // making sure we support "complex allels" from DRAGEN
                {Arrays.asList(DEL_AAA_102, SNP_C_105, SNP_C_106, INS_GGG_106), SNP_C_105,  false ,  "AAAAACAGGAAA", new byte[]{0,0,0,2,4,0,17,2,4,0,0,0}, "7M2I3M"},
        };
    }
    @Test(dataProvider = "testGeneratePDHaplotypeDataProvider")
    public void testGeneratePDHaplotypeFromVariants(List<VariantContext> variants, VariantContext targetVariant, boolean useRefBase, String expectedBases, byte[] expectedAltArray, String expectedCigar) {
        Haplotype ref = new Haplotype("AAAAAAAAAA".getBytes(), true, 500, TextCigarCodec.decode("10M"));
        ref.setGenomeLocation(new SimpleInterval("20", 100, 110));

        PartiallyDeterminedHaplotype result = PartiallyDeterminedHaplotypeComputationEngine.createNewPDHaplotypeFromEvents(ref, targetVariant, useRefBase, variants);
        Assert.assertEquals(new String(result.getBases()), expectedBases);
        Assert.assertEquals(result.getAlternateBases(), expectedAltArray);
        Assert.assertEquals(result.getCigar(), TextCigarCodec.decode(expectedCigar));
    }


    @Test
    public void testConstructArtificialHaplotypeAndTestEquivalentEvents() {
        Haplotype ref = new Haplotype("AAAAAAAAAA".getBytes(), true, 500, TextCigarCodec.decode("10M"));
        //TODO better tests of this messy case

    }
}