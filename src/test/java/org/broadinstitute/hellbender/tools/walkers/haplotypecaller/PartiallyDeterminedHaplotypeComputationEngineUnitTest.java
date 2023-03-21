package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
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

    VariantContext SNP_C_90 = new VariantContextBuilder("a","20",90, 90, Arrays.asList(Allele.REF_A,Allele.ALT_C)).make();
    VariantContext DEL_AAAAAAA_98 = new VariantContextBuilder("a","20",98, 104, Arrays.asList(Allele.create("AAAAAAA", true),Allele.ALT_A)).make();
    VariantContext SNP_C_100 = new VariantContextBuilder("a","20",100, 100, Arrays.asList(Allele.REF_A,Allele.ALT_C)).make();
    VariantContext SNP_G_101 = new VariantContextBuilder("a","20",101, 101, Arrays.asList(Allele.REF_A,Allele.ALT_G)).make();
    VariantContext SNP_G_102 = new VariantContextBuilder("a","20",102, 102, Arrays.asList(Allele.REF_A,Allele.ALT_G)).make();
    VariantContext SNP_C_104 = new VariantContextBuilder("a","20",104, 104, Arrays.asList(Allele.REF_A,Allele.ALT_C)).make();
    VariantContext SNP_C_105 = new VariantContextBuilder("a","20",105, 105, Arrays.asList(Allele.REF_A,Allele.ALT_C)).make();
    VariantContext SNP_G_105 = new VariantContextBuilder("a","20",105, 105, Arrays.asList(Allele.REF_A,Allele.ALT_G)).make();
    VariantContext SNP_C_106 = new VariantContextBuilder("a","20",106, 106, Arrays.asList(Allele.REF_A,Allele.ALT_C)).make();
    VariantContext SNP_T_106 = new VariantContextBuilder("a","20",106, 106, Arrays.asList(Allele.REF_A,Allele.ALT_T)).make();
    VariantContext SNP_C_109 = new VariantContextBuilder("a","20",109, 109, Arrays.asList(Allele.REF_A,Allele.ALT_C)).make();
    VariantContext SNP_C_107 = new VariantContextBuilder("a","20",107, 107, Arrays.asList(Allele.REF_A,Allele.ALT_C)).make();

    VariantContext DEL_AA_105 = new VariantContextBuilder("a","20",105, 106, Arrays.asList(Allele.create("AA", true),Allele.ALT_A)).make();
    VariantContext DEL_AA_100 = new VariantContextBuilder("a","20",100, 101, Arrays.asList(Allele.create("AA", true),Allele.ALT_A)).make();
    VariantContext DEL_AAA_102 = new VariantContextBuilder("a","20",102, 104, Arrays.asList(Allele.create("AAA", true),Allele.ALT_A)).make();
    VariantContext DEL_AAAAAAA_102 = new VariantContextBuilder("a","20",102, 108, Arrays.asList(Allele.create("AAAAAAA", true),Allele.ALT_A)).make();


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
                { Collections.emptyList(), "AAAAAAAAAA", "10M", 0 },
                { Arrays.asList(SNP_C_105), "AAAAACAAAA", "5M1X4M", 0 },
                { Arrays.asList(SNP_C_100), "CAAAAAAAAA", "1X9M", 0 },
                { Arrays.asList(SNP_C_109), "AAAAAAAAAC", "9M1X", 0 },
                { Arrays.asList(SNP_C_105, SNP_C_106), "AAAAACCAAA", "5M2X3M", 0 },

                { Arrays.asList(DEL_AA_105), "AAAAAAAAA", "6M1D3M", 0 },
                { Arrays.asList(DEL_AA_100), "AAAAAAAAA", "1M1D8M", 0 },
                { Arrays.asList(DEL_AA_105, SNP_C_109), "AAAAAAAAC", "6M1D2M1X", 0 },
                { Arrays.asList(DEL_AA_105, SNP_C_107, SNP_C_109), "AAAAAACAC", "6M1D1X1M1X", 0 },

                { Arrays.asList(INS_TT_105),  "AAAAAATAAAA", "6M1I4M", 0 },
                { Arrays.asList(INS_GGG_106), "AAAAAAAGGAAA", "7M2I3M", 0 },
                { Arrays.asList(DEL_AA_100, INS_GGG_106, SNP_C_109), "AAAAAAGGAAC", "1M1D5M2I2M1X", 0 },

                //this tests that SNPS can be inserted immediately prior to (and following) indels
                { Arrays.asList( SNP_C_105, DEL_AA_105 ), "AAAAACAAA", "5M1X1D3M", 1 },
                { Arrays.asList( SNP_C_100, DEL_AA_100 ), "CAAAAAAAA", "1X1D8M", 1 },
                { Arrays.asList( SNP_C_100, DEL_AA_100, SNP_G_102 ), "CGAAAAAAA", "1X1D1X7M", 1 },
                { Arrays.asList( SNP_C_105, INS_TT_105 ), "AAAAACTAAAA", "5M1X1I4M", 1 },
                { Arrays.asList( SNP_C_100, INS_TT_100, SNP_G_101 ), "CTGAAAAAAAA", "1X1I1X8M", 1 },
                { Arrays.asList( SNP_C_100, INS_TT_100, SNP_G_101, SNP_C_105, DEL_AA_105 ), "CTGAAACAAA", "1X1I1X3M1X1D3M", 2 },

                //testing that the logic around anchor bases isn't resulting in variants being dropped accidntally
                { Arrays.asList( SNP_C_104, DEL_AA_105 ), "AAAACAAAA", "4M1X1M1D3M", 0 },
                { Arrays.asList( SNP_C_104, INS_TT_105 ), "AAAACATAAAA", "4M1X1M1I4M", 0 },

        };
    }
    @Test(dataProvider = "testConstructHaplotypeFromVariantsDataProvider")
    public void basicConstructHaplotypeFromVariants(List<VariantContext> variants, String expectedBases, String expectedCigar, int numberOfCompounds) {
        Haplotype ref = new Haplotype("AAAAAAAAAA".getBytes(), true, 500, TextCigarCodec.decode("10M"));
        ref.setGenomeLocation(new SimpleInterval("20", 100, 110));

        Haplotype result = PartiallyDeterminedHaplotypeComputationEngine.constructHaplotypeFromVariants(ref, variants, true);
        Assert.assertEquals(result.getBases(), expectedBases.getBytes());
        Assert.assertEquals(result.getCigar(), TextCigarCodec.decode(expectedCigar));

        // Assert that the resulting event map matches the input variants:
        EventMap resultEMap = result.getEventMap();
        // NOTE, because of representation in VCF lines, the compound alleles get compressed into a single in the event map, here we assert that this is correct.
        Assert.assertEquals(resultEMap.getNumberOfEvents(), variants.size() - numberOfCompounds);
    }

    @Test(expectedExceptions = GATKException.class)
    public void TestOutOfOrderInputs() {
        Haplotype ref = new Haplotype("AAAAAAAAAA".getBytes(), true, 500, TextCigarCodec.decode("10M"));
        ref.setGenomeLocation(new SimpleInterval("20", 100, 110));
        List<VariantContext> variants = Arrays.asList(SNP_C_105, SNP_G_105);

        Haplotype result = PartiallyDeterminedHaplotypeComputationEngine.constructHaplotypeFromVariants(ref, variants, true);
    }

    @Test(expectedExceptions = GATKException.class)
    public void TestSNPsOverlapping() {
        Haplotype ref = new Haplotype("AAAAAAAAAA".getBytes(), true, 500, TextCigarCodec.decode("10M"));
        ref.setGenomeLocation(new SimpleInterval("20", 100, 110));
        List<VariantContext> variants = Arrays.asList(SNP_C_109, DEL_AA_100);

        Haplotype result = PartiallyDeterminedHaplotypeComputationEngine.constructHaplotypeFromVariants(ref, variants, true);
    }

    @Test(expectedExceptions = GATKException.class)
    public void TestVariantNotOverlappingHap() {
        Haplotype ref = new Haplotype("AAAAAAAAAA".getBytes(), true, 500, TextCigarCodec.decode("10M"));
        ref.setGenomeLocation(new SimpleInterval("20", 100, 110));
        List<VariantContext> variants = Arrays.asList(SNP_C_90);

        Haplotype result = PartiallyDeterminedHaplotypeComputationEngine.constructHaplotypeFromVariants(ref, variants, true);
    }

    @Test(expectedExceptions = GATKException.class)
    public void TestVariantIndelPartiallyOverlapping() {
        Haplotype ref = new Haplotype("AAAAAAAAAA".getBytes(), true, 500, TextCigarCodec.decode("10M"));
        ref.setGenomeLocation(new SimpleInterval("20", 100, 110));
        List<VariantContext> variants = Arrays.asList(DEL_AAAAAAA_98);

        Haplotype result = PartiallyDeterminedHaplotypeComputationEngine.constructHaplotypeFromVariants(ref, variants, true);
    }

    //This is a test asserting that a real edge case that was prone to cause failures in the PDHMM is handled properly when compound variants are taken into account.
    //(62,Rlen=1,[C])->(82,Rlen=1,[C])->(84,Rlen=13,[C])
    @Test
    public void testMessyAlignemntSite() {
        Haplotype ref = new Haplotype("AAGAAAGATGGAGGCCCAGCCAGATCTGGACCCCACAGGCCGTCTCCCCACACAGCCATTCATGTGGTCTACTTCCAGCCATTCATGTGGTCTATTTCCAAGAAAATAGCCCATCCCCCCAAGATAACACCTTCTCAAAAACTTTACAGCTTTGTGTCTACACTGATATTTAGGTATTTTCTTTCTTTTTTTTTTATGATTAACACATCTAATTCAAGAATATCTTGGCAGGATATTCCCCGCTTAGGAAATG".getBytes(), true, 575, TextCigarCodec.decode("253M"));
        ref.setGenomeLocation(new SimpleInterval("20", 24152646, 24152898));

        VariantContext VC1 = new VariantContextBuilder("a", "20", 24152708, 24152708, Arrays.asList(Allele.REF_T, Allele.ALT_C)).make();
        VariantContext VC2 = new VariantContextBuilder("a", "20", 24152728, 24152728, Arrays.asList(Allele.REF_T, Allele.ALT_C)).make();
        VariantContext VC3 = new VariantContextBuilder("a", "20", 24152729, 24152741, Arrays.asList(Allele.create("CATGTGGTCTATT", true), Allele.ALT_C)).make();

        List<VariantContext> variants = Arrays.asList(VC1, VC2, VC3);

        Haplotype result = PartiallyDeterminedHaplotypeComputationEngine.constructHaplotypeFromVariants(ref, variants, true);
        Assert.assertEquals(result.getCigar(), TextCigarCodec.decode("62M1X19M1X1M12D157M"));

        // Assert that the resulting event map matches the input variants:
        EventMap resultEMap = result.getEventMap();
        Assert.assertEquals(resultEMap.getNumberOfEvents(), variants.size());
        for (VariantContext v : variants) {
            VariantContext actualVC = resultEMap.get(v.getStart());
            Assert.assertNotNull(actualVC);
            Assert.assertEquals(actualVC.getAlleles(), v.getAlleles());
        }
    }


    @DataProvider
    public Object[][] testGeneratePDHaplotypeDataProvider() {
        return new Object[][] {
                {Arrays.asList(SNP_C_105, SNP_C_106), SNP_C_106, false, "AAAAAACAAA", new byte[]{0,0,0,0,0,17,0,0,0,0}, "6M1X3M"},
                {Arrays.asList(SNP_C_105, SNP_C_106), SNP_C_106, true , "AAAAAAAAAA", new byte[]{0,0,0,0,0,17,0,0,0,0}, "10M"},

                {Arrays.asList(INS_TT_103, SNP_C_105, SNP_C_106), INS_TT_103, false, "AAAATAAAAAA", new byte[]{0,0,0,0,0,0,17,17,0,0,0}, "4M1I6M"},
                {Arrays.asList(INS_TT_103, SNP_C_105, SNP_C_106), INS_TT_103, true , "AAAAAAAAAA",  new byte[]{0,0,0,0,0,17,17,0,0,0}, "10M"},
                {Arrays.asList(INS_TT_103, SNP_C_105, SNP_C_106), SNP_C_105,  false, "AAAATACAAAA", new byte[]{0,0,0,0,6,0,0,17,0,0,0}, "4M1I1M1X4M"},
                {Arrays.asList(INS_TT_103, SNP_C_105, SNP_C_106), SNP_C_105,  true , "AAAATAAAAAA", new byte[]{0,0,0,0,6,0,0,17,0,0,0}, "4M1I6M"},

                {Arrays.asList(DEL_AAA_102, SNP_C_105, SNP_C_106), DEL_AAA_102, false, "AAAAAAAA"  , new byte[]{0,0,0,17,17,0,0,0}, "3M2D5M"},
                {Arrays.asList(DEL_AAA_102, SNP_C_105, SNP_C_106), DEL_AAA_102, true , "AAAAAAAAAA", new byte[]{0,0,0,0,0,17,17,0,0,0}, "10M"},
                {Arrays.asList(DEL_AAA_102, SNP_C_105, SNP_C_106), SNP_C_105,  false,  "AAAAACAAAA", new byte[]{0,0,0,2,4,0,17,0,0,0}, "5M1X4M"},
                {Arrays.asList(DEL_AAA_102, SNP_C_105, SNP_C_106), SNP_C_105,  true ,  "AAAAAAAAAA", new byte[]{0,0,0,2,4,0,17,0,0,0}, "10M"},
                {Arrays.asList(DEL_AAA_102, SNP_C_105, SNP_C_106), SNP_C_106,  false,  "AAAAAACAAA", new byte[]{0,0,0,2,4,17,0,0,0,0}, "6M1X3M"},
                {Arrays.asList(DEL_AAA_102, SNP_C_105, SNP_C_106), SNP_C_106,  true ,  "AAAAAAAAAA", new byte[]{0,0,0,2,4,17,0,0,0,0}, "10M"},

                // making sure we support "complex allels" from DRAGEN
                {Arrays.asList(DEL_AAA_102, SNP_C_105, SNP_C_106, INS_GGG_106), SNP_C_105,  false ,  "AAAAACAGGAAA", new byte[]{0,0,0,2,4,0,17,2,4,0,0,0}, "5M1X1M2I3M"},
                {Arrays.asList(DEL_AAA_102, SNP_C_105, SNP_C_106, SNP_T_106, INS_GGG_106), SNP_C_105,  true ,  "AAAAAAAGGAAA", new byte[]{0,0,0,2,4,0,81,2,4,0,0,0}, "7M2I3M"},
                {Arrays.asList(DEL_AAA_102, SNP_C_105, SNP_C_106, INS_GGG_106), DEL_AAA_102,  false ,  "AAAAAGGAAA", new byte[]{0,0,0,17,17,2,4,0,0,0}, "3M2D2M2I3M"},
                {Arrays.asList(DEL_AAA_102, SNP_C_105, SNP_C_106, SNP_T_106, INS_GGG_106), DEL_AAA_102,  true ,  "AAAAAAAGGAAA", new byte[]{0,0,0,0,0,17,81,2,4,0,0,0}, "7M2I3M"},
                {Arrays.asList(SNP_G_101, SNP_C_105, DEL_AA_105), SNP_G_101,  false ,  "AGAAAAAAAA", new byte[]{0,0,0,0,0,17,6,0,0,0}, "1M1X8M"},
                {Arrays.asList(SNP_G_101, SNP_C_105, DEL_AA_105), SNP_G_101,  true ,   "AAAAAAAAAA", new byte[]{0,0,0,0,0,17,6,0,0,0}, "10M"},

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
        Assert.assertEquals(result.getDeterminedPosition(), targetVariant.getStart());
    }

    // NOTE: This is an enfocement of a behavior that I consider to be a bug in DRAGEN. Specifically my assumption that we needn't ever concern
    // ourselves with overlapping variants turns out to be false... As it turns out in DRAGEN, they are entirely accepting of construcitng a
    // PD haplotype that is REF at bases that underly a spanning deletion... This means (for example) that if we have a 10 base undetermined
    // deletion from 100-109 and we have a determined ref deletion at position 105-106, that we should STILL construct the halplotype with
    // PD bases from 100-109 even though it means we are assigning that deletion at position 100 to be ref (essentially enforcing that we
    // don't handle spanning deletions). Joint Deteciotion will likely override this behavior in the future.
    @Test
    public void testDeletionUnderlapingDeterminedBases() {
        Haplotype ref = new Haplotype("AAAAAAAAAA".getBytes(), true, 500, TextCigarCodec.decode("10M"));
        ref.setGenomeLocation(new SimpleInterval("20", 100, 110));

        PartiallyDeterminedHaplotype result = PartiallyDeterminedHaplotypeComputationEngine.createNewPDHaplotypeFromEvents(ref, DEL_AA_105, true, Arrays.asList(DEL_AAAAAAA_102, DEL_AA_105));
        Assert.assertEquals(new String(result.getBases()), "AAAAAAAAAA");
        Assert.assertEquals(result.getAlternateBases(), new byte[]{0,0,0,2,0,0,0,0,4,0});
        Assert.assertEquals(result.getCigar(), TextCigarCodec.decode("10M"));
        Assert.assertEquals(result.getDeterminedPosition(), DEL_AA_105.getStart());
    }
}