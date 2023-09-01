package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.haplotype.PartiallyDeterminedHaplotype;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.List;
import java.util.function.BiPredicate;

public class PartiallyDeterminedHaplotypeComputationEngineUnitTest extends GATKBaseTest {

    Event SNP_C_90 = new Event("20",90, Allele.REF_A,Allele.ALT_C);
    Event DEL_AAAAAAA_98 = new Event("20",98, Allele.create("AAAAAAA", true),Allele.ALT_A);
    Event SNP_C_100 = new Event("20",100, Allele.REF_A,Allele.ALT_C);
    Event SNP_G_101 = new Event("20",101, Allele.REF_A,Allele.ALT_G);
    Event SNP_G_102 = new Event("20",102, Allele.REF_A,Allele.ALT_G);
    Event SNP_C_104 = new Event("20",104, Allele.REF_A,Allele.ALT_C);
    Event SNP_C_105 = new Event("20",105, Allele.REF_A,Allele.ALT_C);
    Event SNP_G_105 = new Event("20",105, Allele.REF_A,Allele.ALT_G);
    Event SNP_C_106 = new Event("20",106, Allele.REF_A,Allele.ALT_C);
    Event SNP_T_106 = new Event("20",106, Allele.REF_A,Allele.ALT_T);
    Event SNP_C_109 = new Event("20",109, Allele.REF_A,Allele.ALT_C);
    Event SNP_C_107 = new Event("20",107, Allele.REF_A,Allele.ALT_C);

    Event DEL_AA_105 = new Event("20",105, Allele.create("AA", true),Allele.ALT_A);
    Event DEL_AA_100 = new Event("20",100, Allele.create("AA", true),Allele.ALT_A);
    Event DEL_AAA_102 = new Event("20",102, Allele.create("AAA", true),Allele.ALT_A);
    Event DEL_AAAAAAA_102 = new Event("20",102, Allele.create("AAAAAAA", true),Allele.ALT_A);


    Event INS_TT_105 = new Event("20",105, Allele.REF_A, Allele.create("AT"));
    Event INS_TT_103 = new Event("20",103, Allele.REF_A, Allele.create("AT"));
    Event INS_TT_100 = new Event("20",100, Allele.REF_A, Allele.create("AT"));
    Event INS_GGG_106 = new Event("20",106, Allele.REF_A, Allele.create("AGG"));

    // TODO THESE ARE FOR INVALID TEST CASES
    Event SNP_C_99 = new Event("20",99, Allele.REF_A,Allele.ALT_C);
    Event SNP_C_120 = new Event("20",120, Allele.REF_A,Allele.ALT_C);

    @DataProvider
    public Object[][] testConstructHaplotypeFromVariantsDataProvider() {
        return new Object[][] {
                { Collections.emptyList(), "AAAAAAAAAA", "10M", 0 },
                { List.of(SNP_C_105), "AAAAACAAAA", "5M1X4M", 0 },
                { List.of(SNP_C_100), "CAAAAAAAAA", "1X9M", 0 },
                { List.of(SNP_C_109), "AAAAAAAAAC", "9M1X", 0 },
                { List.of(SNP_C_105, SNP_C_106), "AAAAACCAAA", "5M2X3M", 0 },

                { List.of(DEL_AA_105), "AAAAAAAAA", "6M1D3M", 0 },
                { List.of(DEL_AA_100), "AAAAAAAAA", "1M1D8M", 0 },
                { List.of(DEL_AA_105, SNP_C_109), "AAAAAAAAC", "6M1D2M1X", 0 },
                { List.of(DEL_AA_105, SNP_C_107, SNP_C_109), "AAAAAACAC", "6M1D1X1M1X", 0 },

                { List.of(INS_TT_105),  "AAAAAATAAAA", "6M1I4M", 0 },
                { List.of(INS_GGG_106), "AAAAAAAGGAAA", "7M2I3M", 0 },
                { List.of(DEL_AA_100, INS_GGG_106, SNP_C_109), "AAAAAAGGAAC", "1M1D5M2I2M1X", 0 },

                //this tests that SNPS can be inserted immediately prior to (and following) indels
                { List.of( SNP_C_105, DEL_AA_105 ), "AAAAACAAA", "5M1X1D3M", 1 },
                { List.of( SNP_C_100, DEL_AA_100 ), "CAAAAAAAA", "1X1D8M", 1 },
                { List.of( SNP_C_100, DEL_AA_100, SNP_G_102 ), "CGAAAAAAA", "1X1D1X7M", 1 },
                { List.of( SNP_C_105, INS_TT_105 ), "AAAAACTAAAA", "5M1X1I4M", 1 },
                { List.of( SNP_C_100, INS_TT_100, SNP_G_101 ), "CTGAAAAAAAA", "1X1I1X8M", 1 },
                { List.of( SNP_C_100, INS_TT_100, SNP_G_101, SNP_C_105, DEL_AA_105 ), "CTGAAACAAA", "1X1I1X3M1X1D3M", 2 },

                //testing that the logic around anchor bases isn't resulting in variants being dropped accidentally
                { List.of( SNP_C_104, DEL_AA_105 ), "AAAACAAAA", "4M1X1M1D3M", 0 },
                { List.of( SNP_C_104, INS_TT_105 ), "AAAACATAAAA", "4M1X1M1I4M", 0 },

        };
    }
    @Test(dataProvider = "testConstructHaplotypeFromVariantsDataProvider")
    public void basicConstructHaplotypeFromVariants(List<Event> events, String expectedBases, String expectedCigar, int numberOfCompounds) {
        Haplotype ref = new Haplotype("AAAAAAAAAA".getBytes(), true, 500, TextCigarCodec.decode("10M"));
        ref.setGenomeLocation(new SimpleInterval("20", 100, 110));

        Haplotype result = PartiallyDeterminedHaplotypeComputationEngine.constructHaplotypeFromEvents(ref, events, true);
        Assert.assertEquals(result.getBases(), expectedBases.getBytes());
        Assert.assertEquals(result.getCigar(), TextCigarCodec.decode(expectedCigar));

        // Assert that the resulting event map matches the input variants:
        EventMap resultEMap = result.getEventMap();
        // NOTE, because of representation in VCF lines, the compound alleles get compressed into a single in the event map, here we assert that this is correct.
        Assert.assertEquals(resultEMap.getNumberOfEvents(), events.size() - numberOfCompounds);
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void TestOutOfOrderInputs() {
        Haplotype ref = new Haplotype("AAAAAAAAAA".getBytes(), true, 500, TextCigarCodec.decode("10M"));
        ref.setGenomeLocation(new SimpleInterval("20", 100, 110));
        List<Event> variants = List.of(SNP_C_105, SNP_G_105);

        Haplotype result = PartiallyDeterminedHaplotypeComputationEngine.constructHaplotypeFromEvents(ref, variants, true);
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void TestSNPsOverlapping() {
        Haplotype ref = new Haplotype("AAAAAAAAAA".getBytes(), true, 500, TextCigarCodec.decode("10M"));
        ref.setGenomeLocation(new SimpleInterval("20", 100, 110));
        List<Event> events = List.of(SNP_C_109, DEL_AA_100);

        Haplotype result = PartiallyDeterminedHaplotypeComputationEngine.constructHaplotypeFromEvents(ref, events, true);
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void TestVariantNotOverlappingHap() {
        Haplotype ref = new Haplotype("AAAAAAAAAA".getBytes(), true, 500, TextCigarCodec.decode("10M"));
        ref.setGenomeLocation(new SimpleInterval("20", 100, 110));
        List<Event> events = List.of(SNP_C_90);

        Haplotype result = PartiallyDeterminedHaplotypeComputationEngine.constructHaplotypeFromEvents(ref, events, true);
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void TestVariantIndelPartiallyOverlapping() {
        Haplotype ref = new Haplotype("AAAAAAAAAA".getBytes(), true, 500, TextCigarCodec.decode("10M"));
        ref.setGenomeLocation(new SimpleInterval("20", 100, 110));
        List<Event> events = List.of(DEL_AAAAAAA_98);

        Haplotype result = PartiallyDeterminedHaplotypeComputationEngine.constructHaplotypeFromEvents(ref, events, true);
    }

    //This is a test asserting that a real edge case that was prone to cause failures in the PDHMM is handled properly when compound variants are taken into account.
    //(62,Rlen=1,[C])->(82,Rlen=1,[C])->(84,Rlen=13,[C])
    @Test
    public void testMessyAlignmentSite() {
        Haplotype ref = new Haplotype("AAGAAAGATGGAGGCCCAGCCAGATCTGGACCCCACAGGCCGTCTCCCCACACAGCCATTCATGTGGTCTACTTCCAGCCATTCATGTGGTCTATTTCCAAGAAAATAGCCCATCCCCCCAAGATAACACCTTCTCAAAAACTTTACAGCTTTGTGTCTACACTGATATTTAGGTATTTTCTTTCTTTTTTTTTTATGATTAACACATCTAATTCAAGAATATCTTGGCAGGATATTCCCCGCTTAGGAAATG".getBytes(), true, 575, TextCigarCodec.decode("253M"));
        ref.setGenomeLocation(new SimpleInterval("20", 24152646, 24152898));

        final Event e1 = new Event("20", 24152708, Allele.REF_T, Allele.ALT_C);
        final Event e2 = new Event("20", 24152728, Allele.REF_T, Allele.ALT_C);
        final Event e3 = new Event("20", 24152729, Allele.create("CATGTGGTCTATT", true), Allele.ALT_C);

        final List<Event> events = List.of(e1, e2, e3);

        Haplotype result = PartiallyDeterminedHaplotypeComputationEngine.constructHaplotypeFromEvents(ref, events, true);
        Assert.assertEquals(result.getCigar(), TextCigarCodec.decode("62M1X19M1X1M12D157M"));

        // Assert that the resulting event map matches the input variants:
        EventMap resultEMap = result.getEventMap();
        Assert.assertEquals(resultEMap.getNumberOfEvents(), events.size());
        for (Event e : events) {
            Event actualEvent = resultEMap.get(e.getStart());
            Assert.assertNotNull(actualEvent);
            Assert.assertEquals(actualEvent, e);
        }
    }


    @DataProvider
    public Object[][] testGeneratePDHaplotypeDataProvider() {
        return new Object[][] {
                {List.of(SNP_C_105, SNP_C_106), SNP_C_106, false, "AAAAAACAAA", new byte[]{0,0,0,0,0,17,0,0,0,0}, "6M1X3M"},
                {List.of(SNP_C_105, SNP_C_106), SNP_C_106, true , "AAAAAAAAAA", new byte[]{0,0,0,0,0,17,0,0,0,0}, "10M"},

                {List.of(INS_TT_103, SNP_C_105, SNP_C_106), INS_TT_103, false, "AAAATAAAAAA", new byte[]{0,0,0,0,0,0,17,17,0,0,0}, "4M1I6M"},
                {List.of(INS_TT_103, SNP_C_105, SNP_C_106), INS_TT_103, true , "AAAAAAAAAA",  new byte[]{0,0,0,0,0,17,17,0,0,0}, "10M"},
                {List.of(INS_TT_103, SNP_C_105, SNP_C_106), SNP_C_105,  false, "AAAATACAAAA", new byte[]{0,0,0,0,6,0,0,17,0,0,0}, "4M1I1M1X4M"},
                {List.of(INS_TT_103, SNP_C_105, SNP_C_106), SNP_C_105,  true , "AAAATAAAAAA", new byte[]{0,0,0,0,6,0,0,17,0,0,0}, "4M1I6M"},

                {List.of(DEL_AAA_102, SNP_C_105, SNP_C_106), DEL_AAA_102, false, "AAAAAAAA"  , new byte[]{0,0,0,17,17,0,0,0}, "3M2D5M"},
                {List.of(DEL_AAA_102, SNP_C_105, SNP_C_106), DEL_AAA_102, true , "AAAAAAAAAA", new byte[]{0,0,0,0,0,17,17,0,0,0}, "10M"},
                {List.of(DEL_AAA_102, SNP_C_105, SNP_C_106), SNP_C_105,  false,  "AAAAACAAAA", new byte[]{0,0,0,2,4,0,17,0,0,0}, "5M1X4M"},
                {List.of(DEL_AAA_102, SNP_C_105, SNP_C_106), SNP_C_105,  true ,  "AAAAAAAAAA", new byte[]{0,0,0,2,4,0,17,0,0,0}, "10M"},
                {List.of(DEL_AAA_102, SNP_C_105, SNP_C_106), SNP_C_106,  false,  "AAAAAACAAA", new byte[]{0,0,0,2,4,17,0,0,0,0}, "6M1X3M"},
                {List.of(DEL_AAA_102, SNP_C_105, SNP_C_106), SNP_C_106,  true ,  "AAAAAAAAAA", new byte[]{0,0,0,2,4,17,0,0,0,0}, "10M"},

                // making sure we support "complex alleles" from DRAGEN
                {List.of(DEL_AAA_102, SNP_C_105, SNP_C_106, INS_GGG_106), SNP_C_105,  false ,  "AAAAACAGGAAA", new byte[]{0,0,0,2,4,0,17,2,4,0,0,0}, "5M1X1M2I3M"},
                {List.of(DEL_AAA_102, SNP_C_105, SNP_C_106, SNP_T_106, INS_GGG_106), SNP_C_105,  true ,  "AAAAAAAGGAAA", new byte[]{0,0,0,2,4,0,81,2,4,0,0,0}, "7M2I3M"},
                {List.of(DEL_AAA_102, SNP_C_105, SNP_C_106, INS_GGG_106), DEL_AAA_102,  false ,  "AAAAAGGAAA", new byte[]{0,0,0,17,17,2,4,0,0,0}, "3M2D2M2I3M"},
                {List.of(DEL_AAA_102, SNP_C_105, SNP_C_106, SNP_T_106, INS_GGG_106), DEL_AAA_102,  true ,  "AAAAAAAGGAAA", new byte[]{0,0,0,0,0,17,81,2,4,0,0,0}, "7M2I3M"},
                {List.of(SNP_G_101, SNP_C_105, DEL_AA_105), SNP_G_101,  false ,  "AGAAAAAAAA", new byte[]{0,0,0,0,0,17,6,0,0,0}, "1M1X8M"},
                {List.of(SNP_G_101, SNP_C_105, DEL_AA_105), SNP_G_101,  true ,   "AAAAAAAAAA", new byte[]{0,0,0,0,0,17,6,0,0,0}, "10M"},

        };
    }
    @Test(dataProvider = "testGeneratePDHaplotypeDataProvider")
    public void testGeneratePDHaplotypeFromVariants(List<Event> events, Event targetEvent, boolean useRefBase, String expectedBases, byte[] expectedAltArray, String expectedCigar) {
        Haplotype ref = new Haplotype("AAAAAAAAAA".getBytes(), true, 500, TextCigarCodec.decode("10M"));
        ref.setGenomeLocation(new SimpleInterval("20", 100, 110));

        PartiallyDeterminedHaplotype result = PartiallyDeterminedHaplotypeComputationEngine.createNewPDHaplotypeFromEvents(ref, targetEvent, useRefBase, events);
        Assert.assertEquals(new String(result.getBases()), expectedBases);
        Assert.assertEquals(result.getAlternateBases(), expectedAltArray);
        Assert.assertEquals(result.getCigar(), TextCigarCodec.decode(expectedCigar));
        Assert.assertEquals(result.getDeterminedPosition(), targetEvent.getStart());
    }

    // NOTE: This is an enforcement of a behavior that I consider to be a bug in DRAGEN. Specifically my assumption that we needn't ever concern
    // ourselves with overlapping variants turns out to be false... As it turns out in DRAGEN, they are entirely accepting of constructing a
    // PD haplotype that is REF at bases that underlie a spanning deletion... This means (for example) that if we have a 10 base undetermined
    // deletion from 100-109 and we have a determined ref deletion at position 105-106, that we should STILL construct the haplotype with
    // PD bases from 100-109 even though it means we are assigning that deletion at position 100 to be ref (essentially enforcing that we
    // don't handle spanning deletions). Joint Detection will likely override this behavior in the future.
    @Test
    public void testDeletionUnderlappingDeterminedBases() {
        Haplotype ref = new Haplotype("AAAAAAAAAA".getBytes(), true, 500, TextCigarCodec.decode("10M"));
        ref.setGenomeLocation(new SimpleInterval("20", 100, 110));

        PartiallyDeterminedHaplotype result = PartiallyDeterminedHaplotypeComputationEngine.createNewPDHaplotypeFromEvents(ref, DEL_AA_105, true, List.of(DEL_AAAAAAA_102, DEL_AA_105));
        Assert.assertEquals(new String(result.getBases()), "AAAAAAAAAA");
        Assert.assertEquals(result.getAlternateBases(), new byte[]{0,0,0,2,0,0,0,0,4,0});
        Assert.assertEquals(result.getCigar(), TextCigarCodec.decode("10M"));
        Assert.assertEquals(result.getDeterminedPosition(), DEL_AA_105.getStart());
    }

    @Test
    public void testEventsOverlapForPDHapsCode() {
        final BiPredicate<Event, Event> overlaps = PartiallyDeterminedHaplotypeComputationEngine::eventsOverlapForPDHapsCode;

        // easy SNP cases
        Assert.assertFalse(overlaps.test(SNP_C_100, SNP_G_101));
        Assert.assertFalse(overlaps.test(SNP_C_100, SNP_G_102));
        Assert.assertFalse(overlaps.test(SNP_C_107, SNP_G_105));
        Assert.assertTrue(overlaps.test(SNP_C_105, SNP_G_105));
        Assert.assertTrue(overlaps.test(SNP_T_106, SNP_T_106));

        // overlap of SNP and deletion -- note that we add 1 to deletion start but not to deletion end
        Assert.assertFalse(overlaps.test(DEL_AAA_102, SNP_G_101));
        Assert.assertFalse(overlaps.test(DEL_AAA_102, SNP_G_102));
        Assert.assertTrue(overlaps.test(DEL_AAA_102, SNP_C_104));
        Assert.assertFalse(overlaps.test(DEL_AAA_102, SNP_C_105));

        // overlap of SNP and insertion -- note that we add 0.5 to insertion start and end
        Assert.assertFalse(overlaps.test(SNP_G_102, INS_TT_103));
        Assert.assertFalse(overlaps.test(SNP_C_104, INS_TT_103));
        Assert.assertFalse(overlaps.test(SNP_C_105, INS_TT_105));

        // two insertions should overlap only if they occur at the same position
        Assert.assertTrue(overlaps.test(INS_TT_105, INS_TT_105));
        Assert.assertFalse(overlaps.test(INS_TT_105, INS_GGG_106));

        // two deletions
        Assert.assertTrue(overlaps.test(DEL_AAAAAAA_102, DEL_AAA_102));
        Assert.assertTrue(overlaps.test(DEL_AA_105, DEL_AAAAAAA_102));
        Assert.assertFalse(overlaps.test(DEL_AA_100, DEL_AAA_102));

        // deletion and insertion
        Assert.assertFalse(overlaps.test(INS_TT_105, DEL_AA_105));  // add 1 to deletion start but only 0.5 to insertion end
        Assert.assertFalse(overlaps.test(INS_TT_103, DEL_AA_105));
        Assert.assertTrue(overlaps.test(DEL_AAAAAAA_102, INS_GGG_106));
        Assert.assertTrue(overlaps.test(INS_TT_103, DEL_AAA_102));
    }
}