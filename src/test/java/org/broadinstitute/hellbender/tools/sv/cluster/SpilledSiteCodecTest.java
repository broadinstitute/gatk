package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.*;

/**
 * Unit tests for the package-private {@link SVClusterWalker.SpilledSiteCodec} and
 * {@link SVClusterWalker.SpilledSite} encode/decode round-trip, focusing on genotype FORMAT
 * fields (GQ, DP, AD, PL, FT, phasing) that are not exercised by the integration fixtures
 * (which only carry GT/CN/ECN).
 */
public class SpilledSiteCodecTest extends GATKBaseTest {

    // Use the same hg38 dictionary that all other SV tests use.
    private static final SAMSequenceDictionary DICT = SVTestUtils.hg38Dict;

    // -------------------------------------------------------------------------
    // Helper: round-trip a SpilledSite through the codec and return the result.
    // -------------------------------------------------------------------------
    /** Test view pairing the preserved siteSeq with the record reconstructed from the spilled payload. */
    private static final class DecodedSite {
        final int siteSeq;
        final SVCallRecord record;
        DecodedSite(final int siteSeq, final SVCallRecord record) {
            this.siteSeq = siteSeq;
            this.record = record;
        }
    }

    private static DecodedSite roundTrip(final int siteSeq, final SVCallRecord record) {
        final SVClusterWalker.SpilledSiteCodec codec = new SVClusterWalker.SpilledSiteCodec(DICT);
        final ByteArrayOutputStream baos = new ByteArrayOutputStream();
        codec.setOutputStream(baos);
        codec.encode(new SVClusterWalker.SpilledSite(siteSeq,
                SVClusterWalker.SpilledSiteCodec.encodeRecord(record)));

        final SVClusterWalker.SpilledSiteCodec decodeCodec = new SVClusterWalker.SpilledSiteCodec(DICT);
        decodeCodec.setInputStream(new ByteArrayInputStream(baos.toByteArray()));
        final SVClusterWalker.SpilledSite back = decodeCodec.decode();
        if (back == null) {
            return null;
        }
        return new DecodedSite(back.siteSeq,
                SVClusterWalker.SpilledSiteCodec.decodeRecord(back.payload, DICT));
    }

    // -------------------------------------------------------------------------
    // Helper: assert that two genotype lists are identical element-by-element.
    // -------------------------------------------------------------------------
    private static void assertGenotypeLists(final List<Genotype> actual, final List<Genotype> expected) {
        Assert.assertEquals(actual.size(), expected.size(), "Genotype list size");
        for (int i = 0; i < expected.size(); i++) {
            VariantContextTestUtils.assertGenotypesAreEqual(actual.get(i), expected.get(i));
        }
    }

    // -------------------------------------------------------------------------
    // Test 1: siteSeq is preserved
    // -------------------------------------------------------------------------
    @Test
    public void testSiteSeqRoundTrip() {
        final SVCallRecord record = new SVCallRecord(
                "site_seq_test", "chr1", 1000, Boolean.TRUE, "chr1", 2000, Boolean.FALSE,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                1001, Collections.emptyList(),
                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(),
                Collections.emptyMap(), Collections.emptySet(), null, DICT);

        final int siteSeq = 42;
        final DecodedSite back = roundTrip(siteSeq, record);
        Assert.assertNotNull(back);
        Assert.assertEquals(back.siteSeq, siteSeq);
        Assert.assertEquals(back.record.getId(), "site_seq_test");
    }

    // -------------------------------------------------------------------------
    // Test 2: GQ and DP round-trip
    // -------------------------------------------------------------------------
    @Test
    public void testGQandDP() {
        final Genotype g = new GenotypeBuilder("sampleA",
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL))
                .GQ(35)
                .DP(50)
                .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 1)
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                .make();

        final SVCallRecord record = new SVCallRecord(
                "gq_dp_test", "chr1", 5000, Boolean.TRUE, "chr1", 6000, Boolean.FALSE,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                1001, Collections.emptyList(),
                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.singletonList(g),
                Collections.emptyMap(), Collections.emptySet(), null, DICT);

        final DecodedSite back = roundTrip(1, record);
        final Genotype decoded = back.record.getGenotypes().get(0);
        Assert.assertEquals(decoded.getGQ(), 35, "GQ");
        Assert.assertEquals(decoded.getDP(), 50, "DP");
    }

    // -------------------------------------------------------------------------
    // Test 3: AD and PL round-trip
    // -------------------------------------------------------------------------
    @Test
    public void testADandPL() {
        final int[] ad = {10, 20};
        final int[] pl = {0, 30, 300};
        final Genotype g = new GenotypeBuilder("sampleB",
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL))
                .AD(ad)
                .PL(pl)
                .GQ(20)
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                .make();

        final SVCallRecord record = new SVCallRecord(
                "ad_pl_test", "chr1", 7000, Boolean.TRUE, "chr1", 8000, Boolean.FALSE,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                1001, Collections.emptyList(),
                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.singletonList(g),
                Collections.emptyMap(), Collections.emptySet(), null, DICT);

        final DecodedSite back = roundTrip(2, record);
        final Genotype decoded = back.record.getGenotypes().get(0);
        Assert.assertEquals(decoded.getAD(), ad, "AD array");
        Assert.assertEquals(decoded.getPL(), pl, "PL array");
    }

    // -------------------------------------------------------------------------
    // Test 4: genotype filters (FT field) round-trip
    // -------------------------------------------------------------------------
    @Test
    public void testGenotypeFilters() {
        final Genotype g = new GenotypeBuilder("sampleC",
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL))
                .filters("PASS")
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                .make();

        final SVCallRecord record = new SVCallRecord(
                "ft_test", "chr1", 9000, Boolean.TRUE, "chr1", 10000, Boolean.FALSE,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                1001, Collections.emptyList(),
                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.singletonList(g),
                Collections.emptyMap(), Collections.emptySet(), null, DICT);

        final DecodedSite back = roundTrip(3, record);
        final Genotype decoded = back.record.getGenotypes().get(0);
        // Compare to the source genotype's value (htsjdk normalizes a "PASS" genotype filter to
        // unfiltered/null); the codec must reproduce whatever the original reports.
        Assert.assertEquals(decoded.getFilters(), g.getFilters(), "FT filter string");
    }

    // -------------------------------------------------------------------------
    // Test 5: non-PASS genotype filter
    // -------------------------------------------------------------------------
    @Test
    public void testNonPassGenotypeFilter() {
        final Genotype g = new GenotypeBuilder("sampleD",
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL))
                .filters("LOW_GQ")
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                .make();

        final SVCallRecord record = new SVCallRecord(
                "ft_fail_test", "chr1", 11000, Boolean.TRUE, "chr1", 12000, Boolean.FALSE,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                1001, Collections.emptyList(),
                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.singletonList(g),
                Collections.emptyMap(), Collections.emptySet(), null, DICT);

        final DecodedSite back = roundTrip(4, record);
        Assert.assertEquals(back.record.getGenotypes().get(0).getFilters(), "LOW_GQ");
    }

    // -------------------------------------------------------------------------
    // Test 6: phased genotype round-trip
    // -------------------------------------------------------------------------
    @Test
    public void testPhasedGenotype() {
        final Genotype g = new GenotypeBuilder("sampleE",
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL))
                .phased(true)
                .GQ(40)
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                .make();

        final SVCallRecord record = new SVCallRecord(
                "phased_test", "chr1", 13000, Boolean.TRUE, "chr1", 14000, Boolean.FALSE,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                1001, Collections.emptyList(),
                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.singletonList(g),
                Collections.emptyMap(), Collections.emptySet(), null, DICT);

        final DecodedSite back = roundTrip(5, record);
        final Genotype decoded = back.record.getGenotypes().get(0);
        Assert.assertTrue(decoded.isPhased(), "Genotype should be phased");
        Assert.assertEquals(decoded.getGQ(), 40);
    }

    // -------------------------------------------------------------------------
    // Test 7: NO_CALL genotype round-trip
    // -------------------------------------------------------------------------
    @Test
    public void testNoCallGenotype() {
        final Genotype g = new GenotypeBuilder("sampleF",
                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL))
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 2)
                .make();

        final SVCallRecord record = new SVCallRecord(
                "nocall_test", "chr1", 15000, Boolean.TRUE, "chr1", 16000, Boolean.FALSE,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                1001, Collections.emptyList(),
                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.singletonList(g),
                Collections.emptyMap(), Collections.emptySet(), null, DICT);

        final DecodedSite back = roundTrip(6, record);
        final Genotype decoded = back.record.getGenotypes().get(0);
        Assert.assertTrue(decoded.isNoCall(), "Genotype should be NO_CALL");
    }

    // -------------------------------------------------------------------------
    // Test 8: extended attributes of mixed types (Integer CN + String EV)
    // -------------------------------------------------------------------------
    @Test
    public void testExtendedAttributes() {
        final Genotype g = new GenotypeBuilder("sampleG",
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL))
                .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 1)
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                .attribute("EV", "PE,SR")
                .make();

        final SVCallRecord record = new SVCallRecord(
                "ext_attr_test", "chr1", 17000, Boolean.TRUE, "chr1", 18000, Boolean.FALSE,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                1001, Collections.emptyList(),
                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.singletonList(g),
                Collections.emptyMap(), Collections.emptySet(), null, DICT);

        final DecodedSite back = roundTrip(7, record);
        final Genotype decoded = back.record.getGenotypes().get(0);
        Assert.assertEquals(decoded.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT), 1);
        Assert.assertEquals(decoded.getExtendedAttribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT), 2);
        Assert.assertEquals(decoded.getExtendedAttribute("EV"), "PE,SR");
    }

    // -------------------------------------------------------------------------
    // Test 9: multi-allelic CNV record (DEL + DUP ALTs) with several genotypes
    //         covering all FORMAT fields together
    // -------------------------------------------------------------------------
    @Test
    public void testMultiAllelicMultiGenotypeFullFields() {
        // carrier: phased, GQ+DP+AD+PL+FT set, CN=1
        final Genotype g1 = new GenotypeBuilder("carrier",
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL))
                .phased(true)
                .GQ(55)
                .DP(30)
                .AD(new int[]{15, 15, 0})
                .PL(new int[]{0, 100, 200, 300, 400, 500})
                .filters("LOW_GQ")
                .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 1)
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                .attribute("EV", "SR")
                .make();

        // dup carrier: unphased, GQ only, CN=3
        final Genotype g2 = new GenotypeBuilder("dup_carrier",
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP))
                .phased(false)
                .GQ(20)
                .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 3)
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                .make();

        // ref: NO_CALL placeholder with CN only
        final Genotype g3 = new GenotypeBuilder("ref_sample",
                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL))
                .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 2)
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                .make();

        final List<Genotype> origGenotypes = Arrays.asList(g1, g2, g3);

        // Multi-allelic: REF_N / SV_SIMPLE_DEL / SV_SIMPLE_DUP
        final SVCallRecord record = new SVCallRecord(
                "multiallelic_full", "chr1", 20000, null, "chr1", 30000, null,
                GATKSVVCFConstants.StructuralVariantAnnotationType.CNV, null, Collections.emptyList(),
                10001, Collections.emptyList(),
                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                origGenotypes,
                Collections.emptyMap(), Collections.emptySet(), null, DICT);

        final int SEQ = 99;
        final DecodedSite back = roundTrip(SEQ, record);

        Assert.assertNotNull(back, "decode() should not return null");
        Assert.assertEquals(back.siteSeq, SEQ, "siteSeq");

        // Site-level fields
        Assert.assertEquals(back.record.getId(), "multiallelic_full");
        Assert.assertEquals(back.record.getContigA(), "chr1");
        Assert.assertEquals(back.record.getPositionA(), 20000);
        Assert.assertEquals(back.record.getContigB(), "chr1");
        Assert.assertEquals(back.record.getPositionB(), 30000);
        Assert.assertEquals(back.record.getType(), GATKSVVCFConstants.StructuralVariantAnnotationType.CNV);
        Assert.assertEquals(back.record.getAlleles(),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP));

        // Genotype-level assertions via VariantContextTestUtils for completeness
        assertGenotypeLists(back.record.getGenotypes(), origGenotypes);

        // Spot-check individual fields on decoded genotypes
        final Genotype dec1 = back.record.getGenotypes().get(0);
        Assert.assertTrue(dec1.isPhased(), "g1 should be phased");
        Assert.assertEquals(dec1.getGQ(), 55);
        Assert.assertEquals(dec1.getDP(), 30);
        Assert.assertEquals(dec1.getAD(), new int[]{15, 15, 0});
        Assert.assertEquals(dec1.getPL(), new int[]{0, 100, 200, 300, 400, 500});
        Assert.assertEquals(dec1.getFilters(), "LOW_GQ");
        Assert.assertEquals(dec1.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT), 1);
        Assert.assertEquals(dec1.getExtendedAttribute("EV"), "SR");

        final Genotype dec2 = back.record.getGenotypes().get(1);
        Assert.assertFalse(dec2.isPhased(), "g2 should be unphased");
        Assert.assertEquals(dec2.getGQ(), 20);

        final Genotype dec3 = back.record.getGenotypes().get(2);
        Assert.assertTrue(dec3.isNoCall(), "g3 should be NO_CALL");
    }

    // -------------------------------------------------------------------------
    // Test 10: decode() returns null at EOF
    // -------------------------------------------------------------------------
    @Test
    public void testDecodeReturnsNullAtEOF() {
        // A real SortingCollection spill stream always carries the ObjectOutputStream header (written by
        // setOutputStream) followed by >=1 record; an empty byte[] has no header and cannot be opened.
        // So encode one record, then assert the SECOND decode() returns null at end-of-stream.
        final SVCallRecord record = new SVCallRecord(
                "eof_test", "chr1", 1000, Boolean.TRUE, "chr1", 2000, Boolean.FALSE,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                1001, Collections.emptyList(),
                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, DICT);

        final SVClusterWalker.SpilledSiteCodec encoder = new SVClusterWalker.SpilledSiteCodec(DICT);
        final ByteArrayOutputStream baos = new ByteArrayOutputStream();
        encoder.setOutputStream(baos);
        encoder.encode(new SVClusterWalker.SpilledSite(0,
                SVClusterWalker.SpilledSiteCodec.encodeRecord(record)));

        final SVClusterWalker.SpilledSiteCodec decoder = new SVClusterWalker.SpilledSiteCodec(DICT);
        decoder.setInputStream(new ByteArrayInputStream(baos.toByteArray()));
        Assert.assertNotNull(decoder.decode(), "first decode() should return the record");
        Assert.assertNull(decoder.decode(), "decode() should return null at EOF");
    }

    // -------------------------------------------------------------------------
    // Test 11: multiple records encoded into the same stream decode correctly
    //          (validates that the block-data framing does not corrupt a multi-record stream
    //          and that a second round-trip is independent)
    // -------------------------------------------------------------------------
    @Test
    public void testMultipleRecordsInStream() {
        final Genotype g1 = new GenotypeBuilder("s1",
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL))
                .GQ(10).DP(20).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2).make();

        final Genotype g2 = new GenotypeBuilder("s2",
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP))
                .GQ(30).DP(40).phased(true).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2).make();

        final SVCallRecord rec1 = new SVCallRecord(
                "rec1", "chr1", 100, Boolean.TRUE, "chr1", 1100, Boolean.FALSE,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                1001, Collections.emptyList(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.singletonList(g1),
                Collections.emptyMap(), Collections.emptySet(), null, DICT);

        final SVCallRecord rec2 = new SVCallRecord(
                "rec2", "chr1", 2000, Boolean.FALSE, "chr1", 3000, Boolean.TRUE,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, null, Collections.emptyList(),
                1001, Collections.emptyList(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP),
                Collections.singletonList(g2),
                Collections.emptyMap(), Collections.emptySet(), null, DICT);

        final SVClusterWalker.SpilledSiteCodec encCodec = new SVClusterWalker.SpilledSiteCodec(DICT);
        final ByteArrayOutputStream baos = new ByteArrayOutputStream();
        encCodec.setOutputStream(baos);
        encCodec.encode(new SVClusterWalker.SpilledSite(10,
                SVClusterWalker.SpilledSiteCodec.encodeRecord(rec1)));
        encCodec.encode(new SVClusterWalker.SpilledSite(20,
                SVClusterWalker.SpilledSiteCodec.encodeRecord(rec2)));

        final SVClusterWalker.SpilledSiteCodec decCodec = new SVClusterWalker.SpilledSiteCodec(DICT);
        decCodec.setInputStream(new ByteArrayInputStream(baos.toByteArray()));

        final SVClusterWalker.SpilledSite back1 = decCodec.decode();
        Assert.assertNotNull(back1);
        Assert.assertEquals(back1.siteSeq, 10);
        final SVCallRecord back1Rec = SVClusterWalker.SpilledSiteCodec.decodeRecord(back1.payload, DICT);
        Assert.assertEquals(back1Rec.getId(), "rec1");
        Assert.assertEquals(back1Rec.getGenotypes().get(0).getGQ(), 10);
        Assert.assertEquals(back1Rec.getGenotypes().get(0).getDP(), 20);

        final SVClusterWalker.SpilledSite back2 = decCodec.decode();
        Assert.assertNotNull(back2);
        Assert.assertEquals(back2.siteSeq, 20);
        final SVCallRecord back2Rec = SVClusterWalker.SpilledSiteCodec.decodeRecord(back2.payload, DICT);
        Assert.assertEquals(back2Rec.getId(), "rec2");
        Assert.assertTrue(back2Rec.getGenotypes().get(0).isPhased());
        Assert.assertEquals(back2Rec.getGenotypes().get(0).getGQ(), 30);
        Assert.assertEquals(back2Rec.getGenotypes().get(0).getDP(), 40);

        // Third decode should be null (EOF)
        Assert.assertNull(decCodec.decode(), "Third decode should be null at EOF");
    }
}
