package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Unit tests for {@link SVClusterWalker.SpilledSiteCodec} round-trips covering the CPX path
 * (complex-event intervals) and the sites-only (empty genotype list) path.
 */
public class SpilledSiteCodecComplexTest extends GATKBaseTest {

    private static final SAMSequenceDictionary DICT = SVTestUtils.hg38Dict;

    // Alleles used for CPX records: REF_N + the CPX symbolic allele
    private static final Allele CPX_ALT = SVTestUtils.CPX_ALLELE;

    /**
     * Round-trip a {@link SVClusterWalker.SpilledSite} through the codec.
     */
    private static SVClusterWalker.SpilledSite roundTrip(final int siteSeq,
                                                          final SVCallRecord record) throws Exception {
        final SVClusterWalker.SpilledSiteCodec codec = new SVClusterWalker.SpilledSiteCodec(DICT);
        final ByteArrayOutputStream baos = new ByteArrayOutputStream();
        codec.setOutputStream(baos);
        codec.encode(new SVClusterWalker.SpilledSite(siteSeq, record));

        final SVClusterWalker.SpilledSiteCodec decodeCodec = new SVClusterWalker.SpilledSiteCodec(DICT);
        decodeCodec.setInputStream(new ByteArrayInputStream(baos.toByteArray()));
        return decodeCodec.decode();
    }

    /**
     * Test 1: CPX record with a non-empty list of ComplexEventIntervals and two genotypes.
     *
     * <p>The record uses:
     * <ul>
     *   <li>type = CPX, cpxSubtype = delINVdup (a three-interval subtype)</li>
     *   <li>two CPX intervals: a DEL on chr1 and an INV on chr2</li>
     *   <li>two genotypes (one carrier, one non-carrier)</li>
     * </ul>
     */
    @Test
    public void testCpxRoundTrip() throws Exception {
        // Two CPX intervals on dictionary contigs chr1 and chr2
        final SVCallRecord.ComplexEventInterval interval1 = new SVCallRecord.ComplexEventInterval(
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                new SimpleInterval("chr1", 1100, 1500));
        final SVCallRecord.ComplexEventInterval interval2 = new SVCallRecord.ComplexEventInterval(
                GATKSVVCFConstants.StructuralVariantAnnotationType.INV,
                new SimpleInterval("chr2", 500000, 600000));
        final List<SVCallRecord.ComplexEventInterval> cpxIntervals = Arrays.asList(interval1, interval2);

        // Two genotypes: carrier (het) + non-carrier (hom-ref)
        final Genotype carrier = new GenotypeBuilder("sampleA")
                .alleles(Arrays.asList(Allele.REF_N, CPX_ALT))
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                .make();
        final Genotype nonCarrier = new GenotypeBuilder("sampleB")
                .alleles(Arrays.asList(Allele.REF_N, Allele.REF_N))
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                .make();
        final List<Genotype> genotypes = Arrays.asList(carrier, nonCarrier);

        // CPX record built via the public ctor (with dictionary for coordinate validation).
        final SVCallRecord original = new SVCallRecord(
                "cpx_test_1",
                "chr1", 1000,
                null,   // strandA (null for CPX)
                "chr1", 2000,
                null,   // strandB (null for CPX)
                GATKSVVCFConstants.StructuralVariantAnnotationType.CPX,
                GATKSVVCFConstants.ComplexVariantSubtype.delINVdup,
                cpxIntervals,
                1000,   // length
                Collections.emptyList(),
                Collections.singletonList(SVTestUtils.PESR_ALGORITHM),
                Arrays.asList(Allele.REF_N, CPX_ALT),
                genotypes,
                Collections.emptyMap(),
                Collections.emptySet(),
                null,   // log10PError
                SVTestUtils.hg38Dict
        );

        final int siteSeq = 42;
        final SVClusterWalker.SpilledSite back = roundTrip(siteSeq, original);

        // Verify siteSeq preserved
        Assert.assertEquals(back.siteSeq, siteSeq);

        final SVCallRecord decoded = back.record;

        // Type and subtype
        Assert.assertEquals(decoded.getType(), GATKSVVCFConstants.StructuralVariantAnnotationType.CPX);
        Assert.assertEquals(decoded.getComplexSubtype(), GATKSVVCFConstants.ComplexVariantSubtype.delINVdup);

        // Coordinates
        Assert.assertEquals(decoded.getContigA(), "chr1");
        Assert.assertEquals(decoded.getPositionA(), 1000);
        Assert.assertEquals(decoded.getContigB(), "chr1");
        Assert.assertEquals(decoded.getPositionB(), 2000);

        // ID and length
        Assert.assertEquals(decoded.getId(), "cpx_test_1");
        Assert.assertEquals(decoded.getLength(), Integer.valueOf(1000));

        // CPX intervals — ComplexEventInterval.equals() compares type + interval
        final List<SVCallRecord.ComplexEventInterval> decodedIntervals = decoded.getComplexEventIntervals();
        Assert.assertEquals(decodedIntervals.size(), cpxIntervals.size(),
                "Complex event interval count mismatch");
        // The codec canonicalizes (sorts) the list by encode() string; re-sort originals the same way
        final List<SVCallRecord.ComplexEventInterval> sortedOriginals =
                new java.util.ArrayList<>(cpxIntervals);
        sortedOriginals.sort(java.util.Comparator.comparing(SVCallRecord.ComplexEventInterval::encode));
        for (int i = 0; i < sortedOriginals.size(); i++) {
            Assert.assertEquals(decodedIntervals.get(i), sortedOriginals.get(i),
                    "Mismatch at complex interval index " + i);
        }

        // Alleles
        Assert.assertEquals(decoded.getAlleles(), Arrays.asList(Allele.REF_N, CPX_ALT));

        // Genotypes
        Assert.assertEquals(decoded.getGenotypes().size(), 2);
        Assert.assertEquals(decoded.getGenotypes().get("sampleA").getSampleName(), "sampleA");
        Assert.assertEquals(decoded.getGenotypes().get("sampleB").getSampleName(), "sampleB");
    }

    /**
     * Test 2: Sites-only record with an empty genotype list.
     *
     * <p>Uses a simple DEL type (simplest codec path) with no genotypes, exercising the codec's
     * empty-genotype branch and confirming the decoded record has zero genotypes and site fields intact.
     */
    @Test
    public void testSitesOnlyRoundTrip() throws Exception {
        // Sites-only DEL: no genotypes, no CPX intervals
        final SVCallRecord original = new SVCallRecord(
                "del_sites_only",
                "chr1", 5000,
                Boolean.TRUE,   // strandA for DEL
                "chr1", 15000,
                Boolean.FALSE,  // strandB for DEL
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                null,           // cpxSubtype
                Collections.emptyList(),  // cpxIntervals
                10001,          // length
                Collections.emptyList(),
                Collections.singletonList(SVTestUtils.PESR_ALGORITHM),
                Arrays.asList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(),  // empty genotypes → sites-only
                Collections.emptyMap(),
                Collections.emptySet(),
                null,           // log10PError
                SVTestUtils.hg38Dict
        );

        final int siteSeq = 7;
        final SVClusterWalker.SpilledSite back = roundTrip(siteSeq, original);

        // Verify siteSeq preserved
        Assert.assertEquals(back.siteSeq, siteSeq);

        final SVCallRecord decoded = back.record;

        // Type and subtype
        Assert.assertEquals(decoded.getType(), GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
        Assert.assertNull(decoded.getComplexSubtype());

        // Coordinates
        Assert.assertEquals(decoded.getContigA(), "chr1");
        Assert.assertEquals(decoded.getPositionA(), 5000);
        Assert.assertEquals(decoded.getContigB(), "chr1");
        Assert.assertEquals(decoded.getPositionB(), 15000);

        // ID and length
        Assert.assertEquals(decoded.getId(), "del_sites_only");
        Assert.assertEquals(decoded.getLength(), Integer.valueOf(10001));

        // No CPX intervals
        Assert.assertTrue(decoded.getComplexEventIntervals().isEmpty(),
                "Sites-only record should have no complex event intervals");

        // Alleles
        Assert.assertEquals(decoded.getAlleles(), Arrays.asList(Allele.REF_N, Allele.SV_SIMPLE_DEL));

        // Empty genotype list
        Assert.assertEquals(decoded.getGenotypes().size(), 0,
                "Sites-only record should have empty genotype list");
    }
}
