package org.broadinstitute.hellbender.tools.walkers.sv;

import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterWalker;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Unit tests for the package-private static gating helpers extracted from {@link SVCluster}
 * ({@code canDropNonCarrierCnvGenotypes} and {@code canUseGenotypeLightPass1Items}), and for
 * {@link SVCallRecord#getCarrierCount()} / {@link SVCallRecord#getCarrierGenotypeList()}.
 */
public class SVClusterGatingUnitTest extends GATKBaseTest {

    // -------------------------------------------------------------------------
    // canDropNonCarrierCnvGenotypes / canUseGenotypeLightPass1Items truth tables
    // -------------------------------------------------------------------------

    /**
     * Truth-table rows: {algorithm, defragOverlap, depth, mixed, pesr,
     *                    expectedCanDrop, expectedCanUseLight}
     *
     * Reasoning per row:
     *  1. SINGLE_LINKAGE, all overlaps = 0  → both true
     *  2. MAX_CLIQUE,      all overlaps = 0  → both true
     *  3. SINGLE_LINKAGE, pesr > 0           → both false
     *  4. SINGLE_LINKAGE, mixed > 0          → both false
     *  5. SINGLE_LINKAGE, depth > 0          → both false
     *  6. MAX_CLIQUE,     pesr > 0           → both false
     *  7. DEFRAGMENT_CNV, defrag = 0         → canDrop true, canUseLight false
     *     (DEFRAGMENT_CNV always reads carrier state, so light items are unsafe)
     *  8. DEFRAGMENT_CNV, defrag > 0         → canDrop false, canUseLight false
     */
    @DataProvider(name = "gatingTruthTable")
    public Object[][] gatingTruthTable() {
        return new Object[][]{
                // algorithm,                          defrag, depth, mixed, pesr, canDrop, canUseLight
                {SVClusterWalker.CLUSTER_ALGORITHM.SINGLE_LINKAGE,  0.0,   0.0,  0.0,   0.0,  true,  true},
                {SVClusterWalker.CLUSTER_ALGORITHM.MAX_CLIQUE,      0.0,   0.0,  0.0,   0.0,  true,  true},
                {SVClusterWalker.CLUSTER_ALGORITHM.SINGLE_LINKAGE,  0.0,   0.0,  0.0,   0.5,  false, false},
                {SVClusterWalker.CLUSTER_ALGORITHM.SINGLE_LINKAGE,  0.0,   0.0,  0.5,   0.0,  false, false},
                {SVClusterWalker.CLUSTER_ALGORITHM.SINGLE_LINKAGE,  0.0,   0.5,  0.0,   0.0,  false, false},
                {SVClusterWalker.CLUSTER_ALGORITHM.MAX_CLIQUE,      0.0,   0.0,  0.0,   0.5,  false, false},
                {SVClusterWalker.CLUSTER_ALGORITHM.DEFRAGMENT_CNV,  0.0,   0.0,  0.0,   0.0,  true,  false},
                {SVClusterWalker.CLUSTER_ALGORITHM.DEFRAGMENT_CNV,  0.5,   0.0,  0.0,   0.0,  false, false},
        };
    }

    @Test(dataProvider = "gatingTruthTable")
    public void testCanDropNonCarrierCnvGenotypes(
            final SVClusterWalker.CLUSTER_ALGORITHM algorithm,
            final double defragOverlap,
            final double depthOverlap,
            final double mixedOverlap,
            final double pesrOverlap,
            final boolean expectedCanDrop,
            @SuppressWarnings("unused") final boolean expectedCanUseLight) {

        Assert.assertEquals(
                SVCluster.canDropNonCarrierCnvGenotypes(algorithm, defragOverlap, depthOverlap, mixedOverlap, pesrOverlap),
                expectedCanDrop,
                "canDropNonCarrierCnvGenotypes mismatch for " + algorithm + " defrag=" + defragOverlap
                        + " depth=" + depthOverlap + " mixed=" + mixedOverlap + " pesr=" + pesrOverlap);
    }

    @Test(dataProvider = "gatingTruthTable")
    public void testCanUseGenotypeLightPass1Items(
            final SVClusterWalker.CLUSTER_ALGORITHM algorithm,
            @SuppressWarnings("unused") final double defragOverlap,
            final double depthOverlap,
            final double mixedOverlap,
            final double pesrOverlap,
            @SuppressWarnings("unused") final boolean expectedCanDrop,
            final boolean expectedCanUseLight) {

        Assert.assertEquals(
                SVCluster.canUseGenotypeLightPass1Items(algorithm, depthOverlap, mixedOverlap, pesrOverlap),
                expectedCanUseLight,
                "canUseGenotypeLightPass1Items mismatch for " + algorithm
                        + " depth=" + depthOverlap + " mixed=" + mixedOverlap + " pesr=" + pesrOverlap);
    }

    // -------------------------------------------------------------------------
    // SVCallRecord.getCarrierCount() — non-CNV (GT-based carrier detection)
    // -------------------------------------------------------------------------

    /**
     * Build a DEL record with 3 samples where 2 carry the alt allele and 1 is hom-ref.
     * Verify getCarrierCount() == 2 and equals getCarrierGenotypeList().size().
     */
    @Test
    public void testGetCarrierCountNonCnv() {
        // carrier: het del; non-carrier: hom-ref
        final GenotypeBuilder carrierBuilder1 = new GenotypeBuilder("sample1",
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL))
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2);
        final GenotypeBuilder carrierBuilder2 = new GenotypeBuilder("sample2",
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL))
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2);
        final GenotypeBuilder nonCarrierBuilder = new GenotypeBuilder("sample3",
                Lists.newArrayList(Allele.REF_N, Allele.REF_N))
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2);

        final SVCallRecord record = SVTestUtils.makeRecord(
                "del1",
                "chr1", 1000, Boolean.TRUE,
                "chr1", 2000, Boolean.FALSE,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                1001,
                SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Arrays.asList(carrierBuilder1, carrierBuilder2, nonCarrierBuilder));

        final int expectedCarriers = 2;
        Assert.assertEquals(record.getCarrierCount(), expectedCarriers,
                "getCarrierCount() should equal the number of alt-allele carriers");
        Assert.assertEquals(record.getCarrierGenotypeList().size(), expectedCarriers,
                "getCarrierGenotypeList().size() should match getCarrierCount()");
        Assert.assertEquals(record.getCarrierCount(), record.getCarrierGenotypeList().size(),
                "getCarrierCount() must equal getCarrierGenotypeList().size()");
    }

    // -------------------------------------------------------------------------
    // SVCallRecord.getCarrierCount() — CNV (CN vs ECN carrier detection)
    // -------------------------------------------------------------------------

    /**
     * Build a DEL record where genotypes are no-call but have explicit CN values.
     * Carrier = CN < ECN (deletion).  Non-carrier = CN == ECN.
     *
     * <p>We construct the SVCallRecord directly (bypassing SVTestUtils.makeRecord which
     * calls makeGenotypeWithRefAllele and overwrites the ECN attribute based on allele list size)
     * so we can control the CN and ECN attributes independently.</p>
     */
    @Test
    public void testGetCarrierCountCnvViaCopyNumber() {
        // ECN = 2 (diploid); carrier has CN=1 (del carrier); non-carrier has CN=2
        final Genotype cnCarrier = new GenotypeBuilder("sample1")
                .alleles(Collections.singletonList(Allele.NO_CALL))
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 1)
                .make();
        final Genotype cnNonCarrier1 = new GenotypeBuilder("sample2")
                .alleles(Collections.singletonList(Allele.NO_CALL))
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 2)
                .make();
        final Genotype cnNonCarrier2 = new GenotypeBuilder("sample3")
                .alleles(Collections.singletonList(Allele.NO_CALL))
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 2)
                .make();

        final List<Genotype> genotypes = Arrays.asList(cnCarrier, cnNonCarrier1, cnNonCarrier2);
        final List<Allele> alleles = Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL);

        final SVCallRecord record = new SVCallRecord(
                "del_cn", "chr1", 1000, Boolean.TRUE,
                "chr1", 2000, Boolean.FALSE,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                /* cpxSubtype= */ null,
                /* cpxIntervals= */ Collections.emptyList(),
                /* length= */ 1001,
                /* evidence= */ Collections.emptyList(),
                /* algorithms= */ SVTestUtils.DEPTH_ONLY_ALGORITHM_LIST,
                alleles,
                genotypes,
                Collections.emptyMap(),
                Collections.emptySet(),
                /* log10PError= */ null,
                SVTestUtils.hg38Dict);

        final int expectedCarriers = 1;
        Assert.assertEquals(record.getCarrierCount(), expectedCarriers,
                "getCarrierCount() should detect exactly the CN-based carrier");
        Assert.assertEquals(record.getCarrierGenotypeList().size(), expectedCarriers,
                "getCarrierGenotypeList().size() should match CN-based carrier count");
        Assert.assertEquals(record.getCarrierCount(), record.getCarrierGenotypeList().size(),
                "getCarrierCount() must equal getCarrierGenotypeList().size()");
    }
}
