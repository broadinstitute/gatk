package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Unit tests for the genotype-light pass-1 mechanism in {@link SVClusterWalker#lowMemStripAndRegister}.
 *
 * <p>These tests assert the mechanism's invariants (no genotype objects when enabled; carrier
 * genotypes retained when disabled; carrier count always correct) without running the full tool
 * pipeline or touching the traversal engine. The actual peak-heap reduction is only validated by
 * the manual real-shard chr20 smoke test (not in CI).</p>
 */
public class SVClusterLowMemMechanismTest extends GATKBaseTest {

    // ===== Minimal concrete SVClusterWalker subclass for unit testing =====

    /**
     * Minimal concrete subclass of the abstract {@link SVClusterWalker} whose only purpose is to
     * let test code invoke the protected {@link SVClusterWalker#lowMemStripAndRegister} method with
     * controlled gating flags. The tool is NEVER started (no {@code onTraversalStart} or command-line
     * parsing) — only the protected helper method is called directly.
     *
     * <p>Both gating methods ({@link #canUseGenotypeLightPass1Items()} and
     * {@link #canDropNonCarrierCnvGenotypes()}) are overridden to return the flag passed at
     * construction, so a single subclass covers both the enabled and disabled cases.</p>
     *
     * <p>{@link #applyRecord} is a no-op; {@link #routing} is field-initialized to an empty
     * {@link ArrayList} in {@link SVClusterWalker}, so {@code lowMemStripAndRegister} can add
     * routing slots without any additional setup.</p>
     */
    @CommandLineProgramProperties(
            summary = "Test-only SVClusterWalker stub",
            oneLineSummary = "Test-only SVClusterWalker stub",
            programGroup = TestProgramGroup.class
    )
    private static final class TestableWalker extends SVClusterWalker {

        private final boolean genotypeLightEnabled;
        private final boolean dropNonCarrierCnvEnabled;

        TestableWalker(final boolean genotypeLightEnabled, final boolean dropNonCarrierCnvEnabled) {
            this.genotypeLightEnabled = genotypeLightEnabled;
            this.dropNonCarrierCnvEnabled = dropNonCarrierCnvEnabled;
        }

        @Override
        protected boolean canUseGenotypeLightPass1Items() {
            return genotypeLightEnabled;
        }

        @Override
        protected boolean canDropNonCarrierCnvGenotypes() {
            return dropNonCarrierCnvEnabled;
        }

        @Override
        public void applyRecord(final SVCallRecord record) {
            // no-op in unit tests — we only call lowMemStripAndRegister directly
        }
    }

    // ===== Helper: build a non-CNV DEL record with N carriers out of M total samples =====

    /**
     * Build a DEL record with {@code totalSamples} genotypes of which {@code numCarriers} carry the
     * alt allele. Carriers are the first {@code numCarriers} samples ("s0", "s1", ...).
     */
    private static SVCallRecord makeDelRecord(final int totalSamples, final int numCarriers) {
        final List<String> allSamples = new ArrayList<>(totalSamples);
        for (int i = 0; i < totalSamples; i++) {
            allSamples.add("s" + i);
        }
        final Set<String> carrierSamples = allSamples.subList(0, numCarriers)
                .stream().collect(Collectors.toSet());
        return SVTestUtils.makeDeletionRecordWithCarriers(allSamples, carrierSamples);
    }

    // ===== Test A: genotype-light ON, non-CNV =====

    /**
     * When {@code canUseGenotypeLightPass1Items()} returns true and the record is a non-CNV DEL,
     * {@code lowMemStripAndRegister} must return a record with:
     * <ul>
     *   <li>zero genotype objects (engine fed no-genotype item to bound memory), and</li>
     *   <li>{@code getCarrierCount()} equal to the number of actual carriers in the input.</li>
     * </ul>
     */
    @Test
    public void testGenotypeLightOnNonCnv_emptyGenotypesCorrectCarrierCount() {
        final int totalSamples = 5;
        final int numCarriers = 2;

        final SVCallRecord input = makeDelRecord(totalSamples, numCarriers);
        Assert.assertEquals(input.getCarrierCount(), numCarriers,
                "Precondition: input record must have the right carrier count");

        final TestableWalker walker = new TestableWalker(/*genotypeLightEnabled=*/true,
                /*dropNonCarrierCnvEnabled=*/false);
        final SVCallRecord stripped = walker.lowMemStripAndRegister(input);

        // The engine sees no genotype objects — this is the mechanism that bounds pass-1 memory.
        Assert.assertEquals(stripped.getGenotypes().size(), 0,
                "Genotype-light ON: stripped record must carry no genotype objects");

        // The carrier count must still be correct so the REPRESENTATIVE tiebreaker works.
        Assert.assertEquals(stripped.getCarrierCount(), numCarriers,
                "Genotype-light ON: getCarrierCount() must equal the number of input carriers");
    }

    // ===== Test B: genotype-light OFF, non-CNV =====

    /**
     * When {@code canUseGenotypeLightPass1Items()} returns false and the record is a non-CNV DEL,
     * {@code lowMemStripAndRegister} must return a record that retains only the carrier genotypes
     * (same as --fast-mode), with {@code getGenotypes().size() == numCarriers} and all genotypes
     * belonging to carrier samples.
     */
    @Test
    public void testGenotypeLightOffNonCnv_retainsCarrierGenotypes() {
        final int totalSamples = 6;
        final int numCarriers = 3;

        final SVCallRecord input = makeDelRecord(totalSamples, numCarriers);
        // Collect the expected carrier sample names from the input.
        final Set<String> expectedCarrierNames = input.getCarrierGenotypeList()
                .stream().map(Genotype::getSampleName).collect(Collectors.toSet());
        Assert.assertEquals(expectedCarrierNames.size(), numCarriers,
                "Precondition: input must have the expected carrier sample set");

        final TestableWalker walker = new TestableWalker(/*genotypeLightEnabled=*/false,
                /*dropNonCarrierCnvEnabled=*/false);
        final SVCallRecord stripped = walker.lowMemStripAndRegister(input);

        // With genotype-light OFF the walker retains carrier genotypes (like --fast-mode).
        Assert.assertEquals(stripped.getGenotypes().size(), numCarriers,
                "Genotype-light OFF: stripped record must retain exactly the carrier genotypes");

        final Set<String> retainedNames = stripped.getGenotypes().stream()
                .map(Genotype::getSampleName).collect(Collectors.toSet());
        Assert.assertEquals(retainedNames, expectedCarrierNames,
                "Genotype-light OFF: retained genotype samples must match original carriers");

        Assert.assertEquals(stripped.getCarrierCount(), numCarriers,
                "Genotype-light OFF: getCarrierCount() must equal the number of carriers");
    }

    // ===== Test C: genotype-light ON + canDropNonCarrierCnvGenotypes ON, CNV record =====

    /**
     * For a CNV record when both {@code canUseGenotypeLightPass1Items()} and
     * {@code canDropNonCarrierCnvGenotypes()} return true (overlap = 0, non-defrag engine), the
     * stripped record must have:
     * <ul>
     *   <li>zero genotype objects (engine-light path), and</li>
     *   <li>{@code getCarrierCount()} equal to the number of CN-carrier samples.</li>
     * </ul>
     */
    @Test
    public void testGenotypeLightOnCnv_canDropNonCarriers_emptyGenotypesCorrectCarrierCount() {
        final int totalSamples = 8;
        final int numCarriers = 3;

        // makeCnvRecord produces a DEL record so GT-based isCarrier holds; the type is DEL not CNV,
        // which exercises the non-CNV code path in lowMemStripAndRegister. To test the CNV branch
        // we build a CNV-typed record with copy-number attributes that satisfy isCarrier.
        // Build a CNV record using getDiploidCNVGenotypeBuilder: CN=1 (< ECN=2) → carrier for DEL.
        final List<Genotype> cnvGenotypes = new ArrayList<>(totalSamples);
        for (int i = 0; i < totalSamples; i++) {
            final String sample = "cnv_s" + i;
            final int cn = (i < numCarriers) ? 1 : 2; // CN<ECN → carrier; CN==ECN → not carrier
            cnvGenotypes.add(SVTestUtils.getDiploidCNVGenotypeBuilder(sample, cn).make());
        }
        final SVCallRecord input = new SVCallRecord(
                "cnv_test", "chr1", 1000, null, "chr1", 5000, null,
                GATKSVVCFConstants.StructuralVariantAnnotationType.CNV,
                null, Collections.emptyList(),
                4001, Collections.emptyList(),
                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Arrays.asList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                cnvGenotypes,
                Collections.emptyMap(), Collections.emptySet(), null,
                SVTestUtils.hg38Dict);

        Assert.assertEquals(input.getCarrierCount(), numCarriers,
                "Precondition: CNV input must have the expected number of carriers");

        final TestableWalker walker = new TestableWalker(/*genotypeLightEnabled=*/true,
                /*dropNonCarrierCnvEnabled=*/true);
        final SVCallRecord stripped = walker.lowMemStripAndRegister(input);

        Assert.assertEquals(stripped.getGenotypes().size(), 0,
                "CNV + genotype-light ON: stripped record must carry no genotype objects");
        Assert.assertEquals(stripped.getCarrierCount(), numCarriers,
                "CNV + genotype-light ON: getCarrierCount() must equal the number of CN-carriers");
    }
}
