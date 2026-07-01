package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * Differential test locking the compact CNV copy-state path in {@link SVClusterLinkage#computeSampleOverlap}
 * to the genotype-based path. For randomized record pairs it asserts the overlap fraction is IDENTICAL
 * whether a CNV record carries genotypes or the equivalent {@link CnvSampleCopyState} arrays. This is the
 * byte-identical gate for the pass-1 compact-copy-state change (it touches shared static linkage code).
 */
public final class SVClusterLinkageCompactCopyStateTest extends GATKBaseTest {

    private static final Method COMPUTE_SAMPLE_OVERLAP = computeSampleOverlapMethod();

    private static Method computeSampleOverlapMethod() {
        try {
            final Method m = SVClusterLinkage.class.getDeclaredMethod(
                    "computeSampleOverlap", SVCallRecord.class, SVCallRecord.class);
            m.setAccessible(true);
            return m;
        } catch (final NoSuchMethodException e) {
            throw new IllegalStateException(e);
        }
    }

    private static Double overlap(final SVCallRecord a, final SVCallRecord b) {
        try {
            return (Double) COMPUTE_SAMPLE_OVERLAP.invoke(null, a, b);
        } catch (final ReflectiveOperationException e) {
            throw new IllegalStateException(e);
        }
    }

    private static final String CONTIG = "chr1";

    /** Builds a CN-only reduced genotype like lowMemStripAndRegister does (CN/RD_CN/ECN may each be absent). */
    private static Genotype reducedGenotype(final Random rng, final String sample) {
        final GenotypeBuilder gb = new GenotypeBuilder(sample);
        if (rng.nextInt(10) > 0) { // usually has CN
            gb.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, rng.nextInt(5));
        } else if (rng.nextBoolean()) { // sometimes only RD_CN
            gb.attribute(GATKSVVCFConstants.DEPTH_GENOTYPE_COPY_NUMBER_FORMAT, rng.nextInt(5));
        } // else neither -> getCopyState returns -1
        if (rng.nextInt(10) > 0) {
            gb.attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, rng.nextInt(4));
        }
        return gb.make();
    }

    /** Non-CNV carrier genotype (has GT alleles; may also carry CN/RD_CN as PESR depth fields sometimes do). */
    private static Genotype carrierGenotype(final Random rng, final String sample) {
        final List<Allele> alleles = Arrays.asList(
                Allele.create("N", true),
                Allele.create("<DEL>", false));
        final GenotypeBuilder gb = new GenotypeBuilder(sample, alleles);
        if (rng.nextBoolean()) {
            gb.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, rng.nextInt(5));
        }
        gb.attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2); // required by isCarrier
        return gb.make();
    }

    private static SVCallRecord cnvRecord(final String id, final List<Genotype> genotypes) {
        return new SVCallRecord(id, CONTIG, 1000, null, CONTIG, 2000, null,
                GATKSVVCFConstants.StructuralVariantAnnotationType.CNV, null,
                java.util.Collections.emptyList(), 1001, java.util.Collections.emptyList(),
                Arrays.asList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Arrays.asList(Allele.create("N", true), Allele.create("<CNV>", false)),
                genotypes, java.util.Collections.emptyMap(), java.util.Collections.emptySet(), null,
                SVTestUtils.hg38Dict);
    }

    private static SVCallRecord delRecord(final String id, final List<Genotype> genotypes) {
        return new SVCallRecord(id, CONTIG, 1000, true, CONTIG, 2000, false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null,
                java.util.Collections.emptyList(), 1001, java.util.Collections.emptyList(),
                Arrays.asList("pesr"),
                Arrays.asList(Allele.create("N", true), Allele.create("<DEL>", false)),
                genotypes, java.util.Collections.emptyMap(), java.util.Collections.emptySet(), null,
                SVTestUtils.hg38Dict);
    }

    /** Wraps a CNV record so it exposes compact copy state, mirroring the genotype-light pass-1 item. */
    private static SVCallRecord compactCnv(final SVCallRecord genotypeCnv, final CnvSampleCopyState.Dictionary dict) {
        final CnvSampleCopyState compact = CnvSampleCopyState.fromGenotypes(genotypeCnv.getGenotypes(), dict);
        return new CompactCnvRecord(genotypeCnv, compact);
    }

    /** Minimal CNV record carrying compact copy state and NO genotypes, like a low-mem pass-1 CNV item. */
    private static final class CompactCnvRecord extends SVCallRecord implements CnvCopyStateProvider {
        private final CnvSampleCopyState compact;
        CompactCnvRecord(final SVCallRecord base, final CnvSampleCopyState compact) {
            super(base.getId(), base.getContigA(), base.getPositionA(), base.getStrandA(), base.getContigB(),
                    base.getPositionB(), base.getStrandB(), base.getType(), base.getComplexSubtype(),
                    base.getComplexEventIntervals(), base.getLength(), base.getEvidence(), base.getAlgorithms(),
                    base.getAlleles(), GenotypesContext.NO_GENOTYPES, base.getAttributes(), base.getFilters(),
                    base.getLog10PError());
            this.compact = compact;
        }
        @Override
        public CnvSampleCopyState getCnvCopyState() {
            return compact;
        }
    }

    /** Hom-ref (non-carrier) genotype. */
    private static Genotype refGenotype(final String sample) {
        return new GenotypeBuilder(sample,
                Arrays.asList(Allele.create("N", true), Allele.create("N", true)))
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2) // required by isCarrier
                .make();
    }

    /** Minimal non-CNV record carrying compact sorted carrier indices and NO genotypes, like a pass-1 item. */
    private static final class CompactDelRecord extends SVCallRecord implements CarrierIndexProvider {
        private final int[] carrierIndices;
        private final CnvSampleCopyState.Dictionary dict;
        CompactDelRecord(final SVCallRecord base, final int[] carrierIndices, final CnvSampleCopyState.Dictionary dict) {
            super(base.getId(), base.getContigA(), base.getPositionA(), base.getStrandA(), base.getContigB(),
                    base.getPositionB(), base.getStrandB(), base.getType(), base.getComplexSubtype(),
                    base.getComplexEventIntervals(), base.getLength(), base.getEvidence(), base.getAlgorithms(),
                    base.getAlleles(), GenotypesContext.NO_GENOTYPES, base.getAttributes(), base.getFilters(),
                    base.getLog10PError());
            this.carrierIndices = carrierIndices;
            this.dict = dict;
        }
        @Override
        public int[] getCarrierSampleIndices() {
            return carrierIndices;
        }
        @Override
        public java.util.Set<String> getCarrierSampleSet() {
            final java.util.Set<String> s = new java.util.LinkedHashSet<>();
            for (final int i : carrierIndices) {
                s.add(dict.sampleAt(i));
            }
            return s;
        }
    }

    private static int[] sortedCarrierIndices(final SVCallRecord del, final CnvSampleCopyState.Dictionary dict) {
        return del.getCarrierSampleSet().stream().mapToInt(dict::indexOf).sorted().toArray();
    }

    @Test
    public void testNonCnvCompactMatchesGenotypePath() {
        final Random rng = new Random(20260701L);
        final int nSamples = 300;
        final List<String> allSamples = new ArrayList<>();
        for (int i = 0; i < nSamples; i++) {
            allSamples.add("S" + i);
        }
        final CnvSampleCopyState.Dictionary dict = new CnvSampleCopyState.Dictionary(allSamples);

        int nonTrivial = 0;
        for (int trial = 0; trial < 3000; trial++) {
            final double freqA = rng.nextDouble();
            final double freqB = rng.nextDouble();
            final List<Genotype> aGts = new ArrayList<>();
            final List<Genotype> bGts = new ArrayList<>();
            for (final String s : allSamples) {
                aGts.add(rng.nextDouble() < freqA ? carrierGenotype(rng, s) : refGenotype(s));
                bGts.add(rng.nextDouble() < freqB ? carrierGenotype(rng, s) : refGenotype(s));
            }
            final SVCallRecord aGeno = delRecord("A" + trial, aGts);
            final SVCallRecord bGeno = delRecord("B" + trial, bGts);
            final SVCallRecord aCompact = new CompactDelRecord(aGeno, sortedCarrierIndices(aGeno, dict), dict);
            final SVCallRecord bCompact = new CompactDelRecord(bGeno, sortedCarrierIndices(bGeno, dict), dict);

            final Double expected = overlap(aGeno, bGeno);
            Assert.assertEquals(overlap(aCompact, bCompact), expected, 0.0, "compact != genotype, trial " + trial);
            Assert.assertEquals(overlap(bCompact, aCompact), expected, 0.0, "compact swapped, trial " + trial);
            // Mixed: one compact, one genotype (falls back to getCarrierSampleSet, which CompactDelRecord overrides).
            Assert.assertEquals(overlap(aCompact, bGeno), expected, 0.0, "mixed compact/genotype, trial " + trial);

            if (expected != null && expected > 0.0 && expected < 1.0) {
                nonTrivial++;
            }
        }
        Assert.assertTrue(nonTrivial > 500, "data did not exercise partial carrier overlaps: " + nonTrivial);
    }

    @Test
    public void testCompactMatchesGenotypePath() {
        final Random rng = new Random(1234567L);
        final int nSamples = 200;
        final List<String> allSamples = new ArrayList<>();
        for (int i = 0; i < nSamples; i++) {
            allSamples.add("S" + i);
        }
        final CnvSampleCopyState.Dictionary dict = new CnvSampleCopyState.Dictionary(allSamples);

        int comparisons = 0;
        int nonNull = 0;
        int nonTrivial = 0; // overlap strictly between 0 and 1
        for (int trial = 0; trial < 3000; trial++) {
            // CNV A: all samples (typical) or a random subset (defensive coverage of partial records).
            final boolean aAll = rng.nextInt(4) > 0;
            final List<Genotype> aGts = new ArrayList<>();
            for (final String s : allSamples) {
                if (aAll || rng.nextBoolean()) {
                    aGts.add(reducedGenotype(rng, s));
                }
            }
            final SVCallRecord aGeno = cnvRecord("A" + trial, aGts);
            final SVCallRecord aCompact = compactCnv(aGeno, dict);

            final SVCallRecord bGeno;
            final SVCallRecord bForCompact; // compact form if B is CNV, else the same genotype record
            if (rng.nextBoolean()) {
                // B is another CNV (covers compact-compact)
                final boolean bAll = rng.nextInt(4) > 0;
                final List<Genotype> bGts = new ArrayList<>();
                for (final String s : allSamples) {
                    if (bAll || rng.nextBoolean()) {
                        bGts.add(reducedGenotype(rng, s));
                    }
                }
                bGeno = cnvRecord("B" + trial, bGts);
                bForCompact = compactCnv(bGeno, dict);
            } else {
                // B is a carrier-stripped non-CNV (covers compact-vs-genotype mixed path)
                final List<Genotype> bGts = new ArrayList<>();
                for (final String s : allSamples) {
                    if (rng.nextInt(5) == 0) { // sparse carriers
                        bGts.add(carrierGenotype(rng, s));
                    }
                }
                bGeno = delRecord("B" + trial, bGts);
                bForCompact = bGeno; // non-CNV stays as genotypes
            }

            final Double expected = overlap(aGeno, bGeno);
            // Compact on the A side; and (when B is CNV) on the B side too.
            final Double actual = overlap(aCompact, bForCompact);
            final Double actualSwapped = overlap(bForCompact, aCompact); // symmetry

            Assert.assertEquals(actual, expected, 0.0,
                    "compact != genotype for trial " + trial + " (A all=" + aAll + ")");
            Assert.assertEquals(actualSwapped, expected, 0.0,
                    "compact (swapped args) != genotype for trial " + trial);

            comparisons++;
            if (expected != null) {
                nonNull++;
                if (expected > 0.0 && expected < 1.0) {
                    nonTrivial++;
                }
            }
        }
        // Make sure the randomized data actually exercised meaningful overlaps, not just 0/1/null.
        Assert.assertTrue(nonNull > comparisons / 2, "too many null overlaps: " + nonNull + "/" + comparisons);
        Assert.assertTrue(nonTrivial > 100, "data did not exercise partial overlaps: " + nonTrivial);
    }
}
