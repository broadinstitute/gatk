package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * Unit tests for {@link SVClusterWalker.FoldedGenotype}, the memory-compact pass-2 accumulator entry.
 *
 * <ul>
 *   <li>{@link #testCompareMatchesRepresentativeChain} — the cached-key comparator
 *       ({@link SVClusterWalker.FoldedGenotype#compare}) must make the SAME fold decision as the real
 *       {@link CanonicalSVCollapser#getRepresentativeGenotype} chain on randomized genotype pairs.</li>
 *   <li>{@link #testCompactCodecRoundTrip} — the compact (de)serialization must reproduce every
 *       genotype field, so the finalized output genotype is byte-identical.</li>
 * </ul>
 */
public final class SVClusterFoldedGenotypeTest extends GATKBaseTest {

    private static final CanonicalSVCollapser COLLAPSER = new CanonicalSVCollapser(
            SVTestUtils.hg38Reference,
            CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE,
            CanonicalSVCollapser.BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END,
            CanonicalSVCollapser.FlagFieldLogic.OR);

    private static final Allele REF = Allele.create("N", true);
    private static final Allele ALT_DEL = Allele.create("<DEL>", false);
    private static final Allele ALT_DUP = Allele.create("<DUP>", false);
    private static final Allele[] ALLELE_POOL = {Allele.NO_CALL, REF, ALT_DEL, ALT_DUP};

    /** Genotype that exercises the comparator chain: varied alleles + GQ/CNQ/ECN/CN (some absent). */
    private static Genotype randomComparableGenotype(final Random rng, final String sample) {
        final int nAlleles = rng.nextInt(3); // 0, 1, or 2
        final List<Allele> alleles = new ArrayList<>(nAlleles);
        for (int i = 0; i < nAlleles; i++) {
            alleles.add(ALLELE_POOL[rng.nextInt(ALLELE_POOL.length)]);
        }
        final GenotypeBuilder gb = new GenotypeBuilder(sample, alleles);
        if (rng.nextBoolean()) {
            gb.GQ(rng.nextInt(100));
        }
        if (rng.nextBoolean()) {
            gb.attribute(GATKSVVCFConstants.COPY_NUMBER_QUALITY_FORMAT, rng.nextInt(100));
        }
        if (rng.nextBoolean()) {
            gb.attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, rng.nextInt(5));
        }
        if (rng.nextBoolean()) {
            // Mix Integer and String storage to exercise getAttributeAsInt parsing on both sides.
            final int cn = rng.nextInt(6);
            gb.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, rng.nextBoolean() ? cn : String.valueOf(cn));
        }
        return gb.make();
    }

    @Test
    public void testCompareMatchesRepresentativeChain() {
        final Random rng = new Random(8675309L);
        final int n = 4000;
        final List<Genotype> genotypes = new ArrayList<>(n);
        for (int i = 0; i < n; i++) {
            genotypes.add(randomComparableGenotype(rng, "S" + i));
        }
        int decisive = 0;
        for (int i = 0; i < n; i++) {
            final Genotype current = genotypes.get(i);
            final Genotype incoming = genotypes.get((i + 1 + rng.nextInt(n - 1)) % n);

            // Real fold: getRepresentativeGenotype(asList(current, incoming)) keeps the first maximal
            // element, so it returns `incoming` only when `incoming` strictly wins the chain.
            final Genotype representative = COLLAPSER.getRepresentativeGenotype(Arrays.asList(current, incoming));
            final boolean realIncomingWins = representative == incoming;

            final boolean foldIncomingWins = SVClusterWalker.FoldedGenotype.compare(
                    SVClusterWalker.FoldedGenotype.of(incoming),
                    SVClusterWalker.FoldedGenotype.of(current)) > 0;

            Assert.assertEquals(foldIncomingWins, realIncomingWins,
                    "Fold decision diverged from getRepresentativeGenotype.\n  current=" + current
                            + "\n  incoming=" + incoming);
            if (realIncomingWins) {
                decisive++;
            }
        }
        // Sanity: the random data actually exercises both outcomes (not a vacuous all-false pass).
        Assert.assertTrue(decisive > 0 && decisive < n, "Test data did not exercise both fold outcomes");
    }

    /** Genotype that exercises every serialized field, including extended attributes of several types. */
    private static Genotype randomFullGenotype(final Random rng, final String sample) {
        final int nAlleles = 1 + rng.nextInt(2);
        final List<Allele> alleles = new ArrayList<>(nAlleles);
        for (int i = 0; i < nAlleles; i++) {
            alleles.add(ALLELE_POOL[rng.nextInt(ALLELE_POOL.length)]);
        }
        final GenotypeBuilder gb = new GenotypeBuilder(sample, alleles);
        gb.phased(rng.nextBoolean());
        if (rng.nextBoolean()) { gb.GQ(rng.nextInt(100)); }
        if (rng.nextBoolean()) { gb.DP(rng.nextInt(500)); }
        if (rng.nextBoolean()) { gb.AD(new int[]{rng.nextInt(50), rng.nextInt(50)}); }
        if (rng.nextBoolean()) { gb.PL(new int[]{0, rng.nextInt(99), rng.nextInt(99)}); }
        if (rng.nextBoolean()) { gb.filters("LowQual"); }
        // Exercise every attribute-value tag handled by the compact codec.
        if (rng.nextBoolean()) { gb.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, rng.nextInt(6)); } // Integer
        if (rng.nextBoolean()) { gb.attribute("STR_ATTR", "evidence" + rng.nextInt(10)); }              // String
        if (rng.nextBoolean()) { gb.attribute("DBL_ATTR", rng.nextDouble()); }                          // Double
        if (rng.nextBoolean()) { gb.attribute("FLT_ATTR", (float) rng.nextDouble()); }                  // Float
        if (rng.nextBoolean()) { gb.attribute("LNG_ATTR", (long) rng.nextInt()); }                      // Long
        if (rng.nextBoolean()) { gb.attribute("BOOL_ATTR", rng.nextBoolean()); }                        // Boolean
        if (rng.nextBoolean()) { gb.attribute("INTARR_ATTR", new int[]{rng.nextInt(9), rng.nextInt(9)}); } // int[]
        if (rng.nextBoolean()) { gb.attribute("NULL_ATTR", null); }                                     // null
        if (rng.nextBoolean()) { gb.attribute("LIST_ATTR", Arrays.asList(rng.nextInt(10), rng.nextInt(10))); } // List
        if (rng.nextBoolean()) { gb.attribute("STRLIST_ATTR", Arrays.asList("a", "b" + rng.nextInt(5))); }
        return gb.make();
    }

    @Test
    public void testCompactCodecRoundTrip() {
        final Random rng = new Random(424242L);
        for (int i = 0; i < 1000; i++) {
            final Genotype original = randomFullGenotype(rng, "S" + i);
            final Genotype back = SVClusterWalker.FoldedGenotype.of(original).decode();
            assertGenotypeEquals(back, original);
        }
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testUnsupportedAttributeTypeFailsLoud() {
        final GenotypeBuilder gb = new GenotypeBuilder("S1", Arrays.asList(REF));
        gb.attribute("BAD_ATTR", new Object()); // a type the compact codec does not handle
        // Must fail loud rather than silently drop/corrupt the attribute.
        SVClusterWalker.FoldedGenotype.of(gb.make());
    }

    private static void assertGenotypeEquals(final Genotype actual, final Genotype expected) {
        final String ctx = "\n  expected=" + expected + "\n  actual=" + actual;
        Assert.assertEquals(actual.getSampleName(), expected.getSampleName(), "sampleName" + ctx);
        Assert.assertEquals(actual.getAlleles(), expected.getAlleles(), "alleles" + ctx);
        Assert.assertEquals(actual.isPhased(), expected.isPhased(), "phased" + ctx);
        Assert.assertEquals(actual.getGQ(), expected.getGQ(), "GQ" + ctx);
        Assert.assertEquals(actual.getDP(), expected.getDP(), "DP" + ctx);
        Assert.assertEquals(actual.hasAD(), expected.hasAD(), "hasAD" + ctx);
        if (expected.hasAD()) {
            Assert.assertEquals(actual.getAD(), expected.getAD(), "AD" + ctx);
        }
        Assert.assertEquals(actual.hasPL(), expected.hasPL(), "hasPL" + ctx);
        if (expected.hasPL()) {
            Assert.assertEquals(actual.getPL(), expected.getPL(), "PL" + ctx);
        }
        Assert.assertEquals(actual.getFilters(), expected.getFilters(), "filters" + ctx);

        final Map<String, Object> a = actual.getExtendedAttributes();
        final Map<String, Object> e = expected.getExtendedAttributes();
        Assert.assertEquals(a.keySet(), e.keySet(), "attribute keys" + ctx);
        for (final Map.Entry<String, Object> entry : e.entrySet()) {
            final Object ev = entry.getValue();
            final Object av = a.get(entry.getKey());
            if (ev instanceof int[]) {
                Assert.assertEquals((int[]) av, (int[]) ev, "attr " + entry.getKey() + ctx);
            } else {
                Assert.assertEquals(av, ev, "attr " + entry.getKey() + ctx);
            }
        }
    }
}
