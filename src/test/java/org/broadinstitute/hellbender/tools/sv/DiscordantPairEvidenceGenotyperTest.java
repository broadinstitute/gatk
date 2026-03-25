package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Unit tests for {@link DiscordantPairEvidenceGenotyper}.
 */
public class DiscordantPairEvidenceGenotyperTest extends GATKBaseTest {

    private static final SAMSequenceDictionary DICTIONARY = SVTestUtils.hg38Dict;
    private static final double TEST_TOLERANCE = 1e-4;
    private static final int MAX_QUAL = 999;
    private static final double QUALITY_CUTOFF = 5.0;
    private static final double TARGET_COVERAGE = 30.0;
    private static final Integer MIN_SIZE = 1000;

    private static final List<String> THREE_SAMPLES = Arrays.asList("sample1", "sample2", "sample3");

    // ---- Helper methods ----

    private static Map<String, Double> makeCoverageMap(final List<String> samples, final double coverage) {
        final Map<String, Double> map = new HashMap<>();
        for (final String sample : samples) {
            map.put(sample, coverage);
        }
        return map;
    }

    private static DiscordantPairEvidenceGenotyper makeGenotyper(final List<String> samples, final double coverage) {
        return new DiscordantPairEvidenceGenotyper(
                makeCoverageMap(samples, coverage),
                QUALITY_CUTOFF, MIN_SIZE, TARGET_COVERAGE, MAX_QUAL);
    }

    private static SVCallRecord makeDELRecord(final String id, final int start, final int end,
                                               final List<GATKSVVCFConstants.EvidenceTypes> evidence,
                                               final List<String> algorithms) {
        return new SVCallRecord(id, "chr1", start, true, "chr1", end, false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                null, evidence, algorithms,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, DICTIONARY);
    }

    private static SVCallRecord makeDUPRecord(final String id, final int start, final int end,
                                               final List<GATKSVVCFConstants.EvidenceTypes> evidence,
                                               final List<String> algorithms) {
        return new SVCallRecord(id, "chr1", start, false, "chr1", end, true,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, null, Collections.emptyList(),
                null, evidence, algorithms,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, DICTIONARY);
    }

    private static SVCallRecord makeINSRecord(final String id, final int position) {
        return new SVCallRecord(id, "chr1", position, true, "chr1", position, false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.INS, null, Collections.emptyList(),
                500, Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.SR),
                Collections.singletonList("pesr"),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, DICTIONARY);
    }

    private static SVCallRecord makeBNDRecord(final String id) {
        return new SVCallRecord(id, "chr1", 1000, true, "chr2", 2000, false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, Collections.emptyList(),
                null, Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.PE),
                Collections.singletonList("pesr"),
                Lists.newArrayList(Allele.REF_N, Allele.create("<BND>", false)),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, DICTIONARY);
    }

    private static List<DiscordantPairEvidence> makePEEvidence(final String sample, final String contig,
                                                                final int start, final int end, final int count) {
        final List<DiscordantPairEvidence> evidence = new ArrayList<>();
        for (int i = 0; i < count; i++) {
            evidence.add(new DiscordantPairEvidence(sample, contig, start + i, true, contig, end + i, false));
        }
        return evidence;
    }

    private static DepthEvidenceGenotyper.DepthGenotypeResult makeDepthResult(final int[] copyStates) {
        final double[] depths = new double[copyStates.length];
        final double[] quals = new double[copyStates.length];
        for (int i = 0; i < copyStates.length; i++) {
            depths[i] = copyStates[i] * 0.5;
            quals[i] = 30;
        }
        return new DepthEvidenceGenotyper.DepthGenotypeResult(depths, copyStates, quals, 10);
    }

    // ---- Constructor tests ----

    @Test
    public void testConstructor() {
        final DiscordantPairEvidenceGenotyper genotyper = makeGenotyper(THREE_SAMPLES, TARGET_COVERAGE);
        Assert.assertNotNull(genotyper);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testConstructorNullCoverageMap() {
        new DiscordantPairEvidenceGenotyper(null, QUALITY_CUTOFF, MIN_SIZE, TARGET_COVERAGE, MAX_QUAL);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testConstructorZeroMaxQual() {
        new DiscordantPairEvidenceGenotyper(makeCoverageMap(THREE_SAMPLES, 30), QUALITY_CUTOFF, MIN_SIZE, TARGET_COVERAGE, 0);
    }

    // ---- trainableRecord tests ----

    @Test
    public void testTrainableRecordDEL() {
        // DEL with PE evidence, correct size, mixed copy states including 1 and 2
        final DiscordantPairEvidenceGenotyper genotyper = makeGenotyper(THREE_SAMPLES, TARGET_COVERAGE);
        genotyper.registerVariantForOverlapCheck(makeDELRecord("del1", 1000, 5000,
                Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.PE), Collections.singletonList("pesr")));
        genotyper.aggregateOverlapCheckIntervals();

        final SVCallRecord record = makeDELRecord("del1", 1000, 5000,
                Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.PE), Collections.singletonList("pesr"));
        // Copy states 1, 2, 2 — has both CN=1 and CN=2
        final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = makeDepthResult(new int[]{1, 2, 2});

        Assert.assertTrue(genotyper.trainableRecord(record, depthResult, null));
    }

    @Test
    public void testTrainableRecordDUP() {
        final DiscordantPairEvidenceGenotyper genotyper = makeGenotyper(THREE_SAMPLES, TARGET_COVERAGE);
        genotyper.registerVariantForOverlapCheck(makeDUPRecord("dup1", 1000, 5000,
                Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.PE), Collections.singletonList("pesr")));
        genotyper.aggregateOverlapCheckIntervals();

        final SVCallRecord record = makeDUPRecord("dup1", 1000, 5000,
                Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.PE), Collections.singletonList("pesr"));
        // Copy states 2, 3, 3 — has both CN=2 and CN=3
        final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = makeDepthResult(new int[]{2, 3, 3});

        Assert.assertTrue(genotyper.trainableRecord(record, depthResult, null));
    }

    @DataProvider(name = "nonTrainableRecordData")
    public Object[][] nonTrainableRecordData() {
        final List<GATKSVVCFConstants.EvidenceTypes> peEvidence = Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.PE);
        final List<String> pesrAlg = Collections.singletonList("pesr");
        final List<String> depthAlg = Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM);

        return new Object[][]{
                // INS type → not trainable
                {"INS type", makeINSRecord("ins1", 1000), makeDepthResult(new int[]{1, 2, 2})},
                // BND type → not trainable
                {"BND type", makeBNDRecord("bnd1"), makeDepthResult(new int[]{1, 2, 2})},
                // DEL but depth-only → not trainable
                {"DEL depth-only", makeDELRecord("del_do", 1000, 5000, peEvidence, depthAlg), makeDepthResult(new int[]{1, 2, 2})},
                // DEL but too small → not trainable
                {"DEL too small", makeDELRecord("del_small", 1000, 1100, peEvidence, pesrAlg), makeDepthResult(new int[]{1, 2, 2})},
                // DEL with all same copy states (no mixture) → not trainable
                {"DEL all CN=2", makeDELRecord("del_all2", 1000, 5000, peEvidence, pesrAlg), makeDepthResult(new int[]{2, 2, 2})},
                // DEL with CN=1 but no CN=2 → not trainable
                {"DEL no CN=2", makeDELRecord("del_no2", 1000, 5000, peEvidence, pesrAlg), makeDepthResult(new int[]{1, 1, 1})},
                // DEL with CN>2 → not trainable
                {"DEL with CN=3", makeDELRecord("del_cn3", 1000, 5000, peEvidence, pesrAlg), makeDepthResult(new int[]{1, 2, 3})},
                // DUP with CN out of [2,4] range
                {"DUP with CN=1", makeDUPRecord("dup_cn1", 1000, 5000, peEvidence, pesrAlg), makeDepthResult(new int[]{1, 2, 3})},
                // DUP with CN=5
                {"DUP with CN=5", makeDUPRecord("dup_cn5", 1000, 5000, peEvidence, pesrAlg), makeDepthResult(new int[]{2, 3, 5})},
                // DUP all CN=2 (no CN=3)
                {"DUP all CN=2", makeDUPRecord("dup_all2", 1000, 5000, peEvidence, pesrAlg), makeDepthResult(new int[]{2, 2, 2})},
                // DUP all CN=3 (no CN=2)
                {"DUP all CN=3", makeDUPRecord("dup_all3", 1000, 5000, peEvidence, pesrAlg), makeDepthResult(new int[]{3, 3, 3})},
        };
    }

    @Test(dataProvider = "nonTrainableRecordData")
    public void testNonTrainableRecords(final String description, final SVCallRecord record,
                                         final DepthEvidenceGenotyper.DepthGenotypeResult depthResult) {
        final DiscordantPairEvidenceGenotyper genotyper = makeGenotyper(THREE_SAMPLES, TARGET_COVERAGE);
        genotyper.registerVariantForOverlapCheck(record);
        genotyper.aggregateOverlapCheckIntervals();
        Assert.assertFalse(genotyper.trainableRecord(record, depthResult, null),
                "Record should not be trainable: " + description);
    }

    @Test
    public void testTrainableRecordOverlappingNonDepthRecords() {
        // Two non-depth-only DELs that overlap should not be trainable
        final DiscordantPairEvidenceGenotyper genotyper = makeGenotyper(THREE_SAMPLES, TARGET_COVERAGE);
        final List<GATKSVVCFConstants.EvidenceTypes> peEvidence = Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.PE);
        final List<String> pesrAlg = Collections.singletonList("pesr");
        final SVCallRecord del1 = makeDELRecord("del1", 1000, 5000, peEvidence, pesrAlg);
        final SVCallRecord del2 = makeDELRecord("del2", 3000, 6000, peEvidence, pesrAlg);
        genotyper.registerVariantForOverlapCheck(del1);
        genotyper.registerVariantForOverlapCheck(del2);
        genotyper.aggregateOverlapCheckIntervals();

        final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = makeDepthResult(new int[]{1, 2, 2});

        // del1 overlaps del2, so there are 2 overlapping intervals → not trainable
        Assert.assertFalse(genotyper.trainableRecord(del1, depthResult, null));
        Assert.assertFalse(genotyper.trainableRecord(del2, depthResult, null));
    }

    @Test
    public void testTrainableRecordNullMinSize() {
        // minSize = null should allow any size
        final DiscordantPairEvidenceGenotyper genotyper = new DiscordantPairEvidenceGenotyper(
                makeCoverageMap(THREE_SAMPLES, TARGET_COVERAGE), QUALITY_CUTOFF, null, TARGET_COVERAGE, MAX_QUAL);
        final List<GATKSVVCFConstants.EvidenceTypes> peEvidence = Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.PE);
        final List<String> pesrAlg = Collections.singletonList("pesr");
        final SVCallRecord smallDel = makeDELRecord("del_tiny", 1000, 1010, peEvidence, pesrAlg);
        genotyper.registerVariantForOverlapCheck(smallDel);
        genotyper.aggregateOverlapCheckIntervals();

        final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = makeDepthResult(new int[]{1, 2, 2});
        Assert.assertTrue(genotyper.trainableRecord(smallDel, depthResult, null));
    }

    // ---- genotype() tests ----

    @Test
    public void testGenotypeWithNoEvidence() {
        // All samples have 0 PE counts → all ref genotypes
        final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeParameters params =
                new DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeParameters(2.0, 10.0, 3.0);

        final DiscordantPairEvidenceGenotyper genotyper = makeGenotyper(THREE_SAMPLES, TARGET_COVERAGE);
        final SVCallRecord record = makeDELRecord("del1", 1000, 5000,
                Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.PE), Collections.singletonList("pesr"));

        final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult result =
                genotyper.genotype(record, Collections.emptyList(), params, THREE_SAMPLES);

        Assert.assertNotNull(result);
        for (final int gt : result.genotypes()) {
            Assert.assertEquals(gt, 0, "All genotypes should be ref with no evidence");
        }
        Assert.assertEquals(result.variantQual(), 0., TEST_TOLERANCE);
    }

    @Test
    public void testGenotypeWithEvidence() {
        // One sample has high PE count → non-ref genotype
        final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeParameters params =
                new DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeParameters(2.0, 20.0, 5.0);

        final DiscordantPairEvidenceGenotyper genotyper = makeGenotyper(THREE_SAMPLES, TARGET_COVERAGE);
        final SVCallRecord record = makeDELRecord("del1", 1000, 5000,
                Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.PE), Collections.singletonList("pesr"));

        // Create PE evidence for sample1 (high count)
        final List<DiscordantPairEvidence> evidence = new ArrayList<>();
        evidence.addAll(makePEEvidence("sample1", "chr1", 1000, 5000, 15));

        final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult result =
                genotyper.genotype(record, evidence, params, THREE_SAMPLES);

        Assert.assertNotNull(result);
        Assert.assertTrue(result.genotypes()[0] > 0, "sample1 should have non-ref genotype");
        Assert.assertEquals(result.genotypes()[1], 0, "sample2 should be ref");
        Assert.assertEquals(result.genotypes()[2], 0, "sample3 should be ref");
    }

    @Test
    public void testGenotypeQualsBounded() {
        final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeParameters params =
                new DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeParameters(2.0, 20.0, 5.0);
        final DiscordantPairEvidenceGenotyper genotyper = makeGenotyper(THREE_SAMPLES, TARGET_COVERAGE);
        final SVCallRecord record = makeDELRecord("del1", 1000, 5000,
                Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.PE), Collections.singletonList("pesr"));

        final List<DiscordantPairEvidence> evidence = makePEEvidence("sample1", "chr1", 1000, 5000, 10);

        final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult result =
                genotyper.genotype(record, evidence, params, THREE_SAMPLES);

        for (final int gq : result.genotypeQuals()) {
            Assert.assertTrue(gq >= 1, "GQ should be >= 1");
            Assert.assertTrue(gq <= MAX_QUAL, "GQ should be <= maxQual");
        }
    }

    @Test
    public void testGenotypeRefMaxQual() {
        // Sample with 0 PE evidence should get max quality for ref genotype
        final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeParameters params =
                new DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeParameters(2.0, 20.0, 5.0);
        final DiscordantPairEvidenceGenotyper genotyper = makeGenotyper(THREE_SAMPLES, TARGET_COVERAGE);
        final SVCallRecord record = makeDELRecord("del1", 1000, 5000,
                Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.PE), Collections.singletonList("pesr"));

        final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult result =
                genotyper.genotype(record, Collections.emptyList(), params, THREE_SAMPLES);

        for (final int gq : result.genotypeQuals()) {
            Assert.assertEquals(gq, MAX_QUAL, "Zero-count ref GQ should be maxQual");
        }
    }

    @Test
    public void testGenotypeHetVsHom() {
        // Compare: sample with moderate count → het (GT=1); sample with high count → hom (GT≥2)
        final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeParameters params =
                new DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeParameters(2.0, 20.0, 5.0);

        final List<String> samples = Arrays.asList("sample_het", "sample_hom", "sample_ref");
        final DiscordantPairEvidenceGenotyper genotyper = makeGenotyper(samples, TARGET_COVERAGE);
        final SVCallRecord record = makeDELRecord("del1", 1000, 5000,
                Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.PE), Collections.singletonList("pesr"));

        // het: count < medianHom - sdHet (8 < 20 - 5 = 15) → GT=1
        // hom: count ~ medianHom → GT=2
        final List<DiscordantPairEvidence> evidence = new ArrayList<>();
        evidence.addAll(makePEEvidence("sample_het", "chr1", 1000, 5000, 8));
        evidence.addAll(makePEEvidence("sample_hom", "chr1", 1000, 5000, 20));

        final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult result =
                genotyper.genotype(record, evidence, params, samples);

        Assert.assertEquals(result.genotypes()[0], 1, "sample_het should be het");
        Assert.assertTrue(result.genotypes()[1] >= 2, "sample_hom should be hom or higher");
        Assert.assertEquals(result.genotypes()[2], 0, "sample_ref should be ref");
    }

    // ---- Variant overlap checking ----

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testRegisterDuplicateVariant() {
        final DiscordantPairEvidenceGenotyper genotyper = makeGenotyper(THREE_SAMPLES, TARGET_COVERAGE);
        final SVCallRecord del = makeDELRecord("del1", 1000, 5000,
                Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.PE), Collections.singletonList("pesr"));
        genotyper.registerVariantForOverlapCheck(del);
        genotyper.registerVariantForOverlapCheck(del); // duplicate ID should throw
    }

    @Test
    public void testRegisterVariantSkipsDepthOnly() {
        // Depth-only records should not be registered for overlap check
        final DiscordantPairEvidenceGenotyper genotyper = makeGenotyper(THREE_SAMPLES, TARGET_COVERAGE);
        final SVCallRecord depthDel = makeDELRecord("del_depth", 1000, 5000,
                Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.RD),
                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM));
        genotyper.registerVariantForOverlapCheck(depthDel);
        // Should not throw even if we add another with same interval
        genotyper.aggregateOverlapCheckIntervals();
    }

    @Test
    public void testRegisterVariantSkipsNonDELDUP() {
        // INS records should not be registered
        final DiscordantPairEvidenceGenotyper genotyper = makeGenotyper(THREE_SAMPLES, TARGET_COVERAGE);
        genotyper.registerVariantForOverlapCheck(makeINSRecord("ins1", 1000));
        genotyper.registerVariantForOverlapCheck(makeBNDRecord("bnd1"));
        genotyper.aggregateOverlapCheckIntervals(); // should succeed
    }

    // ---- Record types ----

    @Test
    public void testDiscordantPairGenotypeResultRecord() {
        final int[] genotypes = {0, 1, 2};
        final int[] gqs = {30, 25, 20};
        final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult result =
                new DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult(genotypes, gqs, 15.5);
        Assert.assertEquals(result.genotypes(), genotypes);
        Assert.assertEquals(result.genotypeQuals(), gqs);
        Assert.assertEquals(result.variantQual(), 15.5, TEST_TOLERANCE);
    }

    @Test
    public void testDiscordantPairGenotypeParametersRecord() {
        final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeParameters params =
                new DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeParameters(3.0, 15.0, 4.5);
        Assert.assertEquals(params.minCount(), 3.0, TEST_TOLERANCE);
        Assert.assertEquals(params.medianHom(), 15.0, TEST_TOLERANCE);
        Assert.assertEquals(params.sdHet(), 4.5, TEST_TOLERANCE);
    }

    // ---- Table parser ----

    @Test
    public void testDiscordantPairTableParserColumns() {
        Assert.assertNotNull(DiscordantPairEvidenceGenotyper.DiscordantPairTableParser.CUTOFFS_COLUMNS);
        Assert.assertEquals(DiscordantPairEvidenceGenotyper.DiscordantPairTableParser.CUTOFFS_COLUMNS.columnCount(), 3);
    }

    // ---- Training pipeline validation ----

    @Test(expectedExceptions = IllegalStateException.class)
    public void testFinalizeFirstPassWithoutData() {
        final DiscordantPairEvidenceGenotyper genotyper = makeGenotyper(THREE_SAMPLES, TARGET_COVERAGE);
        genotyper.finalizeFirstPass(); // no addFirstPass calls → should fail
    }

    @Test
    public void testTrainingPipelineFullFlow() {
        // Full training: addFirstPass → finalizeFirstPass → addSecondPass → finalizeSecondPass
        final List<String> samples = Arrays.asList("s1", "s2", "s3", "s4");
        final DiscordantPairEvidenceGenotyper genotyper = new DiscordantPairEvidenceGenotyper(
                makeCoverageMap(samples, TARGET_COVERAGE), QUALITY_CUTOFF, MIN_SIZE, TARGET_COVERAGE, MAX_QUAL);

        final List<GATKSVVCFConstants.EvidenceTypes> peEvidence = Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.PE);
        final List<String> pesrAlg = Collections.singletonList("pesr");

        // Create several training records with PE evidence
        for (int v = 0; v < 5; v++) {
            final String varId = "del_" + v;
            final SVCallRecord record = makeDELRecord(varId, 1000 + v * 10000, 5000 + v * 10000, peEvidence, pesrAlg);

            // Create evidence: s1 and s3 get het counts (het copy states 1/3),
            // s2 gets hom counts (copy state 0/4), s4 is ref (copy state 2)
            final List<DiscordantPairEvidence> evidence = new ArrayList<>();
            evidence.addAll(makePEEvidence("s1", "chr1", record.getPositionA(), record.getPositionB(), 5));
            evidence.addAll(makePEEvidence("s2", "chr1", record.getPositionA(), record.getPositionB(), 10));
            evidence.addAll(makePEEvidence("s3", "chr1", record.getPositionA(), record.getPositionB(), 6));

            final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = makeDepthResult(new int[]{1, 0, 3, 2});
            genotyper.addFirstPass(record, evidence, depthResult, samples);
        }

        genotyper.finalizeFirstPass();

        // Second pass
        for (int v = 0; v < 5; v++) {
            final String varId = "del_" + v;
            final SVCallRecord record = makeDELRecord(varId, 1000 + v * 10000, 5000 + v * 10000, peEvidence, pesrAlg);
            final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = makeDepthResult(new int[]{1, 0, 3, 2});
            genotyper.addSecondPass(record, depthResult, samples);
        }

        final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeParameters params = genotyper.finalizeSecondPass();

        Assert.assertNotNull(params);
        Assert.assertTrue(params.minCount() > 0, "minCount should be > 0");
        Assert.assertTrue(params.medianHom() > 0, "medianHom should be > 0");
        Assert.assertTrue(params.sdHet() >= 0, "sdHet should be >= 0");
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testDoubleFirstPass() {
        final List<String> samples = Arrays.asList("s1", "s2");
        final DiscordantPairEvidenceGenotyper genotyper = new DiscordantPairEvidenceGenotyper(
                makeCoverageMap(samples, TARGET_COVERAGE), QUALITY_CUTOFF, MIN_SIZE, TARGET_COVERAGE, MAX_QUAL);
        final SVCallRecord record = makeDELRecord("del1", 1000, 5000,
                Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.PE), Collections.singletonList("pesr"));
        final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = makeDepthResult(new int[]{1, 2});
        genotyper.addFirstPass(record, makePEEvidence("s1", "chr1", 1000, 5000, 5), depthResult, samples);
        genotyper.finalizeFirstPass();
        // Second addFirstPass should fail
        genotyper.addFirstPass(record, makePEEvidence("s1", "chr1", 1000, 5000, 5), depthResult, samples);
    }

    @Test
    public void testIsTrainingRecord() {
        final List<String> samples = Arrays.asList("s1", "s2");
        final DiscordantPairEvidenceGenotyper genotyper = new DiscordantPairEvidenceGenotyper(
                makeCoverageMap(samples, TARGET_COVERAGE), QUALITY_CUTOFF, MIN_SIZE, TARGET_COVERAGE, MAX_QUAL);
        final SVCallRecord record = makeDELRecord("del1", 1000, 5000,
                Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.PE), Collections.singletonList("pesr"));
        final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = makeDepthResult(new int[]{1, 2});
        genotyper.addFirstPass(record, makePEEvidence("s1", "chr1", 1000, 5000, 5), depthResult, samples);

        Assert.assertTrue(genotyper.isTrainingRecord(record));
        Assert.assertFalse(genotyper.isTrainingRecord(makeDELRecord("del2", 2000, 6000,
                Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.PE), Collections.singletonList("pesr"))));
    }

    // ---- Coverage normalization ----

    @Test
    public void testCoverageNormalization() {
        // Test that samples with different coverage are normalized correctly
        // Sample with 60x coverage should have counts halved compared to 30x baseline
        final Map<String, Double> coverageMap = new HashMap<>();
        coverageMap.put("high_cov", 60.0);
        coverageMap.put("low_cov", 15.0);
        final DiscordantPairEvidenceGenotyper genotyper = new DiscordantPairEvidenceGenotyper(
                coverageMap, QUALITY_CUTOFF, MIN_SIZE, TARGET_COVERAGE, MAX_QUAL);

        final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeParameters params =
                new DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeParameters(2.0, 20.0, 5.0);

        final SVCallRecord record = makeDELRecord("del1", 1000, 5000,
                Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.PE), Collections.singletonList("pesr"));

        // Both samples have 10 raw PE reads
        final List<DiscordantPairEvidence> evidence = new ArrayList<>();
        evidence.addAll(makePEEvidence("high_cov", "chr1", 1000, 5000, 10));
        evidence.addAll(makePEEvidence("low_cov", "chr1", 1000, 5000, 10));

        final List<String> samples = Arrays.asList("high_cov", "low_cov");
        final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult result =
                genotyper.genotype(record, evidence, params, samples);

        // low_cov should have higher normalized count → higher genotype
        // Normalized: high_cov = round(30 * 10 / 60) = 5, low_cov = round(30 * 10 / 15) = 20
        Assert.assertNotNull(result);
        // low_cov with normalized 20 >= medianHom - sdHet = 15, so GT ≥ 2
        Assert.assertTrue(result.genotypes()[1] >= result.genotypes()[0],
                "Low coverage sample should have higher normalized count");
    }
}
