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
 * Unit tests for {@link SplitReadEvidenceGenotyper}.
 */
public class SplitReadEvidenceGenotyperTest extends GATKBaseTest {

    private static final SAMSequenceDictionary DICTIONARY = SVTestUtils.hg38Dict;
    private static final double TEST_TOLERANCE = 1e-4;
    private static final int MAX_QUAL = 999;
    private static final double QUALITY_CUTOFF = 5.0;
    private static final double TARGET_COVERAGE = 30.0;
    private static final Integer MIN_SIZE = 1000;
    private static final int NUM_SAMPLES = 100;

    private static final List<String> THREE_SAMPLES = Arrays.asList("sample1", "sample2", "sample3");

    // ---- Helper methods ----

    private static Map<String, Double> makeCoverageMap(final List<String> samples, final double coverage) {
        final Map<String, Double> map = new HashMap<>();
        for (final String sample : samples) {
            map.put(sample, coverage);
        }
        return map;
    }

    private static SplitReadEvidenceGenotyper makeGenotyper(final List<String> samples, final int numSamples) {
        return new SplitReadEvidenceGenotyper(
                makeCoverageMap(samples, TARGET_COVERAGE),
                numSamples, QUALITY_CUTOFF, MIN_SIZE, TARGET_COVERAGE, MAX_QUAL);
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

    private static SVCallRecord makeINSRecord(final String id, final int position) {
        return new SVCallRecord(id, "chr1", position, true, "chr1", position, false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.INS, null, Collections.emptyList(),
                500, Arrays.asList(GATKSVVCFConstants.EvidenceTypes.PE, GATKSVVCFConstants.EvidenceTypes.SR),
                Collections.singletonList("pesr"),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, DICTIONARY);
    }

    private static List<SplitReadEvidence> makeSREvidence(final String sample, final String contig,
                                                           final int position, final int count, final boolean strand) {
        return Collections.singletonList(new SplitReadEvidence(sample, contig, position, count, strand));
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
        final SplitReadEvidenceGenotyper genotyper = makeGenotyper(THREE_SAMPLES, NUM_SAMPLES);
        Assert.assertNotNull(genotyper);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testConstructorNullCoverageMap() {
        new SplitReadEvidenceGenotyper(null, NUM_SAMPLES, QUALITY_CUTOFF, MIN_SIZE, TARGET_COVERAGE, MAX_QUAL);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testConstructorZeroMaxQual() {
        new SplitReadEvidenceGenotyper(makeCoverageMap(THREE_SAMPLES, 30), NUM_SAMPLES, QUALITY_CUTOFF, MIN_SIZE, TARGET_COVERAGE, 0);
    }

    // ---- trainableRecord tests ----

    @Test
    public void testTrainableRecordValid() {
        // A record that was registered as PE training record and has SR evidence and meets minSize
        final List<String> samples = Arrays.asList("s1", "s2", "s3");
        final List<GATKSVVCFConstants.EvidenceTypes> peSrEvidence = Arrays.asList(
                GATKSVVCFConstants.EvidenceTypes.PE, GATKSVVCFConstants.EvidenceTypes.SR);
        final List<String> pesrAlg = Collections.singletonList("pesr");

        // Set up a PE genotyper with a registered training record
        final DiscordantPairEvidenceGenotyper peGenotyper = new DiscordantPairEvidenceGenotyper(
                makeCoverageMap(samples, TARGET_COVERAGE), QUALITY_CUTOFF, MIN_SIZE, TARGET_COVERAGE, MAX_QUAL);
        final SVCallRecord record = makeDELRecord("del1", 1000, 5000, peSrEvidence, pesrAlg);
        final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = makeDepthResult(new int[]{1, 2, 2});
        peGenotyper.addFirstPass(record, makePEEvidence("s1", "chr1", 1000, 5000, 5), depthResult, samples);

        // SR genotyper
        final SplitReadEvidenceGenotyper srGenotyper = makeGenotyper(samples, NUM_SAMPLES);
        Assert.assertTrue(srGenotyper.trainableRecord(record, peGenotyper, null));
    }

    @DataProvider(name = "nonTrainableSRRecordData")
    public Object[][] nonTrainableSRRecordData() {
        final List<GATKSVVCFConstants.EvidenceTypes> peSrEvidence = Arrays.asList(
                GATKSVVCFConstants.EvidenceTypes.PE, GATKSVVCFConstants.EvidenceTypes.SR);
        final List<GATKSVVCFConstants.EvidenceTypes> peOnlyEvidence = Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.PE);
        final List<String> pesrAlg = Collections.singletonList("pesr");

        return new Object[][]{
                // Not a PE training record → not trainable
                {"not PE training record", makeDELRecord("del_new", 8000, 12000, peSrEvidence, pesrAlg), false},
                // Too small
                {"too small", makeDELRecord("del1", 1000, 1100, peSrEvidence, pesrAlg), true},
                // No SR evidence
                {"no SR evidence", makeDELRecord("del_nosr", 1000, 5000, peOnlyEvidence, pesrAlg), true},
        };
    }

    @Test(dataProvider = "nonTrainableSRRecordData")
    public void testNonTrainableSRRecords(final String description, final SVCallRecord record, final boolean addToPE) {
        final List<String> samples = Arrays.asList("s1", "s2", "s3");
        final DiscordantPairEvidenceGenotyper peGenotyper = new DiscordantPairEvidenceGenotyper(
                makeCoverageMap(samples, TARGET_COVERAGE), QUALITY_CUTOFF, MIN_SIZE, TARGET_COVERAGE, MAX_QUAL);
        if (addToPE) {
            final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = makeDepthResult(new int[]{1, 2, 2});
            peGenotyper.addFirstPass(record, makePEEvidence("s1", "chr1", record.getPositionA(), record.getPositionB(), 5), depthResult, samples);
        }

        final SplitReadEvidenceGenotyper srGenotyper = makeGenotyper(samples, NUM_SAMPLES);
        Assert.assertFalse(srGenotyper.trainableRecord(record, peGenotyper, null),
                "Record should not be trainable: " + description);
    }

    @Test
    public void testTrainableRecordNullMinSize() {
        // minSize = null should allow any size
        final List<String> samples = Arrays.asList("s1", "s2", "s3");
        final List<GATKSVVCFConstants.EvidenceTypes> peSrEvidence = Arrays.asList(
                GATKSVVCFConstants.EvidenceTypes.PE, GATKSVVCFConstants.EvidenceTypes.SR);
        final List<String> pesrAlg = Collections.singletonList("pesr");

        final DiscordantPairEvidenceGenotyper peGenotyper = new DiscordantPairEvidenceGenotyper(
                makeCoverageMap(samples, TARGET_COVERAGE), QUALITY_CUTOFF, null, TARGET_COVERAGE, MAX_QUAL);
        final SVCallRecord smallRecord = makeDELRecord("del1", 1000, 1010, peSrEvidence, pesrAlg);
        final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = makeDepthResult(new int[]{1, 2, 2});
        peGenotyper.addFirstPass(smallRecord, makePEEvidence("s1", "chr1", 1000, 1010, 5), depthResult, samples);

        final SplitReadEvidenceGenotyper srGenotyper = new SplitReadEvidenceGenotyper(
                makeCoverageMap(samples, TARGET_COVERAGE), NUM_SAMPLES, QUALITY_CUTOFF, null, TARGET_COVERAGE, MAX_QUAL);
        Assert.assertTrue(srGenotyper.trainableRecord(smallRecord, peGenotyper, null));
    }

    // ---- genotype() tests ----

    @Test
    public void testGenotypeNoEvidence() {
        final SplitReadEvidenceGenotyper genotyper = makeGenotyper(THREE_SAMPLES, NUM_SAMPLES);
        final SVCallRecord record = makeDELRecord("del1", 1000, 5000,
                Arrays.asList(GATKSVVCFConstants.EvidenceTypes.PE, GATKSVVCFConstants.EvidenceTypes.SR),
                Collections.singletonList("pesr"));

        final SplitReadEvidenceGenotyper.SplitReadGenotypeParameters params =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeParameters(2.0, 20.0, 5.0);
        final SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs cutoffs =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs(
                        new SplitReadEvidenceGenotyper.CutoffResult(0.5, 0.5, 10, 2, 0, 2),
                        new SplitReadEvidenceGenotyper.CutoffResult(0.5, 0.5, 20, 5, 2, 100)
                );
        final SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics metrics =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics(params, cutoffs);

        final SplitReadEvidenceGenotyper.SplitReadGenotypeResult result =
                genotyper.genotype(record, Collections.emptyList(), Collections.emptyList(), metrics, 15, 2.0, THREE_SAMPLES);

        Assert.assertNotNull(result);
        for (final int gt : result.genotypes()) {
            Assert.assertEquals(gt, 0, "All genotypes should be ref with no evidence");
        }
    }

    @Test
    public void testGenotypeWithStartAndEndEvidence() {
        final SplitReadEvidenceGenotyper genotyper = makeGenotyper(THREE_SAMPLES, NUM_SAMPLES);
        final SVCallRecord record = makeDELRecord("del1", 1000, 5000,
                Arrays.asList(GATKSVVCFConstants.EvidenceTypes.PE, GATKSVVCFConstants.EvidenceTypes.SR),
                Collections.singletonList("pesr"));

        final SplitReadEvidenceGenotyper.SplitReadGenotypeParameters params =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeParameters(2.0, 20.0, 5.0);
        final SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs cutoffs =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs(
                        new SplitReadEvidenceGenotyper.CutoffResult(0.0, 0.0, 10, 2, 0, 2),
                        new SplitReadEvidenceGenotyper.CutoffResult(0.0, 0.0, 20, 5, 2, 100)
                );
        final SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics metrics =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics(params, cutoffs);

        // sample1 has both start and end SR evidence
        final List<SplitReadEvidence> startEvidence = makeSREvidence("sample1", "chr1", 1000, 8, true);
        final List<SplitReadEvidence> endEvidence = makeSREvidence("sample1", "chr1", 5000, 8, false);

        final SplitReadEvidenceGenotyper.SplitReadGenotypeResult result =
                genotyper.genotype(record, startEvidence, endEvidence, metrics, 15, 2.0, THREE_SAMPLES);

        Assert.assertNotNull(result);
        // sample1 should have non-ref genotype due to SR evidence
        Assert.assertTrue(result.genotypes()[0] > 0, "sample1 should have non-ref genotype");
        // sample2 and sample3 have no evidence
        Assert.assertEquals(result.genotypes()[1], 0);
        Assert.assertEquals(result.genotypes()[2], 0);
    }

    @Test
    public void testGenotypeRefMaxQual() {
        // Sample with 0 SR evidence should get max quality
        final SplitReadEvidenceGenotyper genotyper = makeGenotyper(THREE_SAMPLES, NUM_SAMPLES);
        final SVCallRecord record = makeDELRecord("del1", 1000, 5000,
                Arrays.asList(GATKSVVCFConstants.EvidenceTypes.PE, GATKSVVCFConstants.EvidenceTypes.SR),
                Collections.singletonList("pesr"));

        final SplitReadEvidenceGenotyper.SplitReadGenotypeParameters params =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeParameters(2.0, 20.0, 5.0);
        final SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs cutoffs =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs(
                        new SplitReadEvidenceGenotyper.CutoffResult(0.0, 0.0, 10, 2, 0, 2),
                        new SplitReadEvidenceGenotyper.CutoffResult(0.0, 0.0, 20, 5, 2, 100)
                );
        final SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics metrics =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics(params, cutoffs);

        final SplitReadEvidenceGenotyper.SplitReadGenotypeResult result =
                genotyper.genotype(record, Collections.emptyList(), Collections.emptyList(), metrics, 15, 2.0, THREE_SAMPLES);

        // All samples have zero counts → should get maxQual GQ
        for (final int gq : result.genotypeQuals()) {
            Assert.assertEquals(gq, MAX_QUAL);
        }
    }

    @Test
    public void testGenotypeHetVsHom() {
        // Low count → het (GT=1), high count → hom (GT≥2)
        final SplitReadEvidenceGenotyper genotyper = makeGenotyper(THREE_SAMPLES, NUM_SAMPLES);
        final SVCallRecord record = makeDELRecord("del1", 1000, 5000,
                Arrays.asList(GATKSVVCFConstants.EvidenceTypes.PE, GATKSVVCFConstants.EvidenceTypes.SR),
                Collections.singletonList("pesr"));

        final SplitReadEvidenceGenotyper.SplitReadGenotypeParameters params =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeParameters(2.0, 40.0, 10.0);
        final SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs cutoffs =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs(
                        new SplitReadEvidenceGenotyper.CutoffResult(0.0, 0.0, 10, 2, 0, 2),
                        new SplitReadEvidenceGenotyper.CutoffResult(0.0, 0.0, 20, 5, 2, 100)
                );
        final SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics metrics =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics(params, cutoffs);

        // sample1: both-sided, sum=10+10=20, medianHom/2=20, bothsideHomCutoff=2*20=40 → het
        // sample2: both-sided, sum=30+30=60, > bothsideHomCutoff → hom
        final List<SplitReadEvidence> startEvidence = new ArrayList<>();
        startEvidence.addAll(makeSREvidence("sample1", "chr1", 1000, 10, true));
        startEvidence.addAll(makeSREvidence("sample2", "chr1", 1000, 30, true));
        final List<SplitReadEvidence> endEvidence = new ArrayList<>();
        endEvidence.addAll(makeSREvidence("sample1", "chr1", 5000, 10, false));
        endEvidence.addAll(makeSREvidence("sample2", "chr1", 5000, 30, false));

        final SplitReadEvidenceGenotyper.SplitReadGenotypeResult result =
                genotyper.genotype(record, startEvidence, endEvidence, metrics, 15, 2.0, THREE_SAMPLES);

        Assert.assertEquals(result.genotypes()[0], 1, "sample1 should be het (both-sided but moderate count)");
        Assert.assertTrue(result.genotypes()[1] >= 2, "sample2 should be hom (high count)");
        Assert.assertEquals(result.genotypes()[2], 0, "sample3 should be ref");
    }

    @Test
    public void testGenotypeINSUsesMedianHomIns() {
        // For INS type, the medianHomIns parameter should be used instead of params.medianHom
        final SplitReadEvidenceGenotyper genotyper = makeGenotyper(THREE_SAMPLES, NUM_SAMPLES);
        final SVCallRecord insRecord = makeINSRecord("ins1", 1000);

        // Setting high medianHom (100) but low medianHomIns (6) means that INS will use low threshold
        final SplitReadEvidenceGenotyper.SplitReadGenotypeParameters params =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeParameters(1.0, 100.0, 5.0);
        final SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs cutoffs =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs(
                        new SplitReadEvidenceGenotyper.CutoffResult(0.0, 0.0, 10, 2, 0, 2),
                        new SplitReadEvidenceGenotyper.CutoffResult(0.0, 0.0, 20, 5, 2, 100)
                );
        final SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics metrics =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics(params, cutoffs);

        // Evidence with count 4 at start and end → countSum = 8
        // With medianHomIns = 6: medianHet = 3.0, bothsideHomCutoff = 2.0 * 3.0 = 6.0
        // 8 > 6.0 → should be hom
        final List<SplitReadEvidence> startEvidence = makeSREvidence("sample1", "chr1", 1000, 4, true);
        final List<SplitReadEvidence> endEvidence = makeSREvidence("sample1", "chr1", 1000, 4, false);

        final SplitReadEvidenceGenotyper.SplitReadGenotypeResult result =
                genotyper.genotype(insRecord, startEvidence, endEvidence, metrics, 6, 2.0, THREE_SAMPLES);

        Assert.assertTrue(result.genotypes()[0] >= 2,
                "INS with sufficient evidence using medianHomIns should be hom");
    }

    // ---- Record type tests ----

    @Test
    public void testSplitReadGenotypeResultRecord() {
        final int[] genotypes = {0, 1, 2};
        final int[] gqs = {30, 25, 20};
        final SplitReadEvidenceGenotyper.SplitReadGenotypeResult result =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeResult(genotypes, gqs, 15.5, true, false);
        Assert.assertEquals(result.genotypes(), genotypes);
        Assert.assertEquals(result.genotypeQuals(), gqs);
        Assert.assertEquals(result.variantQual(), 15.5, TEST_TOLERANCE);
        Assert.assertTrue(result.bothsidePass());
        Assert.assertFalse(result.backgroundFail());
    }

    @Test
    public void testSplitReadGenotypeParametersRecord() {
        final SplitReadEvidenceGenotyper.SplitReadGenotypeParameters params =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeParameters(3.0, 15.0, 4.5);
        Assert.assertEquals(params.minCount(), 3.0, TEST_TOLERANCE);
        Assert.assertEquals(params.medianHom(), 15.0, TEST_TOLERANCE);
        Assert.assertEquals(params.sdHet(), 4.5, TEST_TOLERANCE);
    }

    @Test
    public void testCutoffResultRecord() {
        final SplitReadEvidenceGenotyper.CutoffResult cutoff =
                new SplitReadEvidenceGenotyper.CutoffResult(0.3, 0.5, 10, 2, 0, 50);
        Assert.assertEquals(cutoff.fracSingle(), 0.3, TEST_TOLERANCE);
        Assert.assertEquals(cutoff.fracBoth(), 0.5, TEST_TOLERANCE);
        Assert.assertEquals(cutoff.countPass(), 10);
        Assert.assertEquals(cutoff.countFail(), 2);
        Assert.assertEquals(cutoff.freqMin(), 0);
        Assert.assertEquals(cutoff.freqMax(), 50);
    }

    @Test
    public void testSplitReadGenotypeFrequencyCutoffsRecord() {
        final SplitReadEvidenceGenotyper.CutoffResult rare =
                new SplitReadEvidenceGenotyper.CutoffResult(0.3, 0.5, 10, 2, 0, 5);
        final SplitReadEvidenceGenotyper.CutoffResult common =
                new SplitReadEvidenceGenotyper.CutoffResult(0.4, 0.6, 20, 5, 5, 100);
        final SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs cutoffs =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs(rare, common);
        Assert.assertEquals(cutoffs.rare(), rare);
        Assert.assertEquals(cutoffs.common(), common);
    }

    @Test
    public void testSplitReadGenotypeMetricsRecord() {
        final SplitReadEvidenceGenotyper.SplitReadGenotypeParameters params =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeParameters(2.0, 15.0, 4.0);
        final SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs cutoffs =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs(
                        new SplitReadEvidenceGenotyper.CutoffResult(0.3, 0.5, 10, 2, 0, 5),
                        new SplitReadEvidenceGenotyper.CutoffResult(0.4, 0.6, 20, 5, 5, 100));
        final SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics metrics =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics(params, cutoffs);
        Assert.assertEquals(metrics.parameters(), params);
        Assert.assertEquals(metrics.cutoffs(), cutoffs);
    }

    // ---- Table parser tests ----

    @Test
    public void testSplitReadTableParserColumns() {
        Assert.assertNotNull(SplitReadEvidenceGenotyper.SplitReadTableParser.CUTOFFS_COLUMNS);
        Assert.assertEquals(SplitReadEvidenceGenotyper.SplitReadTableParser.CUTOFFS_COLUMNS.columnCount(), 15);
    }

    // ---- Training pipeline validation ----

    @Test(expectedExceptions = IllegalStateException.class)
    public void testFinalizeFirstPassWithoutData() {
        final SplitReadEvidenceGenotyper genotyper = makeGenotyper(THREE_SAMPLES, NUM_SAMPLES);
        genotyper.finalizeFirstPass();
    }

    @Test
    public void testTrainingPipelineFullFlow() {
        // Full training: addFirstPass → finalizeFirstPass → addSecondPass → finalizeSecondPass
        final List<String> samples = Arrays.asList("s1", "s2", "s3", "s4");
        final SplitReadEvidenceGenotyper genotyper = new SplitReadEvidenceGenotyper(
                makeCoverageMap(samples, TARGET_COVERAGE),
                samples.size(), QUALITY_CUTOFF, MIN_SIZE, TARGET_COVERAGE, MAX_QUAL);

        final List<GATKSVVCFConstants.EvidenceTypes> peSrEvidence = Arrays.asList(
                GATKSVVCFConstants.EvidenceTypes.PE, GATKSVVCFConstants.EvidenceTypes.SR);
        final List<String> pesrAlg = Collections.singletonList("pesr");

        // Create several training records with SR evidence at both ends
        for (int v = 0; v < 5; v++) {
            final String varId = "del_" + v;
            final int start = 1000 + v * 10000;
            final int end = 5000 + v * 10000;
            final SVCallRecord record = makeDELRecord(varId, start, end, peSrEvidence, pesrAlg);

            // SR evidence for s1 (het) and s2 (hom) at both breakpoints
            final List<SplitReadEvidence> startEvidence = new ArrayList<>();
            startEvidence.addAll(makeSREvidence("s1", "chr1", start, 5, true));
            startEvidence.addAll(makeSREvidence("s2", "chr1", start, 10, true));
            startEvidence.addAll(makeSREvidence("s3", "chr1", start, 4, true));

            final List<SplitReadEvidence> endEvidence = new ArrayList<>();
            endEvidence.addAll(makeSREvidence("s1", "chr1", end, 5, false));
            endEvidence.addAll(makeSREvidence("s2", "chr1", end, 10, false));
            endEvidence.addAll(makeSREvidence("s3", "chr1", end, 4, false));

            final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = makeDepthResult(new int[]{1, 0, 3, 2});
            genotyper.addFirstPass(record, startEvidence, endEvidence, depthResult, samples);
        }

        genotyper.finalizeFirstPass();

        // Second pass
        for (int v = 0; v < 5; v++) {
            final String varId = "del_" + v;
            final SVCallRecord record = makeDELRecord(varId, 1000 + v * 10000, 5000 + v * 10000, peSrEvidence, pesrAlg);
            final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = makeDepthResult(new int[]{1, 0, 3, 2});
            genotyper.addSecondPass(record, depthResult, samples);
        }

        final SplitReadEvidenceGenotyper.SplitReadGenotypeParameters params = genotyper.finalizeSecondPass();

        Assert.assertNotNull(params);
        Assert.assertTrue(params.minCount() > 0, "minCount should be > 0");
        Assert.assertTrue(params.medianHom() > 0, "medianHom should be > 0");
        Assert.assertTrue(params.sdHet() >= 0, "sdHet should be >= 0");
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testDoubleFirstPass() {
        final List<String> samples = Arrays.asList("s1", "s2");
        final SplitReadEvidenceGenotyper genotyper = new SplitReadEvidenceGenotyper(
                makeCoverageMap(samples, TARGET_COVERAGE), samples.size(), QUALITY_CUTOFF, MIN_SIZE, TARGET_COVERAGE, MAX_QUAL);
        final SVCallRecord record = makeDELRecord("del1", 1000, 5000,
                Arrays.asList(GATKSVVCFConstants.EvidenceTypes.PE, GATKSVVCFConstants.EvidenceTypes.SR),
                Collections.singletonList("pesr"));
        final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = makeDepthResult(new int[]{1, 2});
        genotyper.addFirstPass(record,
                makeSREvidence("s1", "chr1", 1000, 5, true),
                makeSREvidence("s1", "chr1", 5000, 5, false),
                depthResult, samples);
        genotyper.finalizeFirstPass();
        // Second addFirstPass should fail
        genotyper.addFirstPass(record,
                makeSREvidence("s1", "chr1", 1000, 5, true),
                makeSREvidence("s1", "chr1", 5000, 5, false),
                depthResult, samples);
    }

    // ---- Coverage normalization ----

    @Test
    public void testCoverageNormalization() {
        // Samples with different coverage should be normalized
        final Map<String, Double> coverageMap = new HashMap<>();
        coverageMap.put("high_cov", 60.0);
        coverageMap.put("low_cov", 15.0);
        coverageMap.put("normal_cov", 30.0);
        final SplitReadEvidenceGenotyper genotyper = new SplitReadEvidenceGenotyper(
                coverageMap, 3, QUALITY_CUTOFF, MIN_SIZE, TARGET_COVERAGE, MAX_QUAL);

        final SplitReadEvidenceGenotyper.SplitReadGenotypeParameters params =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeParameters(2.0, 20.0, 5.0);
        final SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs cutoffs =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs(
                        new SplitReadEvidenceGenotyper.CutoffResult(0.0, 0.0, 10, 2, 0, 2),
                        new SplitReadEvidenceGenotyper.CutoffResult(0.0, 0.0, 20, 5, 2, 100)
                );
        final SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics metrics =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics(params, cutoffs);

        final SVCallRecord record = makeDELRecord("del1", 1000, 5000,
                Arrays.asList(GATKSVVCFConstants.EvidenceTypes.PE, GATKSVVCFConstants.EvidenceTypes.SR),
                Collections.singletonList("pesr"));

        // All samples have 10 raw SR reads, but different coverage
        // Normalized: high_cov = round(30 * 10 / 60) = 5, low_cov = round(30 * 10 / 15) = 20, normal = 10
        final List<SplitReadEvidence> startEvidence = new ArrayList<>();
        startEvidence.addAll(makeSREvidence("high_cov", "chr1", 1000, 10, true));
        startEvidence.addAll(makeSREvidence("low_cov", "chr1", 1000, 10, true));
        startEvidence.addAll(makeSREvidence("normal_cov", "chr1", 1000, 10, true));

        final List<String> samples = Arrays.asList("high_cov", "low_cov", "normal_cov");
        final SplitReadEvidenceGenotyper.SplitReadGenotypeResult result =
                genotyper.genotype(record, startEvidence, Collections.emptyList(), metrics, 15, 2.0, samples);

        Assert.assertNotNull(result);
        // low_cov should have higher normalized count (20) than high_cov (5)
        Assert.assertTrue(result.genotypes()[1] >= result.genotypes()[0],
                "Low coverage sample should have equal or higher genotype due to normalization");
    }

    // ---- backgroundFail and bothsidePass flags ----

    @Test
    public void testBothsidePassWithGoodEvidence() {
        // Records with good bothside SR evidence should get bothsidePass=true
        final SplitReadEvidenceGenotyper genotyper = makeGenotyper(THREE_SAMPLES, NUM_SAMPLES);
        final SVCallRecord record = makeDELRecord("del1", 1000, 5000,
                Arrays.asList(GATKSVVCFConstants.EvidenceTypes.PE, GATKSVVCFConstants.EvidenceTypes.SR),
                Collections.singletonList("pesr"));

        final SplitReadEvidenceGenotyper.SplitReadGenotypeParameters params =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeParameters(2.0, 20.0, 5.0);
        // Set very low cutoffs so everything passes
        final SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs cutoffs =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs(
                        new SplitReadEvidenceGenotyper.CutoffResult(0.0, 0.0, 10, 2, 0, 2),
                        new SplitReadEvidenceGenotyper.CutoffResult(0.0, 0.0, 20, 5, 2, 100)
                );
        final SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics metrics =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics(params, cutoffs);

        // All 3 samples have both-sided evidence
        final List<SplitReadEvidence> startEvidence = new ArrayList<>();
        final List<SplitReadEvidence> endEvidence = new ArrayList<>();
        for (final String sample : THREE_SAMPLES) {
            startEvidence.addAll(makeSREvidence(sample, "chr1", 1000, 5, true));
            endEvidence.addAll(makeSREvidence(sample, "chr1", 5000, 5, false));
        }

        final SplitReadEvidenceGenotyper.SplitReadGenotypeResult result =
                genotyper.genotype(record, startEvidence, endEvidence, metrics, 15, 2.0, THREE_SAMPLES);

        Assert.assertNotNull(result);
        Assert.assertTrue(result.bothsidePass(), "Should have bothsidePass with good both-sided evidence");
        Assert.assertFalse(result.backgroundFail(), "Should not have backgroundFail when passes");
    }
}
