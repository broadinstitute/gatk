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

        private static Map<String, Double> makeCoverageMap(final Map<String, Double> coverageBySample) {
                return new HashMap<>(coverageBySample);
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

        private static void addSREvidence(final List<SplitReadEvidence> evidence, final Collection<String> samples,
                                                                          final String contig, final int position, final int count, final boolean strand) {
                for (final String sample : samples) {
                        evidence.addAll(makeSREvidence(sample, contig, position, count, strand));
                }
        }

        private static List<String> combineSamples(final List<String> firstGroup, final List<String> secondGroup) {
                final List<String> samples = new ArrayList<>();
                samples.addAll(firstGroup);
                samples.addAll(secondGroup);
                return samples;
        }

        private static Map<String, Double> makeCoverageMap(final List<String> lowCoverageSamples,
                                                                                                           final double lowCoverage,
                                                                                                           final List<String> highCoverageSamples,
                                                                                                           final double highCoverage) {
                final Map<String, Double> coverageMap = new HashMap<>();
                lowCoverageSamples.forEach(sample -> coverageMap.put(sample, lowCoverage));
                highCoverageSamples.forEach(sample -> coverageMap.put(sample, highCoverage));
                return coverageMap;
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

        @Test
        public void testSplitReadTableParserSupportsLargePassFailCounts() {
                final SplitReadEvidenceGenotyper.SplitReadTableParser parser = new SplitReadEvidenceGenotyper.SplitReadTableParser();
                final SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics metrics = new SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics(
                                new SplitReadEvidenceGenotyper.SplitReadGenotypeParameters(5.0, 98.0, 19.511015999999998),
                                new SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs(
                                                new SplitReadEvidenceGenotyper.CutoffResult(0.0, 0.0, 3_000_000_000L, 4_000_000_000L, 0, 5),
                                                new SplitReadEvidenceGenotyper.CutoffResult(0.1, 0.0, 5_000_000_000L, 6_000_000_000L, 5, 516)
                                )
                );

                final org.broadinstitute.hellbender.utils.tsv.DataLine line = new org.broadinstitute.hellbender.utils.tsv.DataLine(
                        SplitReadEvidenceGenotyper.SplitReadTableParser.CUTOFFS_COLUMNS,
                        IllegalArgumentException::new);
                parser.composeCutoffsLine(metrics, line);
                final SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics parsed = parser.parseTableLine(line);

                Assert.assertEquals(parsed.cutoffs().rare().countPass(), 3_000_000_000L);
                Assert.assertEquals(parsed.cutoffs().rare().countFail(), 4_000_000_000L);
                Assert.assertEquals(parsed.cutoffs().common().countPass(), 5_000_000_000L);
                Assert.assertEquals(parsed.cutoffs().common().countFail(), 6_000_000_000L);
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

    @Test
    public void testGenotypeUsesNormalizedSummedSupportForSingleSidedBackgroundRatio() {
        final List<String> carrierSamples = Arrays.asList("c1", "c2", "c3", "c4");
        final List<String> backgroundSamples = Arrays.asList("b1", "b2", "b3", "b4", "b5", "b6");
                final List<String> samples = combineSamples(carrierSamples, backgroundSamples);
                final Map<String, Double> coverageMap = makeCoverageMap(carrierSamples, 10.0, backgroundSamples, 60.0);

        final SplitReadEvidenceGenotyper genotyper = new SplitReadEvidenceGenotyper(
                makeCoverageMap(coverageMap), 100, QUALITY_CUTOFF, MIN_SIZE, TARGET_COVERAGE, MAX_QUAL);
        final SVCallRecord record = makeDELRecord("del_normalized_background", 1000, 5000,
                Arrays.asList(GATKSVVCFConstants.EvidenceTypes.PE, GATKSVVCFConstants.EvidenceTypes.SR),
                Collections.singletonList("pesr"));

        final SplitReadEvidenceGenotyper.SplitReadGenotypeParameters params =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeParameters(3.0, 20.0, 5.0);
        final SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs cutoffs =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs(
                        new SplitReadEvidenceGenotyper.CutoffResult(1.1, 1.1, 0, 0, 0, 0),
                        new SplitReadEvidenceGenotyper.CutoffResult(0.5, 1.1, 0, 0, 1, 100)
                );
        final SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics metrics =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics(params, cutoffs);

        final List<SplitReadEvidence> startEvidence = new ArrayList<>();
        addSREvidence(startEvidence, carrierSamples, "chr1", 1000, 2, true);
        addSREvidence(startEvidence, backgroundSamples, "chr1", 1000, 2, true);

        final SplitReadEvidenceGenotyper.SplitReadGenotypeResult result =
                genotyper.genotype(record, startEvidence, Collections.emptyList(), metrics, 15, 2.0, samples);

        Assert.assertFalse(result.backgroundFail(),
                "Single-sided background filtering must use normalized summed support so low-coverage carriers are not drowned out by high-coverage background samples");
        Assert.assertFalse(result.bothsidePass(), "One-sided evidence should not trigger bothsidePass");
        for (int i = 0; i < carrierSamples.size(); i++) {
            Assert.assertTrue(result.genotypes()[i] > 0, "Carrier sample should be non-ref: " + carrierSamples.get(i));
        }
        for (int i = carrierSamples.size(); i < samples.size(); i++) {
            Assert.assertEquals(result.genotypes()[i], 0, "High-coverage background sample should remain ref");
        }
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

    @Test
    public void testAccumulateHistogramOnlyUsesNormalizedSummedSupportForCutoffOptimization() {
        final List<String> rareCarrierSamples = Arrays.asList("r1", "r2");
        final List<String> commonCarrierSamples = Arrays.asList("c1", "c2", "c3", "c4");
        final List<String> backgroundSamples = Arrays.asList("b1", "b2", "b3", "b4", "b5", "b6");
                final List<String> carrierSamples = combineSamples(rareCarrierSamples, commonCarrierSamples);
                final List<String> samples = combineSamples(carrierSamples, backgroundSamples);
                final Map<String, Double> coverageMap = makeCoverageMap(carrierSamples, 10.0, backgroundSamples, 60.0);
        final SplitReadEvidenceGenotyper.SplitReadGenotypeParameters params =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeParameters(3.0, 20.0, 5.0);

        final SVCallRecord rarePassRecord = makeDELRecord("rare_pass", 1000, 5000,
                Arrays.asList(GATKSVVCFConstants.EvidenceTypes.PE, GATKSVVCFConstants.EvidenceTypes.SR),
                Collections.singletonList("pesr"));
        final SVCallRecord commonPassRecord = makeDELRecord("common_pass", 6000, 12000,
                Arrays.asList(GATKSVVCFConstants.EvidenceTypes.PE, GATKSVVCFConstants.EvidenceTypes.SR),
                Collections.singletonList("pesr"));
        final SVCallRecord rareFailRecord = makeINSRecord("rare_fail", 20000);
        final SVCallRecord commonFailRecord = makeINSRecord("common_fail", 30000);

        final List<SplitReadEvidence> rarePassEvidence = new ArrayList<>();
        addSREvidence(rarePassEvidence, rareCarrierSamples, "chr1", 1000, 2, true);
        addSREvidence(rarePassEvidence, backgroundSamples, "chr1", 1000, 2, true);

        final List<SplitReadEvidence> commonPassEvidence = new ArrayList<>();
        addSREvidence(commonPassEvidence, commonCarrierSamples, "chr1", 6000, 2, true);
        addSREvidence(commonPassEvidence, backgroundSamples, "chr1", 6000, 2, true);

        final List<SplitReadEvidence> rareFailEvidence = new ArrayList<>();
        addSREvidence(rareFailEvidence, rareCarrierSamples, "chr1", 20000, 2, true);
        addSREvidence(rareFailEvidence, backgroundSamples, "chr1", 20000, 4, true);
        final List<SplitReadEvidence> commonFailEvidence = new ArrayList<>();
        addSREvidence(commonFailEvidence, commonCarrierSamples, "chr1", 30000, 2, true);
        addSREvidence(commonFailEvidence, backgroundSamples, "chr1", 30000, 4, true);

        final SplitReadEvidenceGenotyper cutoffGenotyper = new SplitReadEvidenceGenotyper(
                makeCoverageMap(coverageMap), 100, QUALITY_CUTOFF, MIN_SIZE, TARGET_COVERAGE, MAX_QUAL);
        cutoffGenotyper.accumulateHistogramOnly(rarePassRecord, rarePassEvidence, Collections.emptyList(), true, params, samples);
        cutoffGenotyper.accumulateHistogramOnly(commonPassRecord, commonPassEvidence, Collections.emptyList(), true, params, samples);
        cutoffGenotyper.accumulateHistogramOnly(rareFailRecord, rareFailEvidence, Collections.emptyList(), false, params, samples);
        cutoffGenotyper.accumulateHistogramOnly(commonFailRecord, commonFailEvidence, Collections.emptyList(), false, params, samples);

        final SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs cutoffs = cutoffGenotyper.finalizeThirdPass();

        Assert.assertEquals(cutoffs.rare().fracSingle(), 0.3, TEST_TOLERANCE,
                "Rare single-sided cutoff should separate normalized pass ratio 1.0 from fail ratio 0.2");
        Assert.assertEquals(cutoffs.common().fracSingle(), 0.5, TEST_TOLERANCE,
                "Common single-sided cutoff should separate normalized pass ratio 1.0 from fail ratio 0.4");
    }

    @Test
    public void testGenotypeTrainingUsesNormalizedSummedSupportForCutoffOptimization() {
        final List<String> rareCarrierSamples = Arrays.asList("r1", "r2");
        final List<String> commonCarrierSamples = Arrays.asList("c1", "c2", "c3", "c4");
        final List<String> backgroundSamples = Arrays.asList("b1", "b2", "b3", "b4", "b5", "b6");
                final List<String> carrierSamples = combineSamples(rareCarrierSamples, commonCarrierSamples);
                final List<String> samples = combineSamples(carrierSamples, backgroundSamples);
                final Map<String, Double> passCoverageMap = makeCoverageMap(carrierSamples, 10.0, backgroundSamples, 60.0);
        final SplitReadEvidenceGenotyper genotyper = new SplitReadEvidenceGenotyper(
                makeCoverageMap(passCoverageMap), 100, QUALITY_CUTOFF, MIN_SIZE, TARGET_COVERAGE, MAX_QUAL);
        final SplitReadEvidenceGenotyper.SplitReadGenotypeParameters params =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeParameters(3.0, 20.0, 5.0);
        final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = makeDepthResult(new int[samples.size()]);
        final int[] peGenotypeQuals = new int[samples.size()];
        Arrays.fill(peGenotypeQuals, 1);
        final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult peResult =
                new DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult(new int[samples.size()], peGenotypeQuals, 10.0);

        final SVCallRecord rarePassRecord = makeDELRecord("rare_pass_training", 1000, 5000,
                Arrays.asList(GATKSVVCFConstants.EvidenceTypes.PE, GATKSVVCFConstants.EvidenceTypes.SR),
                Collections.singletonList("pesr"));
        final SVCallRecord commonPassRecord = makeDELRecord("common_pass_training", 6000, 12000,
                Arrays.asList(GATKSVVCFConstants.EvidenceTypes.PE, GATKSVVCFConstants.EvidenceTypes.SR),
                Collections.singletonList("pesr"));
        final SVCallRecord rareFailRecord = makeINSRecord("rare_fail_training", 20000);
        final SVCallRecord commonFailRecord = makeINSRecord("common_fail_training", 30000);

        final List<SplitReadEvidence> rarePassEvidence = new ArrayList<>();
        addSREvidence(rarePassEvidence, rareCarrierSamples, "chr1", 1000, 2, true);
        addSREvidence(rarePassEvidence, backgroundSamples, "chr1", 1000, 2, true);

        final List<SplitReadEvidence> commonPassEvidence = new ArrayList<>();
        addSREvidence(commonPassEvidence, commonCarrierSamples, "chr1", 6000, 2, true);
        addSREvidence(commonPassEvidence, backgroundSamples, "chr1", 6000, 2, true);

        final List<SplitReadEvidence> rareFailEvidence = new ArrayList<>();
        addSREvidence(rareFailEvidence, rareCarrierSamples, "chr1", 20000, 2, true);
        addSREvidence(rareFailEvidence, backgroundSamples, "chr1", 20000, 4, true);

        final List<SplitReadEvidence> commonFailEvidence = new ArrayList<>();
        addSREvidence(commonFailEvidence, commonCarrierSamples, "chr1", 30000, 2, true);
        addSREvidence(commonFailEvidence, backgroundSamples, "chr1", 30000, 4, true);

        genotyper.genotypeTraining(rarePassRecord, rarePassEvidence, Collections.emptyList(), depthResult, peResult, params, samples);
        genotyper.genotypeTraining(commonPassRecord, commonPassEvidence, Collections.emptyList(), depthResult, peResult, params, samples);
        genotyper.genotypeTraining(rareFailRecord, rareFailEvidence, Collections.emptyList(), null, peResult, params, samples);
        genotyper.genotypeTraining(commonFailRecord, commonFailEvidence, Collections.emptyList(), null, peResult, params, samples);

        final SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs cutoffs = genotyper.finalizeThirdPass();

        Assert.assertEquals(cutoffs.rare().fracSingle(), 0.3, TEST_TOLERANCE,
                "Full training path should match normalized rare single-sided cutoff");
        Assert.assertEquals(cutoffs.common().fracSingle(), 0.5, TEST_TOLERANCE,
                "Full training path should match normalized common single-sided cutoff");
    }

    // ---- Regression tests for v1.1 porting correctness ----

    /**
     * Regression test: bothsidePass frequency binning must use twoSidedPassCount (not
     * bothsideNonZeroCount) for the rare/common frequency check, matching v1.1's
     * recover.bothsides.txt column $2 which is the count of samples passing the two-sided
     * SR threshold.
     *
     * <p>Scenario: 5 samples, one has strong both-sided evidence (passes the threshold),
     * three have weak both-sided evidence (both sides nonzero, but below threshold), and one
     * has no evidence. This gives twoSidedPassCount=1 but bothsideNonZeroCount=4.</p>
     *
     * <p>With cutoffs set so that rare freqMax=2: twoSidedPassCount=1 fits in rare
     * (1&le;2), but bothsideNonZeroCount=4 does not (4&gt;2). Using the wrong count causes
     * bothsidePass=false when it should be true.</p>
     */
    @Test
    public void testBothsidePassUsesPassingCountForFrequencyBinning() {
        // 5 samples, all 30x coverage. numSamples=100 → rareMax=2, commonMin=2.
        final List<String> samples = Arrays.asList("s1", "s2", "s3", "s4", "s5");
        final SplitReadEvidenceGenotyper genotyper = makeGenotyper(samples, NUM_SAMPLES);

        final SVCallRecord record = makeDELRecord("del1", 1000, 5000,
                Arrays.asList(GATKSVVCFConstants.EvidenceTypes.PE, GATKSVVCFConstants.EvidenceTypes.SR),
                Collections.singletonList("pesr"));

        // s1: start=3, end=3 → countSum=6, non-ref, twoSided (3>threshold=1 on both sides)
        // s2-s4: start=1, end=1 → countSum=2, non-ref, NOT twoSidedPass (1 not > 1) but
        //        nonzero on both sides (1 > 0)
        // s5: no evidence → ref
        //
        // Result: twoSidedPassCount=1, bothsideNonZeroCount=4, nonRefCount=4
        //         backgroundRatio = 1/4 = 0.25
        final List<SplitReadEvidence> startEvidence = new ArrayList<>();
        final List<SplitReadEvidence> endEvidence = new ArrayList<>();
        startEvidence.addAll(makeSREvidence("s1", "chr1", 1000, 3, true));
        endEvidence.addAll(makeSREvidence("s1", "chr1", 5000, 3, false));
        for (final String sample : Arrays.asList("s2", "s3", "s4")) {
            startEvidence.addAll(makeSREvidence(sample, "chr1", 1000, 1, true));
            endEvidence.addAll(makeSREvidence(sample, "chr1", 5000, 1, false));
        }

        // Construct cutoffs:
        // rare:   fracBoth=0.2 (so 0.25 passes), freqMax=2 (twoSidedPassCount=1 fits; bothsideNonZeroCount=4 does NOT)
        // common: fracBoth=0.5 (so 0.25 fails), freqMin=2
        // fracSingle set above 1.0 so onesidePass is guaranteed false — test depends entirely on bothsidePass
        final SplitReadEvidenceGenotyper.SplitReadGenotypeParameters params =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeParameters(2.0, 20.0, 5.0);
        final SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs cutoffs =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs(
                        new SplitReadEvidenceGenotyper.CutoffResult(1.1, 0.2, 10, 2, 0, 2),
                        new SplitReadEvidenceGenotyper.CutoffResult(1.1, 0.5, 20, 5, 2, 100)
                );
        final SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics metrics =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics(params, cutoffs);

        final SplitReadEvidenceGenotyper.SplitReadGenotypeResult result =
                genotyper.genotype(record, startEvidence, endEvidence, metrics, 15, 2.0, samples);

        Assert.assertNotNull(result);
        // twoSidedPassCount=1 ≤ rareFreqMax=2 and backgroundRatio=0.25 ≥ rareFracBoth=0.2
        // → bothsidePass must be true (v1.1 uses twoSidedPassCount for frequency binning)
        Assert.assertTrue(result.bothsidePass(),
                "bothsidePass should use twoSidedPassCount (1) not bothsideNonZeroCount (4) " +
                "for rare/common frequency binning — this matches v1.1's recover.bothsides.txt $2 column");
        Assert.assertFalse(result.backgroundFail(),
                "backgroundFail should be false when bothsidePass is true");

        // Verify non-ref genotypes: s1-s4 should be non-ref, s5 should be ref
        Assert.assertTrue(result.genotypes()[0] > 0, "s1 should be non-ref");
        Assert.assertTrue(result.genotypes()[1] > 0, "s2 should be non-ref");
        Assert.assertTrue(result.genotypes()[2] > 0, "s3 should be non-ref");
        Assert.assertTrue(result.genotypes()[3] > 0, "s4 should be non-ref");
        Assert.assertEquals(result.genotypes()[4], 0, "s5 should be ref");
    }

    /**
     * Complementary test: when twoSidedPassCount genuinely exceeds the rare frequency max
     * and the ratio doesn't meet the common threshold, bothsidePass should correctly be false.
     */
    @Test
    public void testBothsidePassCorrectlyFalseWhenPassingCountExceedsRareMax() {
        // 5 samples, all 30x. numSamples=100 → rareMax=2, commonMin=2.
        final List<String> samples = Arrays.asList("s1", "s2", "s3", "s4", "s5");
        final SplitReadEvidenceGenotyper genotyper = makeGenotyper(samples, NUM_SAMPLES);

        final SVCallRecord record = makeDELRecord("del1", 1000, 5000,
                Arrays.asList(GATKSVVCFConstants.EvidenceTypes.PE, GATKSVVCFConstants.EvidenceTypes.SR),
                Collections.singletonList("pesr"));

        // All 4 non-ref samples have strong both-sided evidence (count=3 on each side).
        // This gives twoSidedPassCount=4 (all pass > 1 on both sides),
        // bothsideNonZeroCount=4, backgroundRatio=4/4=1.0.
        // twoSidedPassCount=4 > rareMax=2, so rare check fails.
        // For common: ratio=1.0 >= commonBoth=0.5 and twoSidedPassCount=4 >= commonMin=2 → passes common.
        final List<SplitReadEvidence> startEvidence = new ArrayList<>();
        final List<SplitReadEvidence> endEvidence = new ArrayList<>();
        for (final String sample : Arrays.asList("s1", "s2", "s3", "s4")) {
            startEvidence.addAll(makeSREvidence(sample, "chr1", 1000, 3, true));
            endEvidence.addAll(makeSREvidence(sample, "chr1", 5000, 3, false));
        }

        // Set commonBoth very high (1.1) so common check also fails
        final SplitReadEvidenceGenotyper.SplitReadGenotypeParameters params =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeParameters(2.0, 20.0, 5.0);
        final SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs cutoffs =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs(
                        new SplitReadEvidenceGenotyper.CutoffResult(1.1, 0.2, 10, 2, 0, 2),
                        new SplitReadEvidenceGenotyper.CutoffResult(1.1, 1.1, 20, 5, 2, 100)
                );
        final SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics metrics =
                new SplitReadEvidenceGenotyper.SplitReadGenotypeMetrics(params, cutoffs);

        final SplitReadEvidenceGenotyper.SplitReadGenotypeResult result =
                genotyper.genotype(record, startEvidence, endEvidence, metrics, 15, 2.0, samples);

        // twoSidedPassCount=4 > rareMax=2 → rare fails. ratio=1.0 but commonBoth=1.1 → common fails.
        Assert.assertFalse(result.bothsidePass(),
                "bothsidePass should be false when twoSidedPassCount exceeds rare max and common ratio threshold is not met");
        Assert.assertTrue(result.backgroundFail(),
                "backgroundFail should be true when both onesidePass and bothsidePass are false");
    }

    /**
     * Regression test: hasBothSideSupport must use strict {@code >} (not {@code >=}) to be
     * consistent with countBothSideSupport. Both methods port the same v1.1 awk pattern
     * {@code $NF>(sr_count/2)}. Using {@code >=} in hasBothSideSupport but {@code >} in
     * countBothSideSupport would cause a sample at exactly the threshold to be counted as
     * having both-side support in Phase 1 (affecting het/hom calibration) but NOT counted
     * in Phase 2/3 training or genotyping.
     *
     * <p>Scenario: With QUALITY_CUTOFF=5.0, trainingCountCutoff=1, minCount=max(1/2,1)=1.
     * A sample with raw SR count=1 on each side at 30x coverage → normalized=1.0, which
     * equals the threshold exactly. With strict {@code >}, this sample should NOT qualify
     * as having both-side support. A second sample with count=2 on each side does qualify.</p>
     *
     * <p>If hasBothSideSupport used {@code >=}, the borderline sample would be included in
     * Phase 1 (incorrectly increasing the het/hom calibration data), but countBothSideSupport
     * in Phase 2/3 would exclude it (since it uses {@code >}). The test verifies consistency
     * by checking that a variant with ONLY the borderline sample is not registered as a
     * training variant (firstPassCounts would be empty → finalizeFirstPass throws).</p>
     */
    @Test(expectedExceptions = IllegalStateException.class)
    public void testBothSideSupportThresholdIsStrictGreaterThan() {
        // 2 samples, 30x coverage (= targetCoverage, so normalization is identity).
        final List<String> samples = Arrays.asList("borderline", "noevidence");
        final SplitReadEvidenceGenotyper genotyper = new SplitReadEvidenceGenotyper(
                makeCoverageMap(samples, TARGET_COVERAGE),
                samples.size(), QUALITY_CUTOFF, MIN_SIZE, TARGET_COVERAGE, MAX_QUAL);

        final SVCallRecord record = makeDELRecord("del1", 1000, 5000,
                Arrays.asList(GATKSVVCFConstants.EvidenceTypes.PE, GATKSVVCFConstants.EvidenceTypes.SR),
                Collections.singletonList("pesr"));

        // "borderline" has raw count=1 on each side → normalized=Math.round(30*1/30)=1.0
        // With trainingCountCutoff=1 and minCount=max(1/2,1)=1, threshold=1.
        // Strict > means 1.0 > 1 is FALSE → sample excluded from both-side support.
        // Since "noevidence" has no SR, the variant has zero both-side support samples
        // and is NOT added to firstPassCounts.
        final List<SplitReadEvidence> startEvidence = makeSREvidence("borderline", "chr1", 1000, 1, true);
        final List<SplitReadEvidence> endEvidence = makeSREvidence("borderline", "chr1", 5000, 1, false);
        final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = makeDepthResult(new int[]{1, 2});

        genotyper.addFirstPass(record, startEvidence, endEvidence, depthResult, samples);

        // With no variants registered (borderline sample excluded by strict >),
        // firstPassCounts is empty → finalizeFirstPass should throw IllegalStateException.
        // If hasBothSideSupport incorrectly used >= , the borderline sample would pass,
        // the variant would be registered, and this would NOT throw.
        genotyper.finalizeFirstPass();
    }

    /**
     * Complementary test: a sample just above the threshold DOES qualify.
     */
    @Test
    public void testBothSideSupportAboveThresholdPasses() {
        final List<String> samples = Arrays.asList("abovethreshold", "noevidence");
        final SplitReadEvidenceGenotyper genotyper = new SplitReadEvidenceGenotyper(
                makeCoverageMap(samples, TARGET_COVERAGE),
                samples.size(), QUALITY_CUTOFF, MIN_SIZE, TARGET_COVERAGE, MAX_QUAL);

        final SVCallRecord record = makeDELRecord("del1", 1000, 5000,
                Arrays.asList(GATKSVVCFConstants.EvidenceTypes.PE, GATKSVVCFConstants.EvidenceTypes.SR),
                Collections.singletonList("pesr"));

        // Raw count=2 on each side → normalized=2.0 > threshold=1. Passes strict >.
        final List<SplitReadEvidence> startEvidence = makeSREvidence("abovethreshold", "chr1", 1000, 2, true);
        final List<SplitReadEvidence> endEvidence = makeSREvidence("abovethreshold", "chr1", 5000, 2, false);
        final DepthEvidenceGenotyper.DepthGenotypeResult depthResult = makeDepthResult(new int[]{1, 2});

        genotyper.addFirstPass(record, startEvidence, endEvidence, depthResult, samples);

        // Should not throw: the sample above threshold was registered
        genotyper.finalizeFirstPass();
    }

    // ---- Degenerate cutoff grid handling ----
    //
    // Selecting zero cutoffs is not a safe default: at application time
    // rare.freqMax == common.freqMin, so a zero cutoff turns the frequency predicate into a
    // tautology and disables SR background filtering. A grid that cannot distinguish a winner
    // must therefore be reported. Detection does not throw, because Cromwell delocalizes task
    // outputs only on success and a failure here would strand the diagnostics report on the
    // worker; enforcement is the ValidateSRCutoffs task in GenotypeBatch.wdl.

    private static final List<String> SIX_SAMPLES =
            Arrays.asList("s1", "s2", "s3", "s4", "s5", "s6");

    /** minCount 3.0 so that a normalized count of 5 is non-ref and 2 is not. */
    private static SplitReadEvidenceGenotyper.SplitReadGenotypeParameters gridTestParams() {
        return new SplitReadEvidenceGenotyper.SplitReadGenotypeParameters(3.0, 20.0, 5.0);
    }

    private static List<SplitReadEvidence> gridTestEvidence(final int position,
                                                            final List<String> nonRefSamples,
                                                            final List<String> backgroundSamples) {
        final List<SplitReadEvidence> evidence = new ArrayList<>();
        addSREvidence(evidence, nonRefSamples, "chr1", position, 5, true);
        addSREvidence(evidence, backgroundSamples, "chr1", position, 2, true);
        return evidence;
    }

    /**
     * With no entry satisfying the pass criteria the baseline is zero, so every cell scores
     * NaN and the objective is undefined. The cutoffs still fall back to (0.0, 0.0) as before,
     * but the rejection is now reported instead of being silent.
     */
    @Test
    public void testFinalizeThirdPassReportsGridWithNoPassingEntries() {
        final SplitReadEvidenceGenotyper genotyper = makeGenotyper(SIX_SAMPLES, NUM_SAMPLES);
        final SplitReadEvidenceGenotyper.SplitReadGenotypeParameters params = gridTestParams();

        // isCNV = false for every record, so pass is never true
        genotyper.accumulateHistogramOnly(makeINSRecord("fail_rare", 20000),
                gridTestEvidence(20000, Arrays.asList("s1", "s2"), Collections.emptyList()),
                Collections.emptyList(), false, params, SIX_SAMPLES);

        final SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs cutoffs =
                genotyper.finalizeThirdPass();

        Assert.assertEquals(cutoffs.rare().fracSingle(), 0.0, TEST_TOLERANCE);
        Assert.assertEquals(cutoffs.rare().fracBoth(), 0.0, TEST_TOLERANCE);
        assertRejectedWith(genotyper, SplitReadEvidenceGenotyper.SelectionStatus.REJECTED_NO_PASSING_ENTRIES);
        Assert.assertTrue(genotyper.cutoffDiagnosticsReport()
                        .contains("rare_selection_status\tREJECTED_NO_PASSING_ENTRIES"),
                "Diagnostics report must carry the machine-readable status that ValidateSRCutoffs greps");
    }

    /**
     * When every recovery fraction lands in the same histogram bin, no cutoff changes the
     * training partition, all 121 cells score identically, and any choice among them is
     * arbitrary. This is the collapse mode that produced all-zero cutoffs at scale.
     */
    @Test
    public void testFinalizeThirdPassReportsFullyTiedGrid() {
        final SplitReadEvidenceGenotyper genotyper = makeGenotyper(SIX_SAMPLES, NUM_SAMPLES);
        final SplitReadEvidenceGenotyper.SplitReadGenotypeParameters params = gridTestParams();
        final List<String> carriers = Arrays.asList("s1", "s2");

        // Both records put their frac at exactly 1.0 (2 non-ref of 2 samples over one),
        // so the pass and fail mass share a single bin.
        genotyper.accumulateHistogramOnly(
                makeDELRecord("pass_rare", 1000, 5000,
                        Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.SR),
                        Collections.singletonList("pesr")),
                gridTestEvidence(1000, carriers, Collections.emptyList()),
                Collections.emptyList(), true, params, SIX_SAMPLES);
        genotyper.accumulateHistogramOnly(makeINSRecord("fail_rare", 20000),
                gridTestEvidence(20000, carriers, Collections.emptyList()),
                Collections.emptyList(), false, params, SIX_SAMPLES);

        final SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs cutoffs =
                genotyper.finalizeThirdPass();

        Assert.assertEquals(cutoffs.rare().fracSingle(), 0.0, TEST_TOLERANCE);
        Assert.assertEquals(cutoffs.rare().fracBoth(), 0.0, TEST_TOLERANCE);
        assertRejectedWith(genotyper, SplitReadEvidenceGenotyper.SelectionStatus.REJECTED_DEGENERATE_GRID);
        Assert.assertTrue(genotyper.cutoffDiagnosticsReport()
                        .contains("rare_selection_status\tREJECTED_DEGENERATE_GRID"),
                "Diagnostics report must carry the machine-readable status that ValidateSRCutoffs greps");
    }

    private static void assertRejectedWith(final SplitReadEvidenceGenotyper genotyper,
                                           final SplitReadEvidenceGenotyper.SelectionStatus expected) {
        final List<SplitReadEvidenceGenotyper.SelectionOutcome> outcomes = genotyper.cutoffSelectionOutcomes();
        Assert.assertEquals(outcomes.size(), 2, "Expected one outcome per frequency bin");
        Assert.assertEquals(outcomes.get(0).status(), expected);
        Assert.assertTrue(outcomes.get(0).rejected());
        Assert.assertFalse(outcomes.get(0).detail().isEmpty(), "A rejection must explain itself");
    }

    /**
     * A partial tie is a legitimate optimum and must still be accepted, resolving to the
     * lowest tied cutoff. Here fail mass sits in bins 5 and 10 while pass mass sits only in
     * bin 10, so every fracSingle from 0.6 upward separates the training data equally well.
     * The lowest such cutoff is chosen because raising it further would reject ratios that
     * were never observed to be bad. This also guards against the degeneracy checks
     * over-reaching and rejecting a grid that does carry signal.
     */
    @Test
    public void testFinalizeThirdPassAcceptsPartialTieAndTakesLowestCutoff() {
        final SplitReadEvidenceGenotyper genotyper = makeGenotyper(SIX_SAMPLES, NUM_SAMPLES);
        final SplitReadEvidenceGenotyper.SplitReadGenotypeParameters params = gridTestParams();
        final List<String> rareCarriers = Arrays.asList("s1", "s2");
        final List<String> commonCarriers = Arrays.asList("s1", "s2", "s3");

        // Rare bin (non-ref count 2, at rareMax): pass at frac 1.0, fail at 1.0 and 0.5
        genotyper.accumulateHistogramOnly(
                makeDELRecord("pass_rare", 1000, 5000,
                        Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.SR),
                        Collections.singletonList("pesr")),
                gridTestEvidence(1000, rareCarriers, Collections.emptyList()),
                Collections.emptyList(), true, params, SIX_SAMPLES);
        genotyper.accumulateHistogramOnly(makeINSRecord("fail_rare_high", 20000),
                gridTestEvidence(20000, rareCarriers, Collections.emptyList()),
                Collections.emptyList(), false, params, SIX_SAMPLES);
        genotyper.accumulateHistogramOnly(makeINSRecord("fail_rare_low", 30000),
                gridTestEvidence(30000, rareCarriers, Arrays.asList("s3", "s4")),
                Collections.emptyList(), false, params, SIX_SAMPLES);

        // Common bin (non-ref count 3, above commonMin): pass at frac 1.0, fail at 0.5
        genotyper.accumulateHistogramOnly(
                makeDELRecord("pass_common", 40000, 45000,
                        Collections.singletonList(GATKSVVCFConstants.EvidenceTypes.SR),
                        Collections.singletonList("pesr")),
                gridTestEvidence(40000, commonCarriers, Collections.emptyList()),
                Collections.emptyList(), true, params, SIX_SAMPLES);
        genotyper.accumulateHistogramOnly(makeINSRecord("fail_common_low", 50000),
                gridTestEvidence(50000, commonCarriers, Arrays.asList("s4", "s5", "s6")),
                Collections.emptyList(), false, params, SIX_SAMPLES);

        final SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs cutoffs =
                genotyper.finalizeThirdPass();

        Assert.assertEquals(cutoffs.rare().fracSingle(), 0.6, TEST_TOLERANCE,
                "Rare tie should resolve to the lowest fracSingle that drops the fail mass in bin 5");
        Assert.assertEquals(cutoffs.rare().fracBoth(), 0.0, TEST_TOLERANCE,
                "With no both-sided data fracBoth is unconstrained and takes its lowest value");
        Assert.assertEquals(cutoffs.common().fracSingle(), 0.6, TEST_TOLERANCE,
                "Common tie should resolve to the lowest fracSingle that drops the fail mass in bin 5");
        Assert.assertEquals(cutoffs.common().fracBoth(), 0.0, TEST_TOLERANCE,
                "With no both-sided data fracBoth is unconstrained and takes its lowest value");
    }
}
