package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Unit tests for {@link DepthEvidenceGenotyper}.
 */
public class DepthEvidenceGenotyperTest extends GATKBaseTest {

    private static final SAMSequenceDictionary DICTIONARY = SVTestUtils.hg38Dict;
    private static final double TEST_TOLERANCE = 1e-4;
    private static final int MAX_QUAL = 999;

    // Default constants from the source
    private static final double DEFAULT_COPY_STATE_INCREMENT = 0.5;


    private static final List<String> TWO_SAMPLES = Arrays.asList("sample1", "sample2");
    private static final List<String> THREE_SAMPLES = Arrays.asList("sample1", "sample2", "sample3");
    private static final List<String> FIVE_SAMPLES = Arrays.asList("s1", "s2", "s3", "s4", "s5");

    // ---- Helper methods ----

    private static DepthMatrix makeMatrix(final List<String> samples, final double[] sampleMedianValues) {
        // Create a single-bin matrix. With one bin the median equals the value itself.
        final List<SimpleInterval> bins = DepthMatrixTest.makeBins("chr1", 1, 1001, 1000);
        final Map<String, double[]> matrix = new HashMap<>();
        for (int i = 0; i < samples.size(); i++) {
            matrix.put(samples.get(i), new double[]{sampleMedianValues[i]});
        }
        return new DepthMatrix(bins, matrix);
    }

    private static DepthMatrix makeMultiBinMatrix(final List<String> samples, final double[][] sampleValues) {
        final int numBins = sampleValues[0].length;
        final List<SimpleInterval> bins = DepthMatrixTest.makeBins("chr1", 1, 1 + numBins * 1000, 1000);
        final Map<String, double[]> matrix = new HashMap<>();
        for (int i = 0; i < samples.size(); i++) {
            matrix.put(samples.get(i), sampleValues[i]);
        }
        return new DepthMatrix(bins, matrix);
    }

    private static DepthEvidenceGenotyper makeDefaultGenotyper(final List<String> samples) {
        return new DepthEvidenceGenotyper(null, samples, MAX_QUAL, DICTIONARY);
    }

    private static DepthEvidenceGenotyper makeGenotyperWithCutoffs(final List<DepthEvidenceGenotyper.CopyStateStats> cutoffs,
                                                                    final List<String> samples) {
        return new DepthEvidenceGenotyper(cutoffs, samples, MAX_QUAL, DICTIONARY);
    }

    // ---- Constructor tests ----

    @Test
    public void testDefaultConstructor() {
        final DepthEvidenceGenotyper genotyper = makeDefaultGenotyper(TWO_SAMPLES);
        Assert.assertNotNull(genotyper);
        Assert.assertEquals(genotyper.getSamplesInOrder(), TWO_SAMPLES);
    }

    @Test
    public void testCustomCutoffsConstructor() {
        final List<DepthEvidenceGenotyper.CopyStateStats> cutoffs = Arrays.asList(
                new DepthEvidenceGenotyper.CopyStateStats(0, 0.0, 0.10, 0.20),
                new DepthEvidenceGenotyper.CopyStateStats(1, 0.45, 0.12, 0.70),
                new DepthEvidenceGenotyper.CopyStateStats(2, 1.0, 0.14, 1.30)
        );
        final DepthEvidenceGenotyper genotyper = makeGenotyperWithCutoffs(cutoffs, TWO_SAMPLES);
        Assert.assertNotNull(genotyper);
        Assert.assertEquals(genotyper.getSamplesInOrder(), TWO_SAMPLES);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testConstructorNullDictionary() {
        new DepthEvidenceGenotyper(null, TWO_SAMPLES, MAX_QUAL, null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testConstructorNullSamples() {
        new DepthEvidenceGenotyper(null, null, MAX_QUAL, DICTIONARY);
    }

    // ---- genotype() tests ----

    @Test
    public void testGenotypeHomRef() {
        // All samples at copy number ~2 (median depth = 1.0 in normalized space)
        final DepthEvidenceGenotyper genotyper = makeDefaultGenotyper(TWO_SAMPLES);
        final DepthMatrix matrix = makeMatrix(TWO_SAMPLES, new double[]{1.0, 1.0});
        final DepthEvidenceGenotyper.DepthGenotypeResult result = genotyper.genotype(matrix);

        Assert.assertNotNull(result);
        Assert.assertEquals(result.copyStates().length, 2);
        // Depth of 1.0 maps to copy state 2 with default cutoffs (0, 0.25, 0.75, 1.25, 1.75, 2.25)
        Assert.assertEquals(result.copyStates()[0], 2);
        Assert.assertEquals(result.copyStates()[1], 2);
        // Variant quality should be 0 when all samples are hom-ref
        Assert.assertEquals(result.variantQual(), 0., TEST_TOLERANCE);
    }

    @Test
    public void testGenotypeHetDel() {
        // sample1 at CN=1 (depth ~0.5), sample2 at CN=2 (depth ~1.0)
        final DepthEvidenceGenotyper genotyper = makeDefaultGenotyper(TWO_SAMPLES);
        final DepthMatrix matrix = makeMatrix(TWO_SAMPLES, new double[]{0.5, 1.0});
        final DepthEvidenceGenotyper.DepthGenotypeResult result = genotyper.genotype(matrix);

        Assert.assertNotNull(result);
        Assert.assertEquals(result.copyStates()[0], 1); // het del
        Assert.assertEquals(result.copyStates()[1], 2); // hom ref
        // Variant quality should be > 0 because we have a carrier
        Assert.assertTrue(result.variantQual() > 0, "Variant quality should be > 0 for variant with carrier");
    }

    @Test
    public void testGenotypeHomDel() {
        // sample1 at CN=0 (depth ~0.0), sample2 at CN=2
        final DepthEvidenceGenotyper genotyper = makeDefaultGenotyper(TWO_SAMPLES);
        final DepthMatrix matrix = makeMatrix(TWO_SAMPLES, new double[]{0.0, 1.0});
        final DepthEvidenceGenotyper.DepthGenotypeResult result = genotyper.genotype(matrix);

        Assert.assertNotNull(result);
        Assert.assertEquals(result.copyStates()[0], 0); // hom del
        Assert.assertEquals(result.copyStates()[1], 2); // hom ref
        Assert.assertTrue(result.variantQual() > 0);
    }

    @Test
    public void testGenotypeDup() {
        // sample1 at CN=3 (depth ~1.5), sample2 at CN=2
        final DepthEvidenceGenotyper genotyper = makeDefaultGenotyper(TWO_SAMPLES);
        final DepthMatrix matrix = makeMatrix(TWO_SAMPLES, new double[]{1.5, 1.0});
        final DepthEvidenceGenotyper.DepthGenotypeResult result = genotyper.genotype(matrix);

        Assert.assertNotNull(result);
        Assert.assertEquals(result.copyStates()[0], 3); // het dup
        Assert.assertEquals(result.copyStates()[1], 2); // hom ref
        Assert.assertTrue(result.variantQual() > 0);
    }

    @Test
    public void testGenotypeHomDup() {
        // sample1 at CN=4 (depth ~2.0), sample2 at CN=2
        final DepthEvidenceGenotyper genotyper = makeDefaultGenotyper(TWO_SAMPLES);
        final DepthMatrix matrix = makeMatrix(TWO_SAMPLES, new double[]{2.0, 1.0});
        final DepthEvidenceGenotyper.DepthGenotypeResult result = genotyper.genotype(matrix);

        Assert.assertNotNull(result);
        Assert.assertEquals(result.copyStates()[0], 4);
        Assert.assertEquals(result.copyStates()[1], 2);
    }

    @Test
    public void testGenotypeHighCopyNumber() {
        // Value beyond the default 6 states (CN>5)
        final DepthEvidenceGenotyper genotyper = makeDefaultGenotyper(TWO_SAMPLES);
        // depth = 4.0 should map to CN=8 with default increment 0.5
        final DepthMatrix matrix = makeMatrix(TWO_SAMPLES, new double[]{4.0, 1.0});
        final DepthEvidenceGenotyper.DepthGenotypeResult result = genotyper.genotype(matrix);

        Assert.assertNotNull(result);
        Assert.assertTrue(result.copyStates()[0] > 5, "High depth should map to high copy state");
        Assert.assertEquals(result.copyStates()[1], 2);
    }

    @Test
    public void testGenotypeEmptyBins() {
        // Edge case: no bins
        final DepthEvidenceGenotyper genotyper = makeDefaultGenotyper(TWO_SAMPLES);
        final DepthMatrix matrix = new DepthMatrix(
                Collections.emptyList(),
                Map.of("sample1", new double[]{}, "sample2", new double[]{})
        );
        final DepthEvidenceGenotyper.DepthGenotypeResult result = genotyper.genotype(matrix);

        Assert.assertNotNull(result);
        // All samples should get NO_DATA values
        Assert.assertEquals(result.copyStates()[0], 2); // NO_DATA_COPY_STATE
        Assert.assertEquals(result.copyStates()[1], 2);
        Assert.assertEquals(result.genotypeQuals()[0], 2., TEST_TOLERANCE); // NO_DATA_GENOTYPE_QUAL
        Assert.assertEquals(result.genotypeQuals()[1], 2., TEST_TOLERANCE);
        Assert.assertEquals(result.variantQual(), 0., TEST_TOLERANCE);
    }

    @Test
    public void testGenotypeSampleDepthsReturned() {
        // Verify that sampleDepths in the result contains the correct medians
        final DepthEvidenceGenotyper genotyper = makeDefaultGenotyper(TWO_SAMPLES);
        final DepthMatrix matrix = makeMatrix(TWO_SAMPLES, new double[]{0.5, 1.5});
        final DepthEvidenceGenotyper.DepthGenotypeResult result = genotyper.genotype(matrix);

        Assert.assertNotNull(result);
        Assert.assertEquals(result.sampleDepths().length, 2);
        // Single-bin matrix, so medians equal the input values
        Assert.assertEquals(result.sampleDepths()[0], 0.5, TEST_TOLERANCE);
        Assert.assertEquals(result.sampleDepths()[1], 1.5, TEST_TOLERANCE);
    }

    @Test
    public void testGenotypeQualitiesPositive() {
        // All genotype qualities should be >= 0
        final DepthEvidenceGenotyper genotyper = makeDefaultGenotyper(THREE_SAMPLES);
        final DepthMatrix matrix = makeMatrix(THREE_SAMPLES, new double[]{0.0, 0.5, 1.0});
        final DepthEvidenceGenotyper.DepthGenotypeResult result = genotyper.genotype(matrix);

        for (final double gq : result.genotypeQuals()) {
            Assert.assertTrue(gq >= 0, "Genotype quality should be >= 0");
        }
    }

    @Test
    public void testGenotypeQualitiesBoundedByMax() {
        // Test with a low maxQual
        final int lowMaxQual = 50;
        final DepthEvidenceGenotyper genotyper = new DepthEvidenceGenotyper(null, TWO_SAMPLES, lowMaxQual, DICTIONARY);
        final DepthMatrix matrix = makeMatrix(TWO_SAMPLES, new double[]{1.0, 1.0});
        final DepthEvidenceGenotyper.DepthGenotypeResult result = genotyper.genotype(matrix);

        for (final double gq : result.genotypeQuals()) {
            Assert.assertTrue(gq <= lowMaxQual, "Genotype quality should be bounded by maxQual");
        }
    }

    @Test
    public void testGenotypeMultipleBins() {
        // Test with multiple bins where median is more meaningful
        final DepthEvidenceGenotyper genotyper = makeDefaultGenotyper(TWO_SAMPLES);
        final DepthMatrix matrix = makeMultiBinMatrix(TWO_SAMPLES, new double[][]{
                {0.4, 0.5, 0.6},   // sample1: median ~0.5 → CN=1
                {0.9, 1.0, 1.1}    // sample2: median ~1.0 → CN=2
        });
        final DepthEvidenceGenotyper.DepthGenotypeResult result = genotyper.genotype(matrix);

        Assert.assertNotNull(result);
        Assert.assertEquals(result.copyStates()[0], 1);
        Assert.assertEquals(result.copyStates()[1], 2);
    }

    @Test
    public void testGenotypeWithCustomCutoffs() {
        // Custom cutoffs with different boundaries
        final List<DepthEvidenceGenotyper.CopyStateStats> cutoffs = Arrays.asList(
                new DepthEvidenceGenotyper.CopyStateStats(0, 0.0, 0.10, 0.30),
                new DepthEvidenceGenotyper.CopyStateStats(1, 0.5, 0.12, 0.80),
                new DepthEvidenceGenotyper.CopyStateStats(2, 1.0, 0.14, 1.40),
                new DepthEvidenceGenotyper.CopyStateStats(3, 1.5, 0.14, 1.90)
        );
        final DepthEvidenceGenotyper genotyper = makeGenotyperWithCutoffs(cutoffs, TWO_SAMPLES);
        // Value 0.75 should be in state 1 (upperBound for state 1 = 0.80)
        final DepthMatrix matrix = makeMatrix(TWO_SAMPLES, new double[]{0.75, 1.0});
        final DepthEvidenceGenotyper.DepthGenotypeResult result = genotyper.genotype(matrix);

        Assert.assertEquals(result.copyStates()[0], 1);
        Assert.assertEquals(result.copyStates()[1], 2);
    }

    // ---- train() tests ----

    @Test
    public void testTrainProducesNonDefaultValues() {
        // Generate some genotype results with known copy states and depths
        final DepthEvidenceGenotyper genotyper = makeDefaultGenotyper(FIVE_SAMPLES);
        final List<DepthEvidenceGenotyper.DepthGenotypeResult> results = new ArrayList<>();

        // Create results where:
        // - Some samples have CN=0 with depth ~0
        // - Some have CN=1 with depth ~0.5
        // - Some have CN=2 with depth ~1.0
        // - Some have CN=3 with depth ~1.5
        results.add(new DepthEvidenceGenotyper.DepthGenotypeResult(
                new double[]{0.02, 0.48, 1.0, 1.52, 0.99},
                new int[]{0, 1, 2, 3, 2},
                new double[]{30, 25, 30, 25, 30},
                10
        ));
        results.add(new DepthEvidenceGenotyper.DepthGenotypeResult(
                new double[]{0.01, 0.51, 1.01, 1.49, 0.50},
                new int[]{0, 1, 2, 3, 1},
                new double[]{30, 25, 30, 25, 25},
                10
        ));
        results.add(new DepthEvidenceGenotyper.DepthGenotypeResult(
                new double[]{0.03, 0.49, 0.98, 1.51, 1.48},
                new int[]{0, 1, 2, 3, 3},
                new double[]{30, 25, 30, 25, 25},
                10
        ));

        final int numStates = 4;
        final List<DepthEvidenceGenotyper.CopyStateStats> trained = genotyper.train(results, numStates);

        Assert.assertNotNull(trained);
        Assert.assertEquals(trained.size(), numStates);

        // Verify copy states are in order
        for (int i = 0; i < numStates; i++) {
            Assert.assertEquals(trained.get(i).copyState(), i);
        }

        // Trained means should be close to the actual depth distributions
        Assert.assertEquals(trained.get(0).mean(), 0.02, 0.05);   // CN=0 cluster around 0
        Assert.assertEquals(trained.get(1).mean(), 0.496, 0.05);  // CN=1 cluster around 0.5
        Assert.assertEquals(trained.get(2).mean(), 0.9967, 0.05); // CN=2 cluster around 1.0
        Assert.assertEquals(trained.get(3).mean(), 1.5, 0.05);    // CN=3 cluster around 1.5

        // Std devs should be > 0 for states with observations
        for (int i = 0; i < numStates; i++) {
            Assert.assertTrue(trained.get(i).stdDev() >= 0, "Std dev should be >= 0");
        }

        // Upper bounds should be monotonically increasing
        for (int i = 0; i < numStates - 1; i++) {
            Assert.assertTrue(trained.get(i).upperBound() < trained.get(i + 1).upperBound(),
                    "Upper bounds should be monotonically increasing");
        }
    }

    @Test
    public void testTrainWithMissingState() {
        // Train with no observations for state 1 — should fall back to defaults
        final DepthEvidenceGenotyper genotyper = makeDefaultGenotyper(TWO_SAMPLES);
        final List<DepthEvidenceGenotyper.DepthGenotypeResult> results = new ArrayList<>();
        results.add(new DepthEvidenceGenotyper.DepthGenotypeResult(
                new double[]{0.0, 1.0},
                new int[]{0, 2},
                new double[]{30, 30},
                10
        ));
        results.add(new DepthEvidenceGenotyper.DepthGenotypeResult(
                new double[]{0.01, 0.99},
                new int[]{0, 2},
                new double[]{30, 30},
                10
        ));

        final int numStates = 3;
        final List<DepthEvidenceGenotyper.CopyStateStats> trained = genotyper.train(results, numStates);
        Assert.assertEquals(trained.size(), numStates);

        // State 1 has no data, should fall back to default mean = 0.5
        Assert.assertEquals(trained.get(1).mean(), DEFAULT_COPY_STATE_INCREMENT, TEST_TOLERANCE);
    }

    @Test
    public void testTrainRoundTrip() {
        // Train, then create a new genotyper with trained cutoffs, and verify it produces
        // correct copy states
        final DepthEvidenceGenotyper genotyper = makeDefaultGenotyper(THREE_SAMPLES);
        final List<DepthEvidenceGenotyper.DepthGenotypeResult> results = new ArrayList<>();
        for (int i = 0; i < 10; i++) {
            results.add(new DepthEvidenceGenotyper.DepthGenotypeResult(
                    new double[]{0.01 + i * 0.001, 0.50 + i * 0.001, 1.0 + i * 0.001},
                    new int[]{0, 1, 2},
                    new double[]{30, 25, 30},
                    10
            ));
        }

        final int numStates = 4;
        final List<DepthEvidenceGenotyper.CopyStateStats> trained = genotyper.train(results, numStates);

        // Build new genotyper with trained cutoffs
        final DepthEvidenceGenotyper trainedGenotyper = makeGenotyperWithCutoffs(trained, THREE_SAMPLES);

        // Now genotype with the trained genotyper — CN=0, CN=1, CN=2
        final DepthMatrix matrix = makeMatrix(THREE_SAMPLES, new double[]{0.01, 0.50, 1.0});
        final DepthEvidenceGenotyper.DepthGenotypeResult result = trainedGenotyper.genotype(matrix);

        Assert.assertEquals(result.copyStates()[0], 0);
        Assert.assertEquals(result.copyStates()[1], 1);
        Assert.assertEquals(result.copyStates()[2], 2);
    }

    // ---- CopyStateStats record tests ----

    @Test
    public void testCopyStateStatsRecord() {
        final DepthEvidenceGenotyper.CopyStateStats stats = new DepthEvidenceGenotyper.CopyStateStats(2, 1.0, 0.15, 1.25);
        Assert.assertEquals(stats.copyState(), 2);
        Assert.assertEquals(stats.mean(), 1.0, TEST_TOLERANCE);
        Assert.assertEquals(stats.stdDev(), 0.15, TEST_TOLERANCE);
        Assert.assertEquals(stats.upperBound(), 1.25, TEST_TOLERANCE);
    }

    // ---- DepthGenotypeResult record tests ----

    @Test
    public void testDepthGenotypeResultRecord() {
        final double[] depths = {0.0, 0.5, 1.0};
        final int[] states = {0, 1, 2};
        final double[] quals = {30, 25, 30};
        final DepthEvidenceGenotyper.DepthGenotypeResult result =
                new DepthEvidenceGenotyper.DepthGenotypeResult(depths, states, quals, 15.0);

        Assert.assertEquals(result.sampleDepths(), depths);
        Assert.assertEquals(result.copyStates(), states);
        Assert.assertEquals(result.genotypeQuals(), quals);
        Assert.assertEquals(result.variantQual(), 15.0, TEST_TOLERANCE);
    }

    // ---- getCutoff (static) ----

    @DataProvider(name = "cutoffData")
    public Object[][] cutoffData() {
        return new Object[][]{
                // Equal std devs → midpoint
                {0.0, 0.1, 1.0, 0.1, 0.5},
                // Different std devs → weighted toward the narrower distribution
                {0.0, 0.1, 1.0, 0.2, 1.0 / 3.0},
                // Very different → strongly weighted
                {0.0, 0.01, 1.0, 1.0, 0.0099, /* tolerance = */ 0.01},
        };
    }

    @Test
    public void testCutoffViaTraining() {
        // The getCutoff method is private, but we can verify its behavior via train()
        // When means are ~0.0 and ~1.0 with similar std devs, the cutoff between them should be ~0.5
        final DepthEvidenceGenotyper genotyper = makeDefaultGenotyper(TWO_SAMPLES);
        final List<DepthEvidenceGenotyper.DepthGenotypeResult> results = new ArrayList<>();
        final Random rng = new Random(42);
        for (int i = 0; i < 50; i++) {
            // Add small variance so std devs are non-zero
            final double noise0 = rng.nextGaussian() * 0.05;
            final double noise1 = rng.nextGaussian() * 0.05;
            results.add(new DepthEvidenceGenotyper.DepthGenotypeResult(
                    new double[]{0.0 + noise0, 1.0 + noise1},
                    new int[]{0, 1},
                    new double[]{30, 30},
                    10
            ));
        }
        final List<DepthEvidenceGenotyper.CopyStateStats> trained = genotyper.train(results, 2);
        // With similar std devs and means at ~0.0 and ~1.0, cutoff should be ~0.5
        Assert.assertEquals(trained.get(0).upperBound(), 0.5, 0.1);
    }

    // ---- Variant quality tests ----

    @Test
    public void testVariantQualAllRef() {
        final DepthEvidenceGenotyper genotyper = makeDefaultGenotyper(THREE_SAMPLES);
        final DepthMatrix matrix = makeMatrix(THREE_SAMPLES, new double[]{1.0, 1.0, 1.0});
        final DepthEvidenceGenotyper.DepthGenotypeResult result = genotyper.genotype(matrix);
        Assert.assertEquals(result.variantQual(), 0., TEST_TOLERANCE,
                "Variant quality should be 0 when all samples are reference");
    }

    @Test
    public void testVariantQualWithCarriers() {
        final DepthEvidenceGenotyper genotyper = makeDefaultGenotyper(THREE_SAMPLES);
        // Two carriers (CN=0 and CN=1) plus one ref
        final DepthMatrix matrix = makeMatrix(THREE_SAMPLES, new double[]{0.0, 0.5, 1.0});
        final DepthEvidenceGenotyper.DepthGenotypeResult result = genotyper.genotype(matrix);
        Assert.assertTrue(result.variantQual() > 0,
                "Variant quality should be > 0 with carriers");
    }

    @Test
    public void testVariantQualBoundedByMaxQual() {
        final int lowMaxQual = 30;
        final DepthEvidenceGenotyper genotyper = new DepthEvidenceGenotyper(null, TWO_SAMPLES, lowMaxQual, DICTIONARY);
        final DepthMatrix matrix = makeMatrix(TWO_SAMPLES, new double[]{0.0, 1.0});
        final DepthEvidenceGenotyper.DepthGenotypeResult result = genotyper.genotype(matrix);
        Assert.assertTrue(result.variantQual() <= lowMaxQual,
                "Variant quality should be bounded by maxQual");
    }

    // ---- DepthTableParser tests ----

    @Test
    public void testDepthTableParserColumns() {
        Assert.assertNotNull(DepthEvidenceGenotyper.DepthTableParser.CUTOFFS_COLUMNS);
        Assert.assertEquals(DepthEvidenceGenotyper.DepthTableParser.CUTOFFS_COLUMNS.columnCount(), 4);
    }

    // ---- Edge case: all samples same depth ----

    @Test
    public void testGenotypeAllSameDepth() {
        // All samples at the same depth
        final DepthEvidenceGenotyper genotyper = makeDefaultGenotyper(THREE_SAMPLES);
        final DepthMatrix matrix = makeMatrix(THREE_SAMPLES, new double[]{0.5, 0.5, 0.5});
        final DepthEvidenceGenotyper.DepthGenotypeResult result = genotyper.genotype(matrix);

        Assert.assertNotNull(result);
        // All should be CN=1
        for (final int cs : result.copyStates()) {
            Assert.assertEquals(cs, 1);
        }
        // With all carriers (CN=1), variant qual should be > 0
        Assert.assertTrue(result.variantQual() > 0);
    }

    // ---- Default copy state boundaries ----

    @DataProvider(name = "defaultCopyStateData")
    public Object[][] defaultCopyStateData() {
        // Test that depths map to the expected copy states with default cutoffs
        // Default cutoffs: state i has upperBound = (i + 0.5) * 0.5
        // State 0: upperBound = 0.25
        // State 1: upperBound = 0.75
        // State 2: upperBound = 1.25
        // State 3: upperBound = 1.75
        // State 4: upperBound = 2.25
        // State 5: upperBound = 2.75
        return new Object[][]{
                {0.0, 0},
                {0.24, 0},
                {0.26, 1},
                {0.74, 1},
                {0.76, 2},
                {1.24, 2},
                {1.26, 3},
                {1.74, 3},
                {1.76, 4},
                {2.24, 4},
                {2.26, 5},
                {2.74, 5},
        };
    }

    @Test(dataProvider = "defaultCopyStateData")
    public void testDefaultCopyStateBoundaries(final double depth, final int expectedCopyState) {
        final List<String> samples = Collections.singletonList("sample");
        final DepthEvidenceGenotyper genotyper = makeDefaultGenotyper(samples);
        final DepthMatrix matrix = makeMatrix(samples, new double[]{depth});
        final DepthEvidenceGenotyper.DepthGenotypeResult result = genotyper.genotype(matrix);
        Assert.assertEquals(result.copyStates()[0], expectedCopyState,
                "Depth " + depth + " should map to copy state " + expectedCopyState);
    }
}
