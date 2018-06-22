package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.metrics.MetricsFile;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;
import java.util.Map;

public class PathSeqScoreIntegrationTest extends CommandLineProgramTest {

    static final int SCORE_TABLE_COLUMNS = 10;
    static final double SCORE_ABSOLUTE_ERROR_TOLERANCE = 1e-6;

    @Override
    public String getTestedClassName() {
        return PathSeqScoreSpark.class.getSimpleName();
    }

    private static boolean testPathogenTaxonScores(final PSPathogenTaxonScore a, final PSPathogenTaxonScore b) {
        if (a.getKingdomTaxonId() != b.getKingdomTaxonId()) return false;
        if (!PathSeqTestUtils.equalWithinTolerance(a.getDescendentScore(),b.getDescendentScore(), SCORE_ABSOLUTE_ERROR_TOLERANCE)) return false;
        if (!PathSeqTestUtils.equalWithinTolerance(a.getTotalReads(),b.getTotalReads(), SCORE_ABSOLUTE_ERROR_TOLERANCE)) return false;
        if (!PathSeqTestUtils.equalWithinTolerance(a.getReferenceLength(),b.getReferenceLength(), SCORE_ABSOLUTE_ERROR_TOLERANCE)) return false;
        if (!PathSeqTestUtils.equalWithinTolerance(a.getSelfScore(),b.getSelfScore(), SCORE_ABSOLUTE_ERROR_TOLERANCE)) return false;
        if (!PathSeqTestUtils.equalWithinTolerance(a.getScoreNormalized(),b.getScoreNormalized(), SCORE_ABSOLUTE_ERROR_TOLERANCE)) return false;
        if (!PathSeqTestUtils.equalWithinTolerance(a.getUnambiguousReads(),b.getUnambiguousReads(), SCORE_ABSOLUTE_ERROR_TOLERANCE)) return false;
        return true;
    }

    private static Map<String, PSPathogenTaxonScore> getScores(final String[] lines) {
        final Map<String, PSPathogenTaxonScore> scores = new HashMap<>(lines.length);
        for (int i = 1; i < lines.length; i++) {
            String[] tok = lines[i].split("\t");
            Assert.assertTrue(tok.length == SCORE_TABLE_COLUMNS, "Expected " + SCORE_TABLE_COLUMNS + " columns, but found " + tok.length);
            final PSPathogenTaxonScore score = new PSPathogenTaxonScore();
            final String taxonId = tok[0];
            score.setKingdomTaxonId(Math.abs(tok[4].hashCode()));
            score.addSelfScore(new Double(tok[5]));
            score.addScoreNormalized(new Double(tok[6]));
            score.addTotalReads(new Integer(tok[7]));
            score.addUnambiguousReads(new Integer(tok[8]));
            score.setReferenceLength(new Long(tok[9]));
            Assert.assertFalse(scores.containsKey(taxonId), "Found more than one entry for taxon ID " + taxonId);
            scores.put(taxonId, score);
        }
        return scores;
    }

    public static void compareScoreTables(final String expected, final String test) {
        final String[] expectedLines = expected.split("\n");
        final String[] testLines = test.split("\n");
        if (expectedLines.length == 0 && testLines.length == 0) return;
        Assert.assertEquals(expectedLines[0], testLines[0], "Headers did not match");

        final Map<String, PSPathogenTaxonScore> expectedScores = getScores(expectedLines);
        final Map<String, PSPathogenTaxonScore> testScores = getScores(testLines);

        for (final String taxonId : testScores.keySet()) {
            Assert.assertTrue(expectedScores.containsKey(taxonId), "Did not expect the taxon ID: " + taxonId);
        }
        for (final String taxonId : expectedScores.keySet()) {
            Assert.assertTrue(testScores.containsKey(taxonId), "Did not find the expected taxon ID: " + taxonId);
        }
        for (final String taxonId : expectedScores.keySet()) {
            Assert.assertTrue(testPathogenTaxonScores(testScores.get(taxonId), expectedScores.get(taxonId)), "Scores for taxon ID " + taxonId + " did not match.");
        }
    }

    @DataProvider(name = "pathseqScoreTestData")
    public Object[][] getTestData() {
        return new Object[][]{
                {"expected_paired.txt", "expected_paired.metrics",
                        "alignment_paired.bam", null, false, false},
                {"expected_unpaired.txt", "expected_unpaired.metrics",
                        null, "alignment_unpaired.bam", false, false},
                {"expected_paired_genome_length.txt", "expected_paired_genome_length.metrics",
                        "alignment_paired.bam", null, true, false},
                {"expected_paired_kingdom_false.txt", "expected_paired_kingdom_false.metrics",
                        "alignment_paired.bam", null, false, true}
        };
    }
    @Test(dataProvider = "pathseqScoreTestData", groups = "spark")
    public void testPathSeqScoreSpark(final String expectedScoresFilename, final String expectedMetricsFilename,
                     final String inputPairedBamFilename, final String inputUnpairedBamFilename,
                     final boolean divideByGenomeLength, final boolean notNormalizedByKingdom) throws IOException {
        final File expectedScoresFile =  getTestFile(expectedScoresFilename);
        final File expectedMetricsFile = getTestFile(expectedMetricsFilename);
        final File inputPairedBamFile = inputPairedBamFilename == null ? null : getTestFile(inputPairedBamFilename);
        final File inputUnpairedBamFile = inputUnpairedBamFilename == null ? null : getTestFile(inputUnpairedBamFilename);
        final File taxFile = getTestFile("tax.db");
        final File outputScoresFile = createTempFile("test", ".txt");
        final File outputBamFile = createTempFile("output", ".bam");
        final File outputMetricsFile = createTempFile("score", ".metrics");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        if (inputPairedBamFile != null) {
            args.addFileArgument(PathSeqScoreSpark.PAIRED_INPUT_LONG_NAME, inputPairedBamFile);
        }
        if (inputUnpairedBamFile != null) {
            args.addFileArgument(PathSeqScoreSpark.UNPAIRED_INPUT_LONG_NAME, inputUnpairedBamFile);
        }
        args.addFileArgument(PSScoreArgumentCollection.SCORE_METRICS_FILE_LONG_NAME, outputMetricsFile);
        args.addFileArgument(PSScoreArgumentCollection.TAXONOMIC_DATABASE_LONG_NAME, taxFile);
        args.addFileArgument(PSScoreArgumentCollection.SCORES_OUTPUT_LONG_NAME, outputScoresFile);
        args.addOutput(outputBamFile);
        args.addBooleanArgument(PSScoreArgumentCollection.DIVIDE_BY_GENOME_LENGTH_LONG_NAME, divideByGenomeLength);
        args.addBooleanArgument(PSScoreArgumentCollection.NOT_NORMALIZED_BY_KINGDOM_LONG_NAME, notNormalizedByKingdom);

        this.runCommandLine(args.getArgsArray());

        final String expectedScoresString = FileUtils.readFileToString(expectedScoresFile, StandardCharsets.UTF_8);
        final String actualScoresString = FileUtils.readFileToString(outputScoresFile, StandardCharsets.UTF_8);
        compareScoreTables(expectedScoresString, actualScoresString);

        Assert.assertTrue(MetricsFile.areMetricsEqual(outputMetricsFile, expectedMetricsFile));
    }

}