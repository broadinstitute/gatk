package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class PathSeqScoreIntegrationTest extends CommandLineProgramTest {

    final static int SCORE_TABLE_COLUMNS = 9;
    final static double SCORE_ABSOLUTE_ERROR_TOLERANCE = 1e-6;

    @Override
    public String getTestedClassName() {
        return PathSeqScoreSpark.class.getSimpleName();
    }

    private static boolean testPathogenTaxonScores(final PSPathogenTaxonScore a, final PSPathogenTaxonScore b) {
        if (!PathSeqTestUtils.equalWithinTolerance(a.descendentScore,b.descendentScore, SCORE_ABSOLUTE_ERROR_TOLERANCE)) return false;
        if (!PathSeqTestUtils.equalWithinTolerance(a.totalReads,b.totalReads, SCORE_ABSOLUTE_ERROR_TOLERANCE)) return false;
        if (!PathSeqTestUtils.equalWithinTolerance(a.referenceLength,b.referenceLength, SCORE_ABSOLUTE_ERROR_TOLERANCE)) return false;
        if (!PathSeqTestUtils.equalWithinTolerance(a.selfScore,b.selfScore, SCORE_ABSOLUTE_ERROR_TOLERANCE)) return false;
        if (!PathSeqTestUtils.equalWithinTolerance(a.scoreNormalized,b.scoreNormalized, SCORE_ABSOLUTE_ERROR_TOLERANCE)) return false;
        if (!PathSeqTestUtils.equalWithinTolerance(a.unambiguousReads,b.unambiguousReads, SCORE_ABSOLUTE_ERROR_TOLERANCE)) return false;
        return true;
    }

    private static Map<String, PSPathogenTaxonScore> getScores(final String[] lines) {
        final Map<String, PSPathogenTaxonScore> scores = new HashMap<>(lines.length);
        for (int i = 1; i < lines.length; i++) {
            String[] tok = lines[i].split("\t");
            Assert.assertTrue(tok.length == SCORE_TABLE_COLUMNS, "Expected " + SCORE_TABLE_COLUMNS + " columns, but found " + tok.length);
            final PSPathogenTaxonScore score = new PSPathogenTaxonScore();
            final String taxonId = tok[0];
            score.selfScore = new Double(tok[4]);
            score.scoreNormalized = new Double(tok[5]);
            score.totalReads = new Integer(tok[6]);
            score.unambiguousReads = new Integer(tok[7]);
            score.referenceLength = new Long(tok[8]);
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

        for (final String taxonId : expectedScores.keySet()) {
            Assert.assertTrue(testScores.containsKey(taxonId), "Did not find taxon ID: " + taxonId);
            Assert.assertTrue(testPathogenTaxonScores(testScores.get(taxonId), expectedScores.get(taxonId)), "Scores for taxon ID " + taxonId + " did not match.");
        }
    }

    @Test(groups = "spark")
    public void test() throws IOException {
        final File expectedFile = getTestFile("expected_paired.txt");
        final File inputFile = getTestFile("alignment_paired.bam");
        final File taxFile = getTestFile("tax.db");
        final File output = createTempFile("test", ".txt");
        final File bamOut = createTempFile("output", ".bam");
        if (!output.delete() || !bamOut.delete()) {
            Assert.fail();
        }
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addArgument("pairedInput", inputFile.getAbsolutePath());
        args.addOutput(bamOut);
        args.addFileArgument("taxonomicDatabasePath", taxFile);
        args.addFileArgument("scoresOutputPath", output);
        this.runCommandLine(args.getArgsArray());

        final String input_expected = FileUtils.readFileToString(expectedFile);
        final String input_test = FileUtils.readFileToString(output);
        compareScoreTables(input_expected, input_test);
    }

    @Test(groups = "spark")
    public void testUnpaired() throws IOException {
        final File expectedFile = getTestFile("expected_unpaired.txt");
        final File inputFile = getTestFile("alignment_unpaired.bam");
        final File taxFile = getTestFile("tax.db");
        final File output = createTempFile("test", ".txt");
        final File warnings = createTempFile("warnings", ".txt");
        if (!output.delete() || !warnings.delete()) {
            Assert.fail();
        }
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addArgument("unpairedInput", inputFile.getAbsolutePath());
        args.addFileArgument("taxonomicDatabasePath", taxFile);
        args.addFileArgument("scoresOutputPath", output);
        args.addFileArgument("scoreWarningsFile",warnings);
        this.runCommandLine(args.getArgsArray());

        final String input_expected = FileUtils.readFileToString(expectedFile);
        final String input_test = FileUtils.readFileToString(output);
        compareScoreTables(input_expected, input_test);
    }

}