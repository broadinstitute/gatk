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

    @Test(groups = "spark")
    public void testDivideByGenomeLength() throws IOException {
        final File expectedFile = getTestFile("expected_paired_genome_length.txt");
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
        args.addBooleanArgument("divideByGenomeLength", true);
        this.runCommandLine(args.getArgsArray());

        final String input_expected = FileUtils.readFileToString(expectedFile);
        final String input_test = FileUtils.readFileToString(output);
        compareScoreTables(input_expected, input_test);
    }

    @Test(groups = "spark")
    public void testKingdomNormalizationFalse() throws IOException {
        final File expectedFile = getTestFile("expected_paired_kingdom_false.txt");
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
        args.addBooleanArgument("notNormalizedByKingdom", true);
        this.runCommandLine(args.getArgsArray());

        final String input_expected = FileUtils.readFileToString(expectedFile);
        final String input_test = FileUtils.readFileToString(output);
        compareScoreTables(input_expected, input_test);
    }

}