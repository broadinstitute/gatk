package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

/**
 * Integration tests for {@link CalculateTargetCoverage}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class CalculateTargetCoverageIntegrationTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "exome");

    private static final File TEST_BAM_NA12872 = new File(TEST_DATA_DIR, "exome-read-counts-NA12872.bam");
    private static final File TEST_BAM_NA12778 = new File(TEST_DATA_DIR, "exome-read-counts-NA12778.bam");
    private static final File TEST_BAM_NA12878 = new File(TEST_DATA_DIR, "exome-read-counts-NA12878.bam");
    private static final File TARGETS_FILE = new File(TEST_DATA_DIR, "exome-read-counts-test-targets.tsv");
    private static final File EXPECTED_COUNTS_FILE = new File(TEST_DATA_DIR, "exome-read-counts.output");
    private static final File EXPECTED_COUNTS_MAX_OF_9_FILE = new File(TEST_DATA_DIR, "exome-read-counts-max-of-9.output");
    private static final File EXPECTED_COUNTS_MIN_MQ_30_FILE = new File(TEST_DATA_DIR, "exome-read-counts-min-MQ-30.output");
    private static final File TARGETS_FILE_WITHOUT_COORDINATES = new File(TEST_DATA_DIR, "exome-read-counts-test-targets-wo-coords.tsv");
    private static final File INEXISTENT_TARGETS_FILE = new File(TEST_DATA_DIR, "fantasy-exome-read-counts-test-targets.tsv");

    @Override
    public String getTestedClassName() {
        return CalculateTargetCoverage.class.getSimpleName();
    }

    @Test
    public void testSimpleRun() throws IOException {
        final File outputFile = createTempFile("ctc-test-",".tsv");
        runCommandLine(
                new String[]{
                        "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME,
                        TARGETS_FILE.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12878.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12778.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12872.toString(),
                        "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
                        outputFile.toString()
                });
        Assert.assertTrue(outputFile.canRead());
        final ReadCountCollection outputCounts = ReadCountCollectionUtils.parse(outputFile);
        final ReadCountCollection expectedCounts = ReadCountCollectionUtils.parse(EXPECTED_COUNTS_FILE);

        Assert.assertEquals(outputCounts.columnNames(), expectedCounts.columnNames());
        Assert.assertEquals(outputCounts.counts().getData(), expectedCounts.counts().getData());
    }

    @Test
    public void testSimpleRunWithMaximumCoverage() throws IOException {
        final File outputFile = createTempFile("ctc-test-",".tsv");
        runCommandLine(
                new String[]{
                        "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME,
                        TARGETS_FILE.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12878.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12778.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12872.toString(),
                        "-" + CalculateTargetCoverage.MAXIMUM_COVERAGE_SHORT_NAME,
                        "9",
                        "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
                        outputFile.toString()
                });
        Assert.assertTrue(outputFile.canRead());
        final ReadCountCollection outputCounts = ReadCountCollectionUtils.parse(outputFile);
        final ReadCountCollection expectedCounts = ReadCountCollectionUtils.parse(EXPECTED_COUNTS_MAX_OF_9_FILE);

        Assert.assertEquals(outputCounts.columnNames(), expectedCounts.columnNames());
        Assert.assertEquals(outputCounts.counts().getData(), expectedCounts.counts().getData());
    }

    @Test
    public void testSimpleRunWithMinimumMQ() throws IOException {
        final File outputFile = createTempFile("ctc-test-",".tsv");
        runCommandLine(
                new String[]{
                        "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME,
                        TARGETS_FILE.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12878.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12778.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12872.toString(),
                        "-" + CalculateTargetCoverage.MINIMUM_MAPPING_QUALITY_SHORT_NAME,
                        "30",
                        "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
                        outputFile.toString()
                });
        Assert.assertTrue(outputFile.canRead());
        final ReadCountCollection outputCounts = ReadCountCollectionUtils.parse(outputFile);
        final ReadCountCollection expectedCounts = ReadCountCollectionUtils.parse(EXPECTED_COUNTS_MIN_MQ_30_FILE);

        Assert.assertEquals(outputCounts.columnNames(), expectedCounts.columnNames());
        Assert.assertEquals(outputCounts.counts().getData(), expectedCounts.counts().getData(), Arrays.toString(outputCounts.counts().getData()[0]));
    }

    @Test(expectedExceptions = UserException.class)
    public void testWithoutTargetFile() throws IOException {
        final File outputFile = createTempFile("ctc-test-", ".tsv");
        runCommandLine(
                new String[]{
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12878.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12778.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12872.toString(),
                        "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
                        outputFile.toString()
                });
   }

    @Test(expectedExceptions = UserException.class)
    public void testWithMissingTargetFile() throws IOException {
        final File outputFile = createTempFile("ctc-test-", ".tsv");
        runCommandLine(
                new String[]{
                        "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME,
                        INEXISTENT_TARGETS_FILE.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12878.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12778.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12872.toString(),
                        "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
                        outputFile.toString()
                });
    }

    @Test(expectedExceptions = UserException.class)
    public void testWithTargetFileWithoutCoordinates() throws IOException {
        final File outputFile = createTempFile("ctc-test-", ".tsv");
        runCommandLine(
                new String[]{
                        "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME,
                        TARGETS_FILE_WITHOUT_COORDINATES.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12878.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12778.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12872.toString(),
                        "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
                        outputFile.toString()
                });
    }
}
