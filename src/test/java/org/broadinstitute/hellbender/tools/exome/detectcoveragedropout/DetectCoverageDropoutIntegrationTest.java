package org.broadinstitute.hellbender.tools.exome.detectcoveragedropout;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * @author lichtens &lt;lichtens@broadinstitute.org&gt;
 */
public class DetectCoverageDropoutIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome/detectcoveragedropout");
    private static final File TEST_TARGETS = new File(TEST_DIR, "test.tn.HCC1143T-100_27M_37M.tsv");
    private static final File TEST_SEGMENTS = new File(TEST_DIR, "HCC1143T-100_27M_37M.seg");

    @Test
    public void testDetectCoverageDropout() throws IOException {
        final File outputFile = createTempFile("test", ".txt");

        final String[] arguments = {
                "-" + DetectCoverageDropout.SEGFILE_SHORT_NAME, TEST_SEGMENTS.getAbsolutePath(),
                "-" + DetectCoverageDropout.TARGET_FILE_SHORT_NAME, TEST_TARGETS.getAbsolutePath(),
                "-" + DetectCoverageDropout.OUTPUT_FILE_SHORT_NAME, outputFile.getAbsolutePath(),
        };
        runCommandLine(arguments);

        List<CoverageDropoutResult> results = CoverageDropoutResult.readCoverageDropoutResultsFromTsv(outputFile);

        Assert.assertTrue(results.size() == 1, "Should have only had one entry for this file.");

        CoverageDropoutResult result = results.get(0);
        Assert.assertTrue(result.getNumSegments() == 10, "This testdata had ten segments, but " + result.getNumSegments() + " seen.");
        Assert.assertTrue(result.getNumSegmentsAfterFiltering() == 4, "This testdata had four segments with min number of targets, but " + result.getNumSegmentsAfterFiltering() + " seen.");
        Assert.assertFalse(result.isCoverageDropout(), "This sample did not have coverage dropout, yet was tagged as such.");
        Assert.assertTrue(result.getNumSegmentsAfterFiltering() < result.getNumSegments(), "This sample should have filtered segments, yet number of filtered and unfiltered segments are equal.");
        Assert.assertFalse(new File(outputFile.getAbsoluteFile().getParent()+"/"+DetectCoverageDropout.WARNING_FILE+outputFile.getName()).exists(), "Warning file exists, but it should not have been written.");
    }

    @Test
    public void testDetectCoverageDropoutFails() throws IOException {
        final File outputFile = createTempFile("testFails", ".txt");

        final String[] arguments = {
                "-" + DetectCoverageDropout.SEGFILE_SHORT_NAME, TEST_SEGMENTS.getAbsolutePath(),
                "-" + DetectCoverageDropout.TARGET_FILE_SHORT_NAME, TEST_TARGETS.getAbsolutePath(),
                "-" + DetectCoverageDropout.OUTPUT_FILE_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + DetectCoverageDropout.SEG_THRESHOLD_SHORT_NAME, Double.toString(1.0),
        };
        runCommandLine(arguments);

        List<CoverageDropoutResult> results = CoverageDropoutResult.readCoverageDropoutResultsFromTsv(outputFile);

        Assert.assertTrue(results.size() == 1, "Should have only had one entry for this file.");

        CoverageDropoutResult result = results.get(0);
        Assert.assertTrue(result.getNumSegments() == 10, "This testdata had ten segments, but " + result.getNumSegments() + " seen.");
        Assert.assertTrue(result.getNumSegmentsAfterFiltering() == 4, "This testdata had four segments with min number of targets, but " + result.getNumSegmentsAfterFiltering() + " seen.");
        Assert.assertTrue(result.isCoverageDropout(), "This sample did not have coverage dropout, yet was tagged as such.");
        Assert.assertTrue(result.getNumSegmentsAfterFiltering() < result.getNumSegments(), "This sample should have filtered segments, yet number of filtered and unfiltered segments are equal.");
        Assert.assertTrue(new File(outputFile.getAbsoluteFile().getParent() + "/" + DetectCoverageDropout.WARNING_FILE + outputFile.getName()).exists(), "Warning file does not exist, but it should have been written.");
    }
}