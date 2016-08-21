package org.broadinstitute.hellbender.utils.plotter;

import org.broadinstitute.hellbender.utils.R.RScriptExecutorException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class ACNVPlotterUnitTest extends BaseTest {
    private static final String TOOLS_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/tools/exome/";
    private static final File TANGENT_NORMALIZED_COVERAGE_FILE = new File(TOOLS_TEST_DIRECTORY, "coverages-for-allelic-integration.tsv");
    private static final File TUMOR_ALLELIC_COUNTS_FILE = new File(TOOLS_TEST_DIRECTORY, "snps-for-allelic-integration.tsv");
    private static final File SEGMENTS_FILE = new File(TOOLS_TEST_DIRECTORY, "acnv-segments-from-allelic-integration.seg");
    private static final String SAMPLE_NAME = "test"; // This is what is in the segments file
    private static final String OUTPUT_PREFIX = "sample";

    @Test
    public void testACNVPlotting() {
        final File tDir = IOUtils.tempDir("Test", "Plotting");
        ACNVPlotter.writeSegmentedAlleleFractionPlot(SAMPLE_NAME, TUMOR_ALLELIC_COUNTS_FILE,
                TANGENT_NORMALIZED_COVERAGE_FILE, SEGMENTS_FILE, tDir.getAbsolutePath() + "/", OUTPUT_PREFIX, false);
        Assert.assertTrue(new File(tDir, OUTPUT_PREFIX + "_ACNV.png").exists());
        Assert.assertTrue(new File(tDir, OUTPUT_PREFIX + "_ACNV.png").length() > 0);
    }

    @Test
    public void testACNVPlottingPerSeg() {
        final File tDir = IOUtils.tempDir("Test", "Plotting");
        ACNVPlotter.writeSegmentedAlleleFractionPlotPerSeg(SAMPLE_NAME, TUMOR_ALLELIC_COUNTS_FILE,
                TANGENT_NORMALIZED_COVERAGE_FILE, SEGMENTS_FILE, tDir.getAbsolutePath() + "/", OUTPUT_PREFIX, false);
        Assert.assertTrue(new File(tDir, OUTPUT_PREFIX + "_Chr_1_30116_62393753.png").exists());
        Assert.assertTrue(new File(tDir, OUTPUT_PREFIX + "_Chr_1_30116_62393753.png").length() > 0);
        Assert.assertTrue(new File(tDir, OUTPUT_PREFIX + "_Chr_10_92745_49997881.png").exists());
        Assert.assertTrue(new File(tDir, OUTPUT_PREFIX + "_Chr_10_92745_49997881.png").length() > 0);
    }

    @Test
    public void testACNVSexChrs() {
        final File tDir = IOUtils.tempDir("Test", "Plotting");
        ACNVPlotter.writeSegmentedAlleleFractionPlot(SAMPLE_NAME, TUMOR_ALLELIC_COUNTS_FILE,
                TANGENT_NORMALIZED_COVERAGE_FILE, SEGMENTS_FILE, tDir.getAbsolutePath() + "/", OUTPUT_PREFIX, true);
        Assert.assertTrue(new File(tDir, OUTPUT_PREFIX + "_ACNV.png").exists());
        Assert.assertTrue(new File(tDir, OUTPUT_PREFIX + "_ACNV.png").length() > 0);
    }

    @Test(expectedExceptions = RScriptExecutorException.class)
    public void testMissingInput() throws IOException {
        final File tDir = IOUtils.tempDir("Test", "Plotting");
        ACNVPlotter.writeSegmentedAlleleFractionPlot(SAMPLE_NAME, new File("Non-existent-file"),
                TANGENT_NORMALIZED_COVERAGE_FILE, SEGMENTS_FILE, tDir.getAbsolutePath() + "/", OUTPUT_PREFIX, true);
    }
}
