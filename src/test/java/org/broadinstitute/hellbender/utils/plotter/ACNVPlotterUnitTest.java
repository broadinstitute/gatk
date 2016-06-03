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
    private static final File SEGMENTS_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/plotter/ACNV_final_segments.seg");
    private static final String SAMPLE_NAME = "test"; // This is what is in the segments file
    @Test
    public void testACNVPlotting() {
        File tDir = IOUtils.tempDir("Test", "Plotting");
        ACNVPlotter.writeSegmentedAlleleFractionPlot(SAMPLE_NAME, TUMOR_ALLELIC_COUNTS_FILE.getAbsolutePath(),
                TANGENT_NORMALIZED_COVERAGE_FILE.getAbsolutePath(),
                SEGMENTS_FILE.getAbsolutePath(), tDir.getAbsolutePath()+"/", false);
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_ACNV.png").exists());
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_ACNV.png").length() > 0);
    }

    @Test
    public void testACNVPlottingPerSeg() {
        File tDir = IOUtils.tempDir("Test", "Plotting");
        ACNVPlotter.writeSegmentedAlleleFractionPlotPerSeg(SAMPLE_NAME, TUMOR_ALLELIC_COUNTS_FILE.getAbsolutePath(),
                TANGENT_NORMALIZED_COVERAGE_FILE.getAbsolutePath(),
                SEGMENTS_FILE.getAbsolutePath(), tDir.getAbsolutePath() + "/", false);
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_Chr_1_30116_62393753.png").exists());
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_Chr_1_30116_62393753.png").length() > 0);
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_Chr_10_92745_49997881.png").exists());
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_Chr_10_92745_49997881.png").length() > 0);
    }

    @Test
    public void testACNVSexChrs() {
        File tDir = IOUtils.tempDir("Test", "Plotting");
        ACNVPlotter.writeSegmentedAlleleFractionPlot(SAMPLE_NAME, TUMOR_ALLELIC_COUNTS_FILE.getAbsolutePath(),
                TANGENT_NORMALIZED_COVERAGE_FILE.getAbsolutePath(),
                SEGMENTS_FILE.getAbsolutePath(), tDir.getAbsolutePath()+"/", true);
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_ACNV.png").exists());
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_ACNV.png").length() > 0);
    }

    @Test(expectedExceptions = RScriptExecutorException.class)
    public void testMissingInput() throws IOException {
        File tDir = IOUtils.tempDir("Test", "Plotting");
        ACNVPlotter.writeSegmentedAlleleFractionPlot(SAMPLE_NAME, "Non-existent-file",
                TANGENT_NORMALIZED_COVERAGE_FILE.getAbsolutePath(),
                SEGMENTS_FILE.getAbsolutePath(), tDir.getAbsolutePath() + "/", true);
    }
}
