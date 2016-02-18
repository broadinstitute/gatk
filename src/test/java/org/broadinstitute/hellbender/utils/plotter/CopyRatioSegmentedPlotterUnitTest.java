package org.broadinstitute.hellbender.utils.plotter;

import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class CopyRatioSegmentedPlotterUnitTest extends BaseTest {
    private static File testDir = new File("src/test/resources/org/broadinstitute/hellbender/utils/plotter");
    private static final File TN_FILE = new File(testDir, "HCC1143.tn.tsv");
    private static final File SEG_FILE = new File(testDir, "HCC1143.seg");
    private static final String SAMPLE_NAME = "HCC1143";

    @Test
    public void testTangentPlotting() throws IOException {
        File tDir = IOUtils.tempDir("Test", "Plotting");
        CopyRatioSegmentedPlotter.writeSegmentedCopyRatioPlot(SAMPLE_NAME, TN_FILE.getAbsolutePath(),
                TN_FILE.getAbsolutePath(), SEG_FILE.getAbsolutePath(), tDir.getAbsolutePath()+"/", false, false);
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_FullGenome.png").exists());
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_FullGenome.png").length() > 0);
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_Before_After_CR_Lim_4.png").exists());
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_Before_After_CR_Lim_4.png").length() > 0);
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_Before_After.png").exists());
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_Before_After.png").length() > 0);
    }

    @Test
    public void testTangentPlottingSexChrs() throws IOException {
        File tDir = IOUtils.tempDir("Test", "Plotting");
        CopyRatioSegmentedPlotter.writeSegmentedCopyRatioPlot(SAMPLE_NAME, TN_FILE.getAbsolutePath(),
                TN_FILE.getAbsolutePath(), SEG_FILE.getAbsolutePath(), tDir.getAbsolutePath()+"/", false, true);
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_FullGenome.png").exists());
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_FullGenome.png").length() > 0);
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_Before_After_CR_Lim_4.png").exists());
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_Before_After_CR_Lim_4.png").length() > 0);
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_Before_After.png").exists());
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_Before_After.png").length() > 0);
    }
}
