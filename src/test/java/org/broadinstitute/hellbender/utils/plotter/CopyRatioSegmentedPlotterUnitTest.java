package org.broadinstitute.hellbender.utils.plotter;

import org.broadinstitute.hellbender.utils.R.RScriptExecutorException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
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
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_preQc.txt").exists());
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_preQc.txt").length() > 0);
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_postQc.txt").exists());
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_postQc.txt").length() > 0);
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_dQc.txt").exists());
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_dQc.txt").length() > 0);
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_scaled_dQc.txt").exists());
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_scaled_dQc.txt").length() > 0);
        final double postQc = ParamUtils.readValuesFromFile(new File(tDir, SAMPLE_NAME + "_preQc.txt"))[0];
        final double preQc = ParamUtils.readValuesFromFile(new File(tDir, SAMPLE_NAME + "_postQc.txt"))[0];
        final double dQc = ParamUtils.readValuesFromFile(new File(tDir, SAMPLE_NAME + "_dQc.txt"))[0];
        final double scaled_dQc = ParamUtils.readValuesFromFile(new File(tDir, SAMPLE_NAME + "_scaled_dQc.txt"))[0];
        Assert.assertEquals(dQc, preQc - postQc);
        Assert.assertEquals(scaled_dQc, (preQc - postQc)/preQc);
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
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_preQc.txt").exists());
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_preQc.txt").length() > 0);
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_postQc.txt").exists());
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_postQc.txt").length() > 0);
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_dQc.txt").exists());
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_dQc.txt").length() > 0);
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_scaled_dQc.txt").exists());
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_scaled_dQc.txt").length() > 0);
        final double postQc = ParamUtils.readValuesFromFile(new File(tDir, SAMPLE_NAME + "_preQc.txt"))[0];
        final double preQc = ParamUtils.readValuesFromFile(new File(tDir, SAMPLE_NAME + "_postQc.txt"))[0];
        final double dQc = ParamUtils.readValuesFromFile(new File(tDir, SAMPLE_NAME + "_dQc.txt"))[0];
        final double scaled_dQc = ParamUtils.readValuesFromFile(new File(tDir, SAMPLE_NAME + "_scaled_dQc.txt"))[0];
        Assert.assertEquals(dQc, preQc - postQc);
        Assert.assertEquals(scaled_dQc, (preQc - postQc)/preQc);
    }

    @Test(expectedExceptions = RScriptExecutorException.class)
    public void testMissingInput() throws IOException {
        File tDir = IOUtils.tempDir("Test", "Plotting");
        CopyRatioSegmentedPlotter.writeSegmentedCopyRatioPlot(SAMPLE_NAME, "Non-existent-file",
                TN_FILE.getAbsolutePath(), SEG_FILE.getAbsolutePath(), tDir.getAbsolutePath()+"/", false, false);
    }
}
