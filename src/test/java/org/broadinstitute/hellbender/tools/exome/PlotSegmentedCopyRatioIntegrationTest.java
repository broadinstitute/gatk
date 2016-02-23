package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class PlotSegmentedCopyRatioIntegrationTest extends CommandLineProgramTest{
    private static File testDir = new File("src/test/resources/org/broadinstitute/hellbender/utils/plotter");
    private static final File TN_FILE = new File(testDir, "HCC1143.tn.tsv");
    private static final File PRE_TN_FILE = new File(testDir, "HCC1143.preTN.tsv");
    private static final File SEG_FILE = new File(testDir, "HCC1143.seg");

    private static final String SAMPLE_NAME = "HCC1143";

    @Test()
    public void testUnLoggedCommandLine() {
        File tDir = IOUtils.tempDir("Test", "Plotting");
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, TN_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, PRE_TN_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, SEG_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, tDir.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.LOG2_SHORT_NAME,
        };
        runCommandLine(arguments);
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_FullGenome.png").exists());
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_FullGenome.png").length() > 0);
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_Before_After_CR_Lim_4.png").exists());
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_Before_After_CR_Lim_4.png").length() > 0);
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_Before_After.png").exists());
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_Before_After.png").length() > 0);
    }

    @Test()
    public void testUnLoggedCommandLineSexChrs() {
        File tDir = IOUtils.tempDir("Test", "Plotting");
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, TN_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, PRE_TN_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, SEG_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, tDir.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.LOG2_SHORT_NAME,
                "-" + ExomeStandardArgumentDefinitions.INCLUDE_SEX_CHROMOSOMES_SHORT_NAME,
        };
        runCommandLine(arguments);
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_FullGenome.png").exists());
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_FullGenome.png").length() > 0);
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_Before_After_CR_Lim_4.png").exists());
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_Before_After_CR_Lim_4.png").length() > 0);
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_Before_After.png").exists());
        Assert.assertTrue(new File(tDir, SAMPLE_NAME + "_Before_After.png").length() > 0);
    }
}
