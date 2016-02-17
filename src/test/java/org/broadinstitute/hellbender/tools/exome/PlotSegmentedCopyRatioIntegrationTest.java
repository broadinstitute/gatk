package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class PlotSegmentedCopyRatioIntegrationTest extends CommandLineProgramTest{

    @Test()
    public void testUnLoggedCommandLine() {
        final File TN_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/plotter/HCC1143.tn.tsv");
        final File PRE_TN_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/plotter/HCC1143.preTN.tsv");
        final File SEG_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/plotter/HCC1143.seg");
        final String sampleName = "HCC1143";
        File tDir = IOUtils.tempDir("Test", "Plotting");
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME, TN_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, PRE_TN_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, SEG_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, tDir.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.LOG2_SHORT_NAME,
        };
        runCommandLine(arguments);
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_FullGenome.png").exists());
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_FullGenome.png").length() > 0);
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_Before_After_CR_Lim_4.png").exists());
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_Before_After_CR_Lim_4.png").length() > 0);
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_Before_After.png").exists());
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_Before_After.png").length() > 0);
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_preQc.txt").exists());
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_preQc.txt").length() > 0);
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_postQc.txt").exists());
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_postQc.txt").length() > 0);
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_dQc.txt").exists());
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_dQc.txt").length() > 0);
    }
    @Test()
    public void testUnLoggedCommandLineSexChrs() {
        final File TN_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/plotter/HCC1143.tn.tsv");
        final File PRE_TN_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/plotter/HCC1143.preTN.tsv");
        final File SEG_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/plotter/HCC1143.seg");
        final String sampleName = "HCC1143";
        File tDir = IOUtils.tempDir("Test", "Plotting");
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME, TN_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, PRE_TN_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, SEG_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, tDir.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.LOG2_SHORT_NAME,
                "-" + ExomeStandardArgumentDefinitions.INCLUDE_SEX_CHROMOSOMES_SHORT_NAME,
        };
        runCommandLine(arguments);
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_FullGenome.png").exists());
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_FullGenome.png").length() > 0);
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_Before_After_CR_Lim_4.png").exists());
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_Before_After_CR_Lim_4.png").length() > 0);
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_Before_After.png").exists());
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_Before_After.png").length() > 0);
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_preQc.txt").exists());
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_preQc.txt").length() > 0);
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_postQc.txt").exists());
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_postQc.txt").length() > 0);
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_dQc.txt").exists());
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_dQc.txt").length() > 0);
    }
}
