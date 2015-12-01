package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.plotter.CopyRatioSegmentedPlotter;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class PlotSegmentedCopyRatioIntegrationTest extends CommandLineProgramTest{

    @Test()
    public void testUnLoggedCommandLine() throws IOException {
        final File TN_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/plotter/HCC1143.tn.tsv");
        final File PRE_TN_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/plotter/HCC1143.preTN.tsv");
        final File SEG_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/plotter/HCC1143.seg");
        final String sampleName = "HCC1143";
        File tDir = IOUtils.tempDir("Test", "Plotting");
        CopyRatioSegmentedPlotter.writeSegmentedCopyRatioPlot(sampleName, TN_FILE.getAbsolutePath(),
                PRE_TN_FILE.getAbsolutePath(), SEG_FILE.getAbsolutePath(), tDir.getAbsolutePath() + "/", true);
        final String[] arguments = {
                "-" + PlotSegmentedCopyRatio.SAMPLE_NAME_SHORT_NAME, sampleName,
                "-" + PlotSegmentedCopyRatio.TARGETS_FILE_SHORT_NAME, TN_FILE.getAbsolutePath(),
                "-" + PlotSegmentedCopyRatio.PRE_TANGENT_TARGETS_FILE_SHORT_NAME, PRE_TN_FILE.getAbsolutePath(),
                "-" + PlotSegmentedCopyRatio.SEGMENT_FILE_SHORT_NAME, SEG_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, tDir.getAbsolutePath(),
                "-" + PlotSegmentedCopyRatio.LOG2_SHORT_NAME,
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
