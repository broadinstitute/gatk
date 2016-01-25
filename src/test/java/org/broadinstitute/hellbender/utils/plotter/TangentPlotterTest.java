package org.broadinstitute.hellbender.utils.plotter;

import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class TangentPlotterTest extends BaseTest {

    @Test
    public void testTangentPlotting() throws IOException {
        final File INPUT_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/input/HCC1143_reduced.tsv");
        final File SEG_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/output/HCC1143_reduced_result.seg");
        String sampleName = "HCC1143";
        File tDir = IOUtils.tempDir("Test", "Plotting");
        System.out.println(tDir);
        CopyRatioSegmentedPlotter.writeSegmentedCopyRatioPlot(sampleName, INPUT_FILE.getAbsolutePath(),
                INPUT_FILE.getAbsolutePath(), SEG_FILE.getAbsolutePath(), tDir.getAbsolutePath()+"/", false, false);
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

    @Test
    public void testTangentPlottingSexChrs() throws IOException {
        final File INPUT_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/input/HCC1143_reduced.tsv");
        final File SEG_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/output/HCC1143_reduced_result.seg");
        String sampleName = "HCC1143";
        File tDir = IOUtils.tempDir("Test", "Plotting");
        System.out.println(tDir);
        CopyRatioSegmentedPlotter.writeSegmentedCopyRatioPlot(sampleName, INPUT_FILE.getAbsolutePath(),
                INPUT_FILE.getAbsolutePath(), SEG_FILE.getAbsolutePath(), tDir.getAbsolutePath()+"/", false, true);
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
