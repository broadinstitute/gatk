package org.broadinstitute.hellbender.utils.plotter;

import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class AlleleFractionPlotterTest extends BaseTest {
    @Test
    public void testAlleleFractionPlotting() throws IOException {
        final File SNP_COUNTS_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/plotter/snps.tsv");
        final File SEGMENTS_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/plotter/allelic.seg");
        String sampleName = "HCC1143";
        File tDir = IOUtils.tempDir("Test", "Plotting");
        AlleleFractionSegmentedPlotter.writeSegmentedAlleleFractionPlot(sampleName, SNP_COUNTS_FILE.getAbsolutePath(),
                SEGMENTS_FILE.getAbsolutePath(), tDir.getAbsolutePath()+"/", false);
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_FullGenome.png").exists());
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_FullGenome.png").length() > 0);
    }

    @Test
    public void testAlleleFractionPlottingSexChrs() throws IOException {
        final File SNP_COUNTS_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/plotter/snps.tsv");
        final File SEGMENTS_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/plotter/allelic.seg");
        String sampleName = "HCC1143";
        File tDir = IOUtils.tempDir("Test", "Plotting");
        AlleleFractionSegmentedPlotter.writeSegmentedAlleleFractionPlot(sampleName, SNP_COUNTS_FILE.getAbsolutePath(),
                SEGMENTS_FILE.getAbsolutePath(), tDir.getAbsolutePath()+"/", true);
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_FullGenome.png").exists());
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_FullGenome.png").length() > 0);
    }
}
