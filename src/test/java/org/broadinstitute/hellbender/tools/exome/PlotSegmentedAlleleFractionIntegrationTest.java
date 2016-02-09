package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class PlotSegmentedAlleleFractionIntegrationTest extends CommandLineProgramTest {

    private static final File SNP_COUNTS_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/plotter/snps.tsv");
    private static final File SEGMENTS_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/plotter/allelic.seg");
    private static final String sampleName = "test"; // This is what is in allelic.seg test file

    @Test
    public void testAlleleFractionPlotting() {
        File tDir = IOUtils.tempDir("Test", "Plotting");
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.SNP_FILE_SHORT_NAME, SNP_COUNTS_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, SEGMENTS_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, tDir.getAbsolutePath()
        };
        runCommandLine(arguments);
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_FullGenome.png").exists());
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_FullGenome.png").length() > 0);
    }

    @Test
    public void testAlleleFractionPlottingSexChrs() {
        File tDir = IOUtils.tempDir("Test", "Plotting");
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.SNP_FILE_SHORT_NAME, SNP_COUNTS_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, SEGMENTS_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, tDir.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.INCLUDE_SEX_CHROMOSOMES_SHORT_NAME,
        };
        runCommandLine(arguments);
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_FullGenome.png").exists());
        Assert.assertTrue(new File(tDir + "/" + sampleName + "_FullGenome.png").length() > 0);
    }
}