package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class PlotACNVResultsIntegrationTest extends CommandLineProgramTest {
    private static final String TOOLS_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/tools/exome/";
    private static final File TANGENT_NORMALIZED_COVERAGE_FILE = new File(TOOLS_TEST_DIRECTORY, "coverages-for-allelic-integration.tsv");
    private static final File TUMOR_ALLELIC_COUNTS_FILE = new File(TOOLS_TEST_DIRECTORY, "snps-for-allelic-integration.tsv");
    private static final File SEGMENTS_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/plotter/ACNV_final_segments.seg");
    private static final String OUTPUT_PREFIX = "sample";

    @Test
    public void testACNVPlotting() {
        File tDir = IOUtils.tempDir("Test", "Plotting");
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.ALLELIC_COUNTS_FILE_SHORT_NAME, TUMOR_ALLELIC_COUNTS_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, TANGENT_NORMALIZED_COVERAGE_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, SEGMENTS_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, tDir.getAbsolutePath(),
                "-" + PlotACNVResults.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX,
        };
        runCommandLine(arguments);
        Assert.assertTrue(new File(tDir, OUTPUT_PREFIX + "_ACNV.png").exists());
        Assert.assertTrue(new File(tDir, OUTPUT_PREFIX + "_ACNV.png").length() > 0);
        Assert.assertTrue(new File(tDir, OUTPUT_PREFIX + "_Chr_1_30116_62393753.png").exists());
        Assert.assertTrue(new File(tDir, OUTPUT_PREFIX + "_Chr_1_30116_62393753.png").length() > 0);
        Assert.assertTrue(new File(tDir, OUTPUT_PREFIX + "_Chr_10_92745_49997881.png").exists());
        Assert.assertTrue(new File(tDir, OUTPUT_PREFIX + "_Chr_10_92745_49997881.png").length() > 0);
    }

    @Test
    public void testACNVPlottingSexChrs() {
        File tDir = IOUtils.tempDir("Test", "Plotting");
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.ALLELIC_COUNTS_FILE_SHORT_NAME, TUMOR_ALLELIC_COUNTS_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, TANGENT_NORMALIZED_COVERAGE_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, SEGMENTS_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, tDir.getAbsolutePath(),
                "-" + PlotACNVResults.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX,
                "-" + ExomeStandardArgumentDefinitions.INCLUDE_SEX_CHROMOSOMES_SHORT_NAME,
        };
        runCommandLine(arguments);
        Assert.assertTrue(new File(tDir, OUTPUT_PREFIX + "_ACNV.png").exists());
        Assert.assertTrue(new File(tDir, OUTPUT_PREFIX + "_ACNV.png").length() > 0);
    }
}