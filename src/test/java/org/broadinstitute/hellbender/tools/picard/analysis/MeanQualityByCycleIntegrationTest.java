package org.broadinstitute.hellbender.tools.picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public final class MeanQualityByCycleIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/analysis/MeanQualityByCycle");

    @Test(groups = {"R"})
    public void test() throws IOException {
        final File input = new File(TEST_DATA_DIR, "first5000a.bam");
        final File expectedFile = new File(TEST_DATA_DIR, "meanqualbycycle.txt");
        final File outfile = BaseTest.createTempFile("testMeanQualityByCycle", ".metrics");
        final File pdf = BaseTest.createTempFile("testMeanQualityByCycle", ".pdf");
        outfile.deleteOnExit();
        pdf.deleteOnExit();
        final String[] args = {
                "--INPUT", input.getAbsolutePath(),
                "--OUTPUT", outfile.getAbsolutePath(),
                "--CHART", pdf.getAbsolutePath(),
                "--PRODUCE_PLOT", "true",
        };
        runCommandLine(args);

        try (final FileReader actualReader = new FileReader(outfile);) {
            final MetricsFile<?,Integer> output = new MetricsFile<>();
            output.read(actualReader);
            Assert.assertEquals(output.getAllHistograms().size(), 1);
            Assert.assertEquals(output.getHistogram().size(), 202);
        }
        Assert.assertTrue(pdf.exists());
        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
    }

    @Test(groups = {"R"})
    public void test_PF_READS_ONLY_false() throws IOException {
        final File input = new File(TEST_DATA_DIR, "example_pfFail_reads.bam");
        final File expectedFile = new File(TEST_DATA_DIR, "pfFailBam.pf.txt");
        final File outfile = BaseTest.createTempFile("pfFailBam.pf", ".metrics");
        final File pdf = BaseTest.createTempFile("pfFailBam.pf", ".pdf");
        outfile.deleteOnExit();
        pdf.deleteOnExit();
        final String[] args = {
                "--INPUT", input.getAbsolutePath(),
                "--OUTPUT", outfile.getAbsolutePath(),
                "--CHART", pdf.getAbsolutePath(),
                "--PF_READS_ONLY", "false",
                "--PRODUCE_PLOT", "true",
        };
        runCommandLine(args);
        Assert.assertTrue(pdf.exists());
        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
    }

    @Test(groups = {"R"})
    public void test_PF_READS_ONLY_true() throws IOException {
        final File input = new File(TEST_DATA_DIR, "example_pfFail_reads.bam");
        final File expectedFile = new File(TEST_DATA_DIR, "pfFailBam.pfOnly.txt");
        final File outfile = BaseTest.createTempFile("pfFailBam.pf", ".metrics");
        final File pdf = BaseTest.createTempFile("pfFailBam.pf", ".pdf");
        outfile.deleteOnExit();
        pdf.deleteOnExit();
        final String[] args = {
                "--INPUT", input.getAbsolutePath(),
                "--OUTPUT", outfile.getAbsolutePath(),
                "--CHART", pdf.getAbsolutePath(),
                "--PF_READS_ONLY", "true",
                "--PRODUCE_PLOT", "true",
        };
        runCommandLine(args);
        Assert.assertTrue(pdf.exists());
        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
    }

    @Test(groups = {"R"})
    public void test_ALIGNED_READS_ONLY_false() throws IOException {
        final File input = new File(TEST_DATA_DIR, "unmapped.bam");
        final File expectedFile = new File(TEST_DATA_DIR, "unmappedBam.ALIGNED_READS_ONLY_false.txt");
        final File outfile = BaseTest.createTempFile("unmappedBam.ALIGNED_READS_ONLY_false", ".metrics");
        final File pdf = BaseTest.createTempFile("unmappedBam.ALIGNED_READS_ONLY_false", ".pdf");
        outfile.deleteOnExit();
        pdf.deleteOnExit();
        final String[] args = {
                "--INPUT", input.getAbsolutePath(),
                "--OUTPUT", outfile.getAbsolutePath(),
                "--CHART", pdf.getAbsolutePath(),
                "--ALIGNED_READS_ONLY", "false",
                "--PRODUCE_PLOT", "true",
        };
        runCommandLine(args);
        Assert.assertTrue(pdf.exists());
        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
    }

    @Test(groups = {"R"})
    public void test_ALIGNED_READS_ONLY_true() throws IOException {
        final File input = new File(TEST_DATA_DIR, "unmapped.bam");
        final File expectedFile = new File(TEST_DATA_DIR, "unmappedBam.ALIGNED_READS_ONLY_true.txt");
        final File outfile = BaseTest.createTempFile("unmappedBam.ALIGNED_READS_ONLY_true", ".metrics");
        final File pdf = BaseTest.createTempFile("unmappedBam.ALIGNED_READS_ONLY_true", ".pdf");
        outfile.deleteOnExit();
        pdf.deleteOnExit();
        final String[] args = {
                "--INPUT", input.getAbsolutePath(),
                "--OUTPUT", outfile.getAbsolutePath(),
                "--CHART", pdf.getAbsolutePath(),
                "--ALIGNED_READS_ONLY", "true",
                "--PRODUCE_PLOT", "true",
        };
        runCommandLine(args);
        Assert.assertTrue(pdf.exists());
        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
    }
}
