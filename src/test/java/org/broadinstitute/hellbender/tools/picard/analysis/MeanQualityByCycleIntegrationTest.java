package org.broadinstitute.hellbender.tools.picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public final class MeanQualityByCycleIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/analysis/MeanQualityByCycle");

    @DataProvider(name="filenames")
    public Object[][] filenames() {
        return new String[][]{
                {"first5000a.bam", null},
                {"first5000a.cram", b37_reference_20_21}
        };
    }

    @Test(dataProvider="filenames", groups = {"R"})
    public void test(final String inputFile, final String referenceName) throws IOException {
        final File input = new File(TEST_DATA_DIR, inputFile);
        final File expectedFile = new File(TEST_DATA_DIR, "meanqualbycycle.txt");
        final File outfile = BaseTest.createTempFile("testMeanQualityByCycle", ".metrics");
        final File pdf = BaseTest.createTempFile("testMeanQualityByCycle", ".pdf");
        outfile.deleteOnExit();
        pdf.deleteOnExit();
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--input");
        args.add(input.getAbsolutePath());
        args.add("--output");
        args.add(outfile.getAbsolutePath());
        args.add("--CHART");
        args.add(pdf.getAbsolutePath());
        args.add("--PRODUCE_PLOT");
        args.add("true");
        if (null != referenceName) {
            final File REF = new File(referenceName);
            args.add("-R");
            args.add(REF.getAbsolutePath());
        }

        runCommandLine(args.getArgsArray());

        try (final FileReader actualReader = new FileReader(outfile);) {
            final MetricsFile<?,Integer> output = new MetricsFile<>();
            output.read(actualReader);
            Assert.assertEquals(output.getAllHistograms().size(), 1);
            Assert.assertEquals(output.getHistogram().size(), 202);
        }
        Assert.assertTrue(pdf.exists(), "exists");
        Assert.assertTrue(pdf.length() > 0, "length");
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
                "--input", input.getAbsolutePath(),
                "--output", outfile.getAbsolutePath(),
                "--CHART", pdf.getAbsolutePath(),
                "--PF_READS_ONLY", "false",
                "--PRODUCE_PLOT", "true",
        };
        runCommandLine(args);
        Assert.assertTrue(pdf.exists(), "exists");
        Assert.assertTrue(pdf.length() > 0, "length");
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
                "--input", input.getAbsolutePath(),
                "--output", outfile.getAbsolutePath(),
                "--CHART", pdf.getAbsolutePath(),
                "--PF_READS_ONLY", "true",
                "--PRODUCE_PLOT", "true",
        };
        runCommandLine(args);
        Assert.assertTrue(pdf.exists(), "exists");
        Assert.assertTrue(pdf.length() > 0, "length");
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
                "--input", input.getAbsolutePath(),
                "--output", outfile.getAbsolutePath(),
                "--CHART", pdf.getAbsolutePath(),
                "--ALIGNED_READS_ONLY", "false",
                "--PRODUCE_PLOT", "true",
        };
        runCommandLine(args);
        Assert.assertTrue(pdf.exists(), "exists");
        Assert.assertTrue(pdf.length() > 0, "length");
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
                "--input", input.getAbsolutePath(),
                "--output", outfile.getAbsolutePath(),
                "--CHART", pdf.getAbsolutePath(),
                "--ALIGNED_READS_ONLY", "true",
                "--PRODUCE_PLOT", "true",
        };
        runCommandLine(args);
        Assert.assertTrue(pdf.exists(), "exists");
        Assert.assertEquals(pdf.length(), 0, "length");   //should be empty because there are no aligned reads here
        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
    }
}
