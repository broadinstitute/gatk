package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public final class MeanQualityByCycleSparkIntegrationTest extends CommandLineProgramTest {

    //NOTE: these tests use the same data and results as the non-spark ones, by design

    private static final File TEST_DATA_DIR = new File(
            "src/test/resources/org/broadinstitute/hellbender/metrics/analysis//MeanQualityByCycle");

    @Override
    public String getTestedClassName() {
        return MeanQualityByCycleSpark.class.getSimpleName();
    }

    @DataProvider(name="filenames")
    public Object[][] filenames() {
        return new String[][]{
                {"first5000a.bam", null},
                {"first5000a.cram", b37_reference_20_21}
        };
    }
    @Test(dataProvider="filenames", groups = {"R", "spark"})
    public void test(final String inputFile, final String referenceName) throws IOException {
        final File input = new File(TEST_DATA_DIR, inputFile);
        final File expectedFile = new File(TEST_DATA_DIR, "meanqualbycycle.txt");
        final File outfile = GATKBaseTest.createTempFile("testMeanQualityByCycle", ".metrics");
        final File pdf = GATKBaseTest.createTempFile("testMeanQualityByCycle", ".pdf");
        outfile.deleteOnExit();
        pdf.deleteOnExit();
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--input");
        args.add(input.getAbsolutePath());
        args.add("--output");
        args.add(outfile.getAbsolutePath());
        args.add("--chart");
        args.add(pdf.getAbsolutePath());
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

    @Test( groups = "spark")
    public void test1() throws IOException {
        //Note we compare to non-spark outputs
        final File unsortedBam = new File(TEST_DATA_DIR, "first5000a.bam");
        final File expectedFile = new File(TEST_DATA_DIR, "meanqualbycycle.txt");
        final File outfile = GATKBaseTest.createTempFile("testMeanQualityByCycle", ".metrics");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(unsortedBam.getCanonicalPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outfile.getCanonicalPath());
        this.runCommandLine(args.getArgsArray());

        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
    }

    //Disabled due to https://github.com/broadinstitute/gatk/issues/1540
    @Test(enabled=false, groups="spark")
    public void test_ADAM() throws IOException {
        //Note we compare to non-spark outputs
        final File adamFile = new File(TEST_DATA_DIR, "first5000a.adam");
        final File expectedFile = new File(TEST_DATA_DIR, "meanqualbycycle.txt");
        final File outfile = GATKBaseTest.createTempFile("testMeanQualityByCycleADAM", ".metrics");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(adamFile.getCanonicalPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outfile.getCanonicalPath());
        this.runCommandLine(args.getArgsArray());

        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
    }

    @Test(groups = "spark")
    public void test_PF_READS_ONLY_false() throws IOException {
        //Note we compare to non-spark outputs
        final File unsortedBam = new File(TEST_DATA_DIR, "example_pfFail_reads.bam");
        final File expectedFile = new File(TEST_DATA_DIR, "pfFailBam.pf.txt");
        final File outfile = GATKBaseTest.createTempFile("pfFailBam.pf.", ".metrics");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(unsortedBam.getCanonicalPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outfile.getCanonicalPath());
        args.add("--" + "pfReadsOnly");
        args.add("false");
        this.runCommandLine(args.getArgsArray());

        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
    }

    @Test(groups = "spark")
    public void test_PF_READS_ONLY_true() throws IOException {
        //Note we compare to non-spark outputs
        final File unsortedBam = new File(TEST_DATA_DIR, "example_pfFail_reads.bam");
        final File expectedFile = new File(TEST_DATA_DIR, "pfFailBam.pfOnly.txt");
        final File outfile = GATKBaseTest.createTempFile("pfFailBam.pf.", ".metrics");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(unsortedBam.getCanonicalPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outfile.getCanonicalPath());
        args.add("--" + "pfReadsOnly");
        args.add("true");
        this.runCommandLine(args.getArgsArray());

        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
    }

    @Test(groups = "spark")
    public void test_ALIGNED_READS_ONLY_false() throws IOException {
        //Note we compare to non-spark outputs
        final File unsortedBam = new File(TEST_DATA_DIR, "unmapped.bam");
        final File expectedFile = new File(TEST_DATA_DIR, "unmappedBam.ALIGNED_READS_ONLY_false.txt");
        final File outfile = GATKBaseTest.createTempFile("unmappedBam.ALIGNED_READS_ONLY_false.", ".metrics");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(unsortedBam.getCanonicalPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outfile.getCanonicalPath());
        args.add("--" + "alignedReadsOnly");
        args.add("false");
        this.runCommandLine(args.getArgsArray());

        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
    }

    @Test(groups = "spark")
    public void test_ALIGNED_READS_ONLY_true() throws IOException {
        //Note we compare to non-spark outputs
        final File unsortedBam = new File(TEST_DATA_DIR, "unmapped.bam");
        final File expectedFile = new File(TEST_DATA_DIR, "unmappedBam.ALIGNED_READS_ONLY_true.txt");
        final File outfile = GATKBaseTest.createTempFile("unmappedBam.ALIGNED_READS_ONLY_true.", ".metrics");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(unsortedBam.getCanonicalPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outfile.getCanonicalPath());
        args.add("--" + "alignedReadsOnly");
        args.add("true");
        this.runCommandLine(args.getArgsArray());

        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
    }

    @Test
    public void testGetRScriptResource() {
        // Make sure the RScript resource can be resolved
        Assert.assertNotNull(MeanQualityByCycleSpark.getMeanQualityByCycleRScriptResource());
    }

}
