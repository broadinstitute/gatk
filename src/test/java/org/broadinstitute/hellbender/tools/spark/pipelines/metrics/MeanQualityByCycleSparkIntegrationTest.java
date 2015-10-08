package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;

public final class MeanQualityByCycleSparkIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/analysis/MeanQualityByCycle");

    @Override
    public String getTestedClassName() {
        return MeanQualityByCycleSpark.class.getSimpleName();
    }

    @Test
    public void test1() throws Exception {
        //Note we compare to non-spark outputs
        final File unsortedBam = new File(TEST_DATA_DIR, "first5000a.bam");
        final File expectedFile = new File(TEST_DATA_DIR, "meanqualbycycle.txt");
        final File outfile = BaseTest.createTempFile("testMeanQualityByCycle", ".metrics");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(unsortedBam.getCanonicalPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outfile.getCanonicalPath());
        this.runCommandLine(args.getArgsArray());

        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
    }

    @Test
    public void test_PF_READS_ONLY_false() throws Exception {
        //Note we compare to non-spark outputs
        final File unsortedBam = new File(TEST_DATA_DIR, "example_pfFail_reads.bam");
        final File expectedFile = new File(TEST_DATA_DIR, "pfFailBam.pf.txt");
        final File outfile = BaseTest.createTempFile("pfFailBam.pf.", ".metrics");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(unsortedBam.getCanonicalPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outfile.getCanonicalPath());
        args.add("--" + "PF_READS_ONLY");
        args.add("false");
        this.runCommandLine(args.getArgsArray());

        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
    }

    @Test
    public void test_PF_READS_ONLY_true() throws Exception {
        //Note we compare to non-spark outputs
        final File unsortedBam = new File(TEST_DATA_DIR, "example_pfFail_reads.bam");
        final File expectedFile = new File(TEST_DATA_DIR, "pfFailBam.pfOnly.txt");
        final File outfile = BaseTest.createTempFile("pfFailBam.pf.", ".metrics");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(unsortedBam.getCanonicalPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outfile.getCanonicalPath());
        args.add("--" + "PF_READS_ONLY");
        args.add("true");
        this.runCommandLine(args.getArgsArray());

        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
    }

    @Test
    public void test_ALIGNED_READS_ONLY_false() throws Exception {
        //Note we compare to non-spark outputs
        final File unsortedBam = new File(TEST_DATA_DIR, "unmapped.bam");
        final File expectedFile = new File(TEST_DATA_DIR, "unmappedBam.ALIGNED_READS_ONLY_false.txt");
        final File outfile = BaseTest.createTempFile("unmappedBam.ALIGNED_READS_ONLY_false.", ".metrics");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(unsortedBam.getCanonicalPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outfile.getCanonicalPath());
        args.add("--" + "ALIGNED_READS_ONLY");
        args.add("false");
        this.runCommandLine(args.getArgsArray());

        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
    }

    @Test
    public void test_ALIGNED_READS_ONLY_true() throws Exception {
        //Note we compare to non-spark outputs
        final File unsortedBam = new File(TEST_DATA_DIR, "unmapped.bam");
        final File expectedFile = new File(TEST_DATA_DIR, "unmappedBam.ALIGNED_READS_ONLY_true.txt");
        final File outfile = BaseTest.createTempFile("unmappedBam.ALIGNED_READS_ONLY_true.", ".metrics");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(unsortedBam.getCanonicalPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outfile.getCanonicalPath());
        args.add("--" + "ALIGNED_READS_ONLY");
        args.add("true");
        this.runCommandLine(args.getArgsArray());

        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
    }
}
