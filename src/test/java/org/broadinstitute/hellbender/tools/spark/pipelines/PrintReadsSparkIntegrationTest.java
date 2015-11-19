package org.broadinstitute.hellbender.tools.spark.pipelines;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

public final class PrintReadsSparkIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return PrintReadsSpark.class.getSimpleName();
    }

    @Test
    public void testCoordinateSorted() throws Exception {
        final File inBam = new File(getTestDataDir(), "print_reads.sorted.bam");
        final File outBam = createTempFile("print_reads_spark", ".bam");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(inBam.getCanonicalPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outBam.getCanonicalPath());

        this.runCommandLine(args.getArgsArray());

        SamAssertionUtils.assertSamsEqual(outBam, inBam);
    }

    @Test
    public void testCoordinateSortedInRegion() throws Exception {
        final File inBam = new File(getTestDataDir(), "print_reads.sorted.bam");
        final File expectedBam = new File(getTestDataDir(), "print_reads.sorted.chr1_1.bam");
        final File outBam = createTempFile("print_reads_spark", ".bam");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(inBam.getCanonicalPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outBam.getCanonicalPath());
        args.add("-L chr7:1-100 -XL chr7:2-100");

        this.runCommandLine(args.getArgsArray());

        SamAssertionUtils.assertSamsEqual(outBam, expectedBam);
    }

    @Test(expectedExceptions = UserException.class)
    public void testNameSorted() throws Exception {
        final File inBam = new File(getTestDataDir(), "print_reads.bam");
        final File outBam = createTempFile("print_reads", ".bam");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(inBam.getCanonicalPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outBam.getCanonicalPath());

        this.runCommandLine(args.getArgsArray());

        SamAssertionUtils.assertSamsEqual(outBam, inBam);
    }

    /**
     * Test that PrintReadsSpark is correctly applying the WellformedReadFilter
     */
    @Test
    public void testReadFiltering() throws IOException {
        final File samWithOneMalformedRead = new File(getTestDataDir(), "print_reads_one_malformed_read.sam");
        final File outBam = createTempFile("print_reads_testReadFiltering", ".bam");

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(samWithOneMalformedRead.getCanonicalPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outBam.getCanonicalPath());

        runCommandLine(args.getArgsArray());
        SamAssertionUtils.assertSamsEqual(outBam, new File(getTestDataDir(), "expected.print_reads_one_malformed_read.bam"));
    }
}