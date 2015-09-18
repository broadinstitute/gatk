package org.broadinstitute.hellbender.tools.spark.pipelines;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.read.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.annotations.Test;

import java.io.File;

public final class SortBamSparkIntegrationTest extends CommandLineProgramTest {

    @Test
    public void test() throws Exception {
        final File unsortedBam = new File(getTestDataDir(), "count_reads.bam");
        final File sortedBam = new File(getTestDataDir(), "count_reads_sorted.bam");
        final File outputBam = createTempFile("sort_bam_spark", ".bam");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--"+ StandardArgumentDefinitions.INPUT_LONG_NAME); args.add(unsortedBam.getCanonicalPath());
        args.add("--"+StandardArgumentDefinitions.OUTPUT_LONG_NAME); args.add(outputBam.getCanonicalPath());
        args.add("--parallelism"); args.add("1");

        this.runCommandLine(args.getArgsArray());

        SamAssertionUtils.assertSamsEqual(outputBam, sortedBam);
    }

}