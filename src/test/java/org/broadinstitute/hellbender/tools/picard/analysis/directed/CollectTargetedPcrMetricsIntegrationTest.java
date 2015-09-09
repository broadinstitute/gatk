package org.broadinstitute.hellbender.tools.picard.analysis.directed;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

/**
 * Created by dkling on 9/3/15.
 */

public final class CollectTargetedPcrMetricsIntegrationTest extends CommandLineProgramTest {


    private static final File TEST_DATA_PATH = new File(getTestDataDir(), "picard/analysis/Collect_Targeted_PCR_Metrics/");


    @Test
    public void testCollect () throws IOException {


        final File input = new File(TEST_DATA_PATH, "first5Kamod_reordered.bam");
        final File amplicon_intervals = new File(TEST_DATA_PATH, "intervals_b37_20_1.interval_list");
        final File target_intervals = new File(TEST_DATA_PATH, "intervals_b37_20_1.interval_list");
        final File expectedFile = new File(TEST_DATA_PATH, "PCRMetrics_new_bam.txt");
        final File outfile = BaseTest.createTempFile("PCRMetrics", ".txt");
        final File reference = new File(largeFileTestDir, "human_g1k_v37.20.21.fasta");
        final File pertargetcoverage = new File(TEST_DATA_PATH, "pcr_metrics_pertarg_new.txt");



        final String[] args = {
                "--INPUT", input.getAbsolutePath(),
                "--AI", amplicon_intervals.getAbsolutePath(),
                "--TI", target_intervals.getAbsolutePath(),
                "--OUTPUT", outfile.getAbsolutePath(),
                "--R", reference.getAbsolutePath(),
                "--N", "intervals_b37_20",
                "--LEVEL", "ALL_READS",
                "--PER_TARGET_COVERAGE", pertargetcoverage.getAbsolutePath(),

        };

        runCommandLine(args);
        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");

    }
}
