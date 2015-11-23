package org.broadinstitute.hellbender.tools.picard.analysis.directed;


import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public final class CollectTargetedPcrMetricsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_PATH = new File(getTestDataDir(), "picard/analysis/directed/CollectTargetedPcrMetrics/");

    //Note: the expected values were generated using picard 1.130
    @Test
    public void testCollect() throws IOException {
        final File input = new File(TEST_DATA_PATH, "microbam.bam");
        final File amplicon_intervals = new File(TEST_DATA_PATH, "lifted_Chr20test_targets.interval_list");
        final File target_intervals = new File(TEST_DATA_PATH, "lifted_Chr20test_regions_t.interval_list");
        final File expectedFile = new File(TEST_DATA_PATH, "PCR_new_Metrics.txt");
        final File outfile = BaseTest.createTempFile("PCRMetrics", ".txt");
        final File reference = new File(largeFileTestDir, "human_g1k_v37.20.21.fasta");
        final File pertargetcoverage = new File(TEST_DATA_PATH, "pcr_metrics_pertarg_new.txt");
        final File pertargetcoverageout = BaseTest.createTempFile("pcr_metrics_pertarg_new_out", ".txt");

        final String[] args = {

                "--input", input.getAbsolutePath(),
                "--AI", amplicon_intervals.getAbsolutePath(),
                "--TI", target_intervals.getAbsolutePath(),
                "--output", outfile.getAbsolutePath(),
                "--R", reference.getAbsolutePath(),
                "--N", "lifted_Chr20test_targets",
                "--LEVEL", "ALL_READS",
                "--PER_TARGET_COVERAGE", pertargetcoverageout.getAbsolutePath(),
        };
        runCommandLine(args);
        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
        IntegrationTestSpec.assertEqualTextFiles(pertargetcoverageout, pertargetcoverage, "#");
    }
}
