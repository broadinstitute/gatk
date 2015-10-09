package org.broadinstitute.hellbender.tools.picard.analysis.directed;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public final class CollectWgsMetricsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/analysis/directed/CollectWgsMetrics");

    @Test
    public void test() throws IOException {
        final File input = new File(NA12878_20_21_WGS_bam);
        final File refFile = new File(b37_reference_20_21);
        final File expectedFile = new File(TEST_DATA_DIR, "CollectWgsMetrics.txt");
        final File outfile = BaseTest.createTempFile("testCollectWgsMetrics", ".metrics");
        final String[] args = new String[]{
                "--INPUT", input.getAbsolutePath(),
                "--OUTPUT", outfile.getAbsolutePath(),
                "--REFERENCE", refFile.getAbsolutePath(),
                "--VALIDATION_STRINGENCY", "LENIENT",
                "--STOP_AFTER", "10000000",
                "--INCLUDE_BQ_HISTOGRAM", "TRUE"
        };
        runCommandLine(args);
        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
    }
}