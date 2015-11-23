package org.broadinstitute.hellbender.tools.picard.analysis.directed;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public final class CollectJumpingLibraryMetricsIntegrationTest extends CommandLineProgramTest {

    private static final File TEST_DATA_PATH = new File(getTestDataDir(), "picard/analysis/directed/CollectJumpingLibraryMetrics/");

    @Override
    public String getTestedClassName() {
        return CollectJumpingLibraryMetrics.class.getSimpleName();
    }

    @Test
    public void testCollect() throws IOException {

        final File input = new File(TEST_DATA_PATH, "first5000a.bam");
        final File expectedFile = new File(TEST_DATA_PATH, "JL_MQ39_T25k_CHIM200k.txt"); //done using picard 1.137(?)
        final File outfile = BaseTest.createTempFile("jumpinglib_1", ".txt");

        final String[] args = {
                "--input", input.getAbsolutePath(),
                "--output", outfile.getAbsolutePath(),
                "--MQ", "39",
                "--T", "25000",
                "--CHIMERA_KB_MIN", "200000"
        };

        runCommandLine(args);
        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
    }
}

