package org.broadinstitute.hellbender.tools.picard.analysis;

/**
 * Created by dkling on 8/26/15.
 */

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.picard.analysis.directed.CollectTargetedPcrMetrics;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;


public final class CollectTargetedPcrMetrics_IT extends CommandLineProgramTest {

    private static final File TEST_DATA_PATH = new File(getTestDataDir(), "picard/analysis/" );

    public String getTestedClassName() {
        return CollectTargetedPcrMetrics.class.getSimpleName();
    }

    @Test
    public void testCollect () throws IOException {

        final File input = new File(TEST_DATA_PATH, "first5000a.bam");
        final File amplicon_intervals = new File(TEST_DATA_PATH, "intervallist_copy.vcf");
        final File target_intervals = new File(TEST_DATA_PATH, "intervallist_copy.vcf");
        final File expectedFile = new File(TEST_DATA_PATH, "PCRMetrics.txt");
        final File outfile = BaseTest.createTempFile("PCRmetrics", ".txt");

        final String[] args = {
                "--INPUT", input.getAbsolutePath(),
                "--Amplicon_Intervals", amplicon_intervals.getAbsolutePath(),
        "--Target_Intervals", target_intervals.getAbsolutePath(),
                "--OUTPUT", outfile.getAbsolutePath(),
        };

        runCommandLine(args);
        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");

    }

}

