package org.broadinstitute.hellbender.tools.picard.analysis.directed;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;
import java.io.File;
import java.io.IOException;

/**
 * Created by dkling on 9/1/15.
 */
public final class CollectTargetedPcrMetrics_IT extends CommandLineProgramTest {

    private static final File TEST_DATA_PATH = new File(getTestDataDir(), "picard/analysis/Collect_Targeted_PCR_Metrics" );

    public String getTestedClassName() {
        return CollectTargetedPcrMetrics.class.getSimpleName();
    }

    @Test
    public void testCollect () throws IOException {

        final File input = new File(TEST_DATA_PATH, "first5000a.bam");
        final File amplicon_intervals = new File(TEST_DATA_PATH, "intervallist.vcf");
        final File target_intervals = new File(TEST_DATA_PATH, "intervallist.vcf");
        final File expectedFile = new File(TEST_DATA_PATH, "PCRMetrics_2.txt");
        final File outfile = BaseTest.createTempFile("PCRMetrics", ".txt");
        final File reference = new File(TEST_DATA_PATH, "human_b37_20.fasta");

        final String[] args = {
                "--INPUT", input.getAbsolutePath(),
                "--AI", amplicon_intervals.getAbsolutePath(),
                "--TI", target_intervals.getAbsolutePath(),
                "--OUTPUT", outfile.getAbsolutePath(),
                "--R", reference.getAbsolutePath() ,
        };

        runCommandLine(args);
        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");

    }

}
