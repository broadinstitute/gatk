package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.tools.PrintReads;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Tests to prove that we can access and query inputs on Google Cloud Storage (GCS) in ReadWalkers.
 */
public class ReadWalkerGCSSupportIntegrationTest extends CommandLineProgramTest {

    private static final String TEST_BAM_ON_GCS = "org/broadinstitute/hellbender/engine/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.10000000-10000020.with.unmapped.bam";
    private static final String EXPECTED_OUTPUT_DIR = publicTestDir + "org/broadinstitute/hellbender/engine/GCSTests/";

    @Override
    public String getTestedToolName() {
        return PrintReads.class.getSimpleName();
    }

    @DataProvider(name = "GCSTestCases")
    public Object[][] getGCSTestCases() {
        final SimpleInterval singleInterval = new SimpleInterval("20", 10000010, 10000016);
        final List<SimpleInterval> multipleIntervals = Arrays.asList(singleInterval, new SimpleInterval("21", 10000005, 10000011));

        final String EXPECTED_WHOLE_FILE_RESULTS = EXPECTED_OUTPUT_DIR + "expected_ReadWalkerGCSSupportIntegrationTest_bam_wholefile.bam";
        final String EXPECTED_SINGLE_INTERVAL_RESULTS = EXPECTED_OUTPUT_DIR + "expected_ReadWalkerGCSSupportIntegrationTest_bam_single_interval.bam";
        final String EXPECTED_MULTIPLE_INTERVALS_RESULTS = EXPECTED_OUTPUT_DIR + "expected_ReadWalkerGCSSupportIntegrationTest_bam_multiple_intervals.bam";
        final String EXPECTED_MULTIPLE_INTERVALS_WITH_UNMAPPED_RESULTS = EXPECTED_OUTPUT_DIR + "expected_ReadWalkerGCSSupportIntegrationTest_bam_multiple_intervals_with_unmapped.bam";
        final String EXPECTED_UNMAPPED_ONLY_RESULTS = EXPECTED_OUTPUT_DIR + "expected_ReadWalkerGCSSupportIntegrationTest_bam_unmapped_only.bam";

        return new Object[][] {
                { TEST_BAM_ON_GCS, null, false, EXPECTED_WHOLE_FILE_RESULTS },
                { TEST_BAM_ON_GCS, Collections.singletonList(singleInterval), false, EXPECTED_SINGLE_INTERVAL_RESULTS },
                { TEST_BAM_ON_GCS, multipleIntervals, false, EXPECTED_MULTIPLE_INTERVALS_RESULTS },
                { TEST_BAM_ON_GCS, multipleIntervals, true, EXPECTED_MULTIPLE_INTERVALS_WITH_UNMAPPED_RESULTS },
                { TEST_BAM_ON_GCS, null, true, EXPECTED_UNMAPPED_ONLY_RESULTS }
        };
    }

    @Test(dataProvider = "GCSTestCases", groups = {"bucket"})
    public void testReadBAMOnGCS( final String bam, final List<SimpleInterval> intervals, final boolean includeUnmapped, final String expectedOutput ) throws IOException {
        final StringBuilder intervalArgBuilder = new StringBuilder("");
        if ( intervals != null ) {
            for ( final SimpleInterval interval : intervals ) {
                intervalArgBuilder.append(" -L ");
                intervalArgBuilder.append(interval.toString());
            }
        }
        if ( includeUnmapped ) {
            intervalArgBuilder.append(" -L ");
            intervalArgBuilder.append("unmapped");
        }
        String intervalArg = intervalArgBuilder.toString();
        
        final IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -I " + getGCPTestInputPath() + bam +
                intervalArg +
                " -O %s",
                Collections.singletonList(expectedOutput)
        );
        testSpec.executeTest("testReadBAMOnGCS", this);
    }

    @Test(groups = {"bucket"})
    public void testBAMReferenceAndIntervalsOnGCS() throws Exception {
        final File output = createTempFile("testBAMReferenceAndIntervalsOnGCS", ".out");
        final File expected = new File(EXPECTED_OUTPUT_DIR, "expected_testBAMReferenceAndIntervalsOnGCS.out");

        final String[] args = new String[] {
            "ExampleReadWalkerWithReference",
            "-I", getGCPTestInputPath() + "large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam",
            "-R", getGCPTestInputPath() + "large/human_g1k_v37.20.21.fasta",
            "-L", getGCPTestInputPath() + "large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.smallIntervalList.intervals",
            "-O", output.getAbsolutePath()
        };
        new Main().instanceMain(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expected);
    }
}
