package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.examples.ExamplePartialReadWalker;
import org.broadinstitute.hellbender.tools.examples.ExampleReadWalkerWithReference;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class PartialReadWalkerIntegrationTest extends CommandLineProgramTest {

    // If true, update the expected outputs in tests that assert an exact match vs. prior output,
    // instead of actually running the tests. Can be used with "./gradlew test -Dtest.single=HaplotypeCallerIntegrationTest"
    // to update all of the exact-match tests at once. After you do this, you should look at the
    // diffs in the new expected outputs in git to confirm that they are consistent with expectations.
    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;

    /*
     * Make sure that someone didn't leave the UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS toggle turned on
     */
    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled");
    }

    @Override
    public String getTestedClassName() {
        return ExamplePartialReadWalker.class.getSimpleName();
    }

    @Test
    public void testPartialReadWalker() throws IOException {
        final String BAM_PATH = publicTestDir + "org/broadinstitute/hellbender/engine/readIndexTest/";
        final String INDEX_PATH = BAM_PATH + "indices/";
        final File expectedFile = new File(publicTestDir + "org/broadinstitute/hellbender/engine/expected_PartialReadWalkerIntegrationTest_testPartialReadWalker.txt");
        final File outFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : createTempFile("testPartialReadWalker", ".txt");

        final String[] args = new String[] {
            "-I", BAM_PATH + "reads_data_source_test1.bam",
            "--read-index", INDEX_PATH + "reads_data_source_test1.bam.bai",
            "-O", outFile.getAbsolutePath(),
            "--stop-on-read-name", "f"
        };
        runCommandLine(args);

        if ( !UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(outFile, expectedFile);
        }
    }

}
