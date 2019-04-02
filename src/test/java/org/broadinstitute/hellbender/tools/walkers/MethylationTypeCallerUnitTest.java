package org.broadinstitute.hellbender.tools.walkers;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;

public class MethylationTypeCallerUnitTest extends CommandLineProgramTest {

    // If true, update the expected outputs in tests that assert an exact match vs. prior output,
    // instead of actually running the tests. Can be used with "./gradlew test -Dtest.single=MethylationTypeCallerUnitTest"
    // to update all of the exact-match tests at once. After you do this, you should look at the
    // diffs in the new expected outputs in git to confirm that they are consistent with expectations.
    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;

    public static final String TEST_FILES_DIR = toolsTestDir + "MethylationTypeCaller/";

    /*
     * Test that in VCF mode we're consistent with past  results
     */
    @Test
    public void testVCFIsConsistentWithPastResults(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File outputVCF = createTempFile("testVCFIsConsistentWithPastResults", ".vcf");
        final File outputVCFIndex = createTempFile("testVCFIndexIsConsistentWithPastResults", ".vcf.idx");
        final File expectedVCF = new File(TEST_FILES_DIR, "chr14.unique_reads.methylC_seq.vcf");
        final File expectedVCFIndex = new File(TEST_FILES_DIR, "chr14.unique_reads.methylC_seq.vcf.idx");

        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedVCF.getAbsolutePath() : outputVCF.getAbsolutePath();

        final String[] inputArgs = {
                "-I", TEST_FILES_DIR + inputFileName,
                "-R", toolsTestDir +  referenceFileName,
                "-O", outputPath
        };

        runCommandLine(inputArgs);

        // Test for an exact match against past results
        if (!UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS) {
            IntegrationTestSpec.assertEqualTextFiles(outputVCF, expectedVCF);
            IntegrationTestSpec.assertEqualTextFiles(outputVCFIndex, expectedVCFIndex);
        }
    }
}
