package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.NativeUtils;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.Test;
import org.broadinstitute.hellbender.tools.walkers.variantrecalling.FlowTestConstants;

import java.io.File;
import java.io.IOException;

public class FlowBasedAlignmentIntegrationTest extends CommandLineProgramTest {

    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;
    private static String    testDir = publicTestDir + "/large";

    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled");
    }
    @Test
    public void testMatrix() throws IOException {
        // Compares the flow-based likelihood matrix exactly; these doubles differ by ~1 ULP across
        // CPU architectures (Java Math.* is not bit-identical across platforms), so the x86-generated
        // reference does not bit-match on other architectures. Run on x86 only.
        if (!NativeUtils.runningOn64BitX86Architecture()) {
            throw new SkipException("Flow-based likelihood matrix is x86-reference-specific (Java Math.* FP differs ~1 ULP across architectures).");
        }

        final File outputDir = createTempDir("testMatrix");
        final File expectedFile = new File(publicTestDir + "/" + FlowTestConstants.FLOW_BASED_ALIGNMENT_DATA_DIR + "/input_jukebox_for_test.expected.alm");
        final File outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : new File(outputDir + "/output.alm");

        final String[] args = new String[] {
                "-R", publicTestDir + "/large/Homo_sapiens_assembly38.fasta.gz",
                "-O", outputDir + "/ignored.vcf",
                "-I", testDir + "/input_jukebox_for_test.bam",
                "--smith-waterman", "FASTEST_AVAILABLE",
                "--likelihood-calculation-engine", "FlowBased",
                "-mbq", "0",
                "--kmer-size", "10",
                "--intervals", "chr9:81148694-81177540",
                "--alm-path", outputFile.getAbsolutePath(),
                "--alm-interval", "chr9"
        };

        runCommandLine(args);  // no assert, just make sure we don't throw

        // verify that output file has been created
        Assert.assertTrue(outputFile.exists());
        if ( !UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedFile);
        }
    }

    @Override
    public String getTestedToolName() {
        return "HaplotypeCaller";
    }

}
