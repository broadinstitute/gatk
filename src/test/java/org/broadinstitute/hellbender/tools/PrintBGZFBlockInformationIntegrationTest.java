package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class PrintBGZFBlockInformationIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testNormalInput() throws IOException {
        final File input = new File(largeFileTestDir, "gvcfs/HG00096.g.vcf.gz");
        final File actualOutput = createTempFile("PrintBGZFBlockInformationIntegrationTest_testNormalInput", ".out");
        final File expectedOutput = new File(toolsTestDir + "diagnostics/expected_PrintBGZFBlockInformationIntegrationTest_testNormalInput.out");

        final String[] args = {
            "--bgzf-file", input.getAbsolutePath(),
            "--"  + StandardArgumentDefinitions.OUTPUT_LONG_NAME, actualOutput.getAbsolutePath()
        };
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(actualOutput, expectedOutput);
    }
}
