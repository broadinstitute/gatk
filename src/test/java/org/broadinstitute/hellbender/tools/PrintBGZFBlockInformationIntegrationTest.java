package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class PrintBGZFBlockInformationIntegrationTest extends CommandLineProgramTest {

    /* Well-formed large BGZF file */
    @Test
    public void testNormalLargeInput() throws IOException {
        final File input = new File(largeFileTestDir, "gvcfs/HG00096.g.vcf.gz");
        final File actualOutput = createTempFile("PrintBGZFBlockInformationIntegrationTest_testNormalLargeInput", ".out");
        final File expectedOutput = new File(toolsTestDir + "PrintBGZFBlockInformation/expected_PrintBGZFBlockInformationIntegrationTest_testNormalLargeInput.out");

        final String[] args = {
            "--bgzf-file", input.getAbsolutePath(),
            "--"  + StandardArgumentDefinitions.OUTPUT_LONG_NAME, actualOutput.getAbsolutePath()
        };
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(actualOutput, expectedOutput);
    }

    /* Well-formed small BGZF file */
    @Test
    public void testNormalSmallInput() throws IOException {
        final File input = new File(toolsTestDir, "PrintBGZFBlockInformation/4featuresHG38Header.vcf.gz");
        final File actualOutput = createTempFile("PrintBGZFBlockInformationIntegrationTest_testNormalSmallInput", ".out");
        final File expectedOutput = new File(toolsTestDir + "PrintBGZFBlockInformation/expected_PrintBGZFBlockInformationIntegrationTest_testNormalSmallInput.out");

        final String[] args = {
                "--bgzf-file", input.getAbsolutePath(),
                "--"  + StandardArgumentDefinitions.OUTPUT_LONG_NAME, actualOutput.getAbsolutePath()
        };
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(actualOutput, expectedOutput);
    }

    /* Malformed BGZF file missing the final 0-byte terminator block */
    @Test
    public void testMissingBGZFTerminatorBlock() throws IOException {
        final File input = new File(toolsTestDir, "PrintBGZFBlockInformation/4featuresHG38Header.NO_BGZF_TERMINATOR_BLOCK.vcf.gz");
        final File actualOutput = createTempFile("PrintBGZFBlockInformationIntegrationTest_testMissingBGZFTerminatorBlock", ".out");
        final File expectedOutput = new File(toolsTestDir + "PrintBGZFBlockInformation/expected_PrintBGZFBlockInformationIntegrationTest_testMissingBGZFTerminatorBlock.out");

        final String[] args = {
                "--bgzf-file", input.getAbsolutePath(),
                "--"  + StandardArgumentDefinitions.OUTPUT_LONG_NAME, actualOutput.getAbsolutePath()
        };
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(actualOutput, expectedOutput);
    }

    /* Malformed BGZF file with an incomplete (truncated) final block */
    @Test(expectedExceptions= UserException.CouldNotReadInputFile.class)
    public void testTruncatedFinalBlock() throws IOException {
        final File input = new File(toolsTestDir, "PrintBGZFBlockInformation/4featuresHG38Header.TRUNCATED_FINAL_BLOCK.vcf.gz");
        final File actualOutput = createTempFile("PrintBGZFBlockInformationIntegrationTest_testTruncatedFinalBlock", ".out");

        final String[] args = {
                "--bgzf-file", input.getAbsolutePath(),
                "--"  + StandardArgumentDefinitions.OUTPUT_LONG_NAME, actualOutput.getAbsolutePath()
        };
        runCommandLine(args);
    }

    /* Malformed BGZF file with an extra 0-byte terminator block in the middle */
    @Test
    public void testExtraTerminatorBlockInMiddle() throws IOException {
        final File input = new File(toolsTestDir, "PrintBGZFBlockInformation/4featuresHG38Header.EXTRA_TERMINATOR_BLOCK_IN_MIDDLE.vcf.gz");
        final File actualOutput = createTempFile("PrintBGZFBlockInformationIntegrationTest_testExtraTerminatorBlockInMiddle", ".out");
        final File expectedOutput = new File(toolsTestDir + "PrintBGZFBlockInformation/expected_PrintBGZFBlockInformationIntegrationTest_testExtraTerminatorBlockInMiddle.out");

        final String[] args = {
                "--bgzf-file", input.getAbsolutePath(),
                "--"  + StandardArgumentDefinitions.OUTPUT_LONG_NAME, actualOutput.getAbsolutePath()
        };
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(actualOutput, expectedOutput);
    }

    /* Regular GZIP file masquerading as a BGZF file */
    @Test(expectedExceptions= UserException.CouldNotReadInputFile.class)
    public void testRegularGzipFile() throws IOException {
        final File input = new File(toolsTestDir, "PrintBGZFBlockInformation/4featuresHG38Header.REGULAR_GZIP.vcf.gz");
        final File actualOutput = createTempFile("PrintBGZFBlockInformationIntegrationTest_testRegularGzipFile", ".out");

        final String[] args = {
                "--bgzf-file", input.getAbsolutePath(),
                "--"  + StandardArgumentDefinitions.OUTPUT_LONG_NAME, actualOutput.getAbsolutePath()
        };
        runCommandLine(args);
    }

    /* We should get an exception for other non-BGZF formats as well */
    @Test(expectedExceptions= UserException.CouldNotReadInputFile.class)
    public void testNonBGZFFile() throws IOException{
        final File input = new File(dbsnp_138_b37_1_65M_vcf);
        final File actualOutput = createTempFile("PrintBGZFBlockInformationIntegrationTest_testNonBGZFFile", ".out");

        final String[] args = {
                "--bgzf-file", input.getAbsolutePath(),
                "--"  + StandardArgumentDefinitions.OUTPUT_LONG_NAME, actualOutput.getAbsolutePath()
        };
        runCommandLine(args);
    }

    /* Make sure that we can handle a standard BAM file */
    @Test
    public void testBamFile() throws IOException {
        final File input = new File(packageRootTestDir + "engine/reads_data_source_test1.bam");
        final File actualOutput = createTempFile("PrintBGZFBlockInformationIntegrationTest_testBamFile", ".out");
        final File expectedOutput = new File(toolsTestDir + "PrintBGZFBlockInformation/expected_PrintBGZFBlockInformationIntegrationTest_testBamFile.out");

        final String[] args = {
                "--bgzf-file", input.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME, actualOutput.getAbsolutePath()
        };
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(actualOutput, expectedOutput);
    }
}
