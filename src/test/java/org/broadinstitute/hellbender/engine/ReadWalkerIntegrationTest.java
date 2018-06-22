package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.examples.ExampleReadWalkerWithReference;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class ReadWalkerIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return ExampleReadWalkerWithReference.class.getSimpleName();
    }

    @Test
    public void testManuallySpecifiedIndices() throws IOException {
        final String BAM_PATH = publicTestDir + "org/broadinstitute/hellbender/engine/readIndexTest/";
        final String INDEX_PATH = BAM_PATH + "indices/";
        final File outFile = createTempFile("testManuallySpecifiedIndices", ".txt");
        final File expectedFile = new File(publicTestDir + "org/broadinstitute/hellbender/engine/expected_ReadWalkerIntegrationTest_testManuallySpecifiedIndices.txt");

        final String[] args = new String[] {
            "-I", BAM_PATH + "reads_data_source_test1.bam",
            "-I", BAM_PATH + "reads_data_source_test2.bam",
            "--read-index", INDEX_PATH + "reads_data_source_test1.bam.bai",
            "--read-index", INDEX_PATH + "reads_data_source_test2.bam.bai",
            "-O", outFile.getAbsolutePath()
        };
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(outFile, expectedFile);
    }

    @Test(expectedExceptions = UserException.class)
    public void testManuallySpecifiedIndicesWrongNumberOfIndices() throws IOException {
        final String BAM_PATH = publicTestDir + "org/broadinstitute/hellbender/engine/readIndexTest/";
        final String INDEX_PATH = BAM_PATH + "indices/";

        final String[] args = new String[] {
                "-I", BAM_PATH + "reads_data_source_test1.bam",
                "-I", BAM_PATH + "reads_data_source_test2.bam",
                "--read-index", INDEX_PATH + "reads_data_source_test1.bam.bai"
        };
        runCommandLine(args);
    }
}
