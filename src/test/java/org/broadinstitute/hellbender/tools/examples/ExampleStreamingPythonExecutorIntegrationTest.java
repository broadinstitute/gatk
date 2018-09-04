package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

public class ExampleStreamingPythonExecutorIntegrationTest extends CommandLineProgramTest {

    private static final String TEST_DATA_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final String TEST_OUTPUT_DIRECTORY = exampleTestDir;

    @DataProvider(name = "streamingBatchSizes")
    public Object[][] getBatchSizes() {
        return new Object[][] {
                { 1 },
                { 2 },
                { 10 },
                { 11 },
                { 12 },
                { 100 }
        };
    }

    @Test(groups = "python", dataProvider = "streamingBatchSizes")
    public void testExampleStreamingPythonExecutor(final int batchSize) throws IOException {
        final IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -I " + TEST_DATA_DIRECTORY + "reads_data_source_test1.bam" +
                        " -O %s" +
                        " --batchSize " + Integer.toString(batchSize),
                Arrays.asList(TEST_OUTPUT_DIRECTORY + "expected_ExampleStreamingPythonExecutorIntegrationTest_output.txt")
        );
        testSpec.executeTest("testExampleStreamingPythonExecutor", this);
    }

}
