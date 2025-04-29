package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
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

    @Test
    public void testInGatkLiteDocker() throws IOException {
        final String gatkLiteDockerProperty = System.getProperty("IN_GATKLITE_DOCKER");

        try {
            System.setProperty("IN_GATKLITE_DOCKER", "true");
            final IntegrationTestSpec testSpec = new IntegrationTestSpec(
                    " -I " + TEST_DATA_DIRECTORY + "reads_data_source_test1.bam" +
                            " -O %s" +
                            " --batchSize " + Integer.toString(2),
                    Arrays.asList(TEST_OUTPUT_DIRECTORY + "expected_ExampleStreamingPythonExecutorIntegrationTest_output.txt")
            );
            
            try {
                testSpec.executeTest("testExampleStreamingPythonExecutor", this);
                Assert.fail("Excepted RuntimeException for running in GATK Lite docker");
            }
            catch(final RuntimeException e) {
                Assert.assertTrue(e.getMessage().contains("Tools requiring Python cannot be run in the Gatk Lite docker image.")); 
            }

        }
        finally {
            if(gatkLiteDockerProperty != null) {
                System.setProperty("IN_GATKLITE_DOCKER", gatkLiteDockerProperty);
            }
            else{
                System.clearProperty("IN_GATKLITE_DOCKER");
            } 
        }
    }

}
