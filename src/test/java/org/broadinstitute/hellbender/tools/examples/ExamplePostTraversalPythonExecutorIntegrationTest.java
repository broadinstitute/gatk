package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

public class ExamplePostTraversalPythonExecutorIntegrationTest extends CommandLineProgramTest {

    private static final String TEST_DATA_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final String TEST_OUTPUT_DIRECTORY = exampleTestDir;

    @Test
    public void testExamplePythonExecutor() throws IOException {
        final IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + hg19MiniReference +
                        " -I " + TEST_DATA_DIRECTORY + "reads_data_source_test1.bam" +
                        " -O %s" +
                        " -pythonOutputFile %s",
                Arrays.asList(
                        TEST_OUTPUT_DIRECTORY + "expected_ExamplePythonExecutorIntegrationTest_Reads_output.txt",
                        TEST_OUTPUT_DIRECTORY + "expected_ExamplePythonExecutorIntegrationTest_Python_output.txt"
                )
        );
        testSpec.executeTest("testExamplePythonExecutor", this);
    }

    @Test
    public void testExamplePythonExecutorInGatkLiteDocker() throws IOException {
        final String gatkLiteDockerProperty = System.getProperty("IN_GATKLITE_DOCKER");

        try {
            System.setProperty("IN_GATKLITE_DOCKER", "true");
            final IntegrationTestSpec testSpec = new IntegrationTestSpec(
                    " -R " + hg19MiniReference +
                            " -I " + TEST_DATA_DIRECTORY + "reads_data_source_test1.bam" +
                            " -O %s" +
                            " -pythonOutputFile %s",
                    Arrays.asList(
                            TEST_OUTPUT_DIRECTORY + "expected_ExamplePythonExecutorIntegrationTest_Reads_output.txt",
                            TEST_OUTPUT_DIRECTORY + "expected_ExamplePythonExecutorIntegrationTest_Python_output.txt"
                    )
            );
            try {
                testSpec.executeTest("testExamplePythonExecutor", this);
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
