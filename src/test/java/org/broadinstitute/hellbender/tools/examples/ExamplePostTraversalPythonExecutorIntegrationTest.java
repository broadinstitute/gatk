package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.EnvironmentTestUtils;
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

    @Test(singleThreaded = true)
    public void testExamplePythonExecutorInGatkLiteDocker() {
        EnvironmentTestUtils.checkWithGATKDockerPropertySet(() -> {
            final IntegrationTestSpec testSpec = new IntegrationTestSpec(
                    " -R " + hg19MiniReference +
                            " -I " + TEST_DATA_DIRECTORY + "reads_data_source_test1.bam" +
                            " -O " + TEST_OUTPUT_DIRECTORY + "expected_ExamplePythonExecutorIntegrationTest_Reads_output.txt" +
                            " -pythonOutputFile " + TEST_OUTPUT_DIRECTORY + "expected_ExamplePythonExecutorIntegrationTest_Python_output.txt",
                    0, UserException.NotAvailableInGatkLiteDocker.class
            );
            try {
                testSpec.executeTest("testExamplePythonExecutorInGatkLiteDocker", this);
            }
            catch (final IOException e) {
                Assert.fail("Failed with IOException: " + e.getMessage());
            }
        });
    }

}
