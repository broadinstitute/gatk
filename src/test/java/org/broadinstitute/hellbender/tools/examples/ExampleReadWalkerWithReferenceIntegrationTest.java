package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

public final class ExampleReadWalkerWithReferenceIntegrationTest extends CommandLineProgramTest {

    private static final String TEST_DATA_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final String TEST_OUTPUT_DIRECTORY = exampleTestDir;

    @Test
    public void testExampleReadWalkerWithReference() throws IOException {
        final IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + hg19MiniReference +
                " -I " + TEST_DATA_DIRECTORY + "reads_data_source_test1.bam" +
                " -O %s",
                Arrays.asList(TEST_OUTPUT_DIRECTORY + "expected_ExampleReadWalkerWithReferenceIntegrationTest_output.txt")
        );
        testSpec.executeTest("testExampleReadWalkerWithReference", this);
    }

}
