package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

public final class ExampleIntervalWalkerIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_DATA_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final String TEST_OUTPUT_DIRECTORY = exampleTestDir;

    @Test
    public void testExampleIntervalWalker() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -L 1:100-200 -L 2:500-600" +
                " -R " + hg19MiniReference +
                " -I " + TEST_DATA_DIRECTORY + "reads_data_source_test1.bam" +
                " -V " + TEST_DATA_DIRECTORY + "feature_data_source_test.vcf" +
                " -O %s",
                Arrays.asList(TEST_OUTPUT_DIRECTORY + "expected_ExampleIntervalWalkerIntegrationTest_output.txt")
        );
        testSpec.executeTest("testExampleIntervalWalker", this);
    }
}
