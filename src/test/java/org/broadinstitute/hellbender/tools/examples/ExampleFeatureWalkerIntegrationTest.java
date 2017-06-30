package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.TestResources;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

public final class ExampleFeatureWalkerIntegrationTest extends CommandLineProgramTest {

    private static final String TEST_DATA_DIRECTORY = TestResources.publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final String TEST_OUTPUT_DIRECTORY = TestResources.publicTestDir + "org/broadinstitute/hellbender/tools/examples/";

    @Test
    public void testExampleFeatureWalker() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + TestResources.hg19MiniReference +
                " -I " + TEST_DATA_DIRECTORY + "reads_data_source_test1.bam" +
                " -F " + TEST_DATA_DIRECTORY + "example_features.bed" +
                " -auxiliaryVariants " + TEST_DATA_DIRECTORY + "feature_data_source_test.vcf" +
                " -O %s",
                Arrays.asList(TEST_OUTPUT_DIRECTORY + "expected_ExampleFeatureWalkerIntegrationTest_output.txt")
        );
        testSpec.executeTest("testExampleIntervalWalker", this);
    }

    @Test
    public void testExampleFeatureWalkerWithIntervals() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + TestResources.hg19MiniReference +
                        " -I " + TEST_DATA_DIRECTORY + "reads_data_source_test1.bam" +
                        " -F " + TEST_DATA_DIRECTORY + "example_features.bed" +
                        " -L 1 " +
                        " -auxiliaryVariants " + TEST_DATA_DIRECTORY + "feature_data_source_test.vcf" +
                        " -O %s",
                Arrays.asList(TEST_OUTPUT_DIRECTORY + "expected_ExampleFeatureWalkerIntegrationTestWithIntervals_output.txt")
        );
        testSpec.executeTest("testExampleIntervalWalker", this);
    }
}
