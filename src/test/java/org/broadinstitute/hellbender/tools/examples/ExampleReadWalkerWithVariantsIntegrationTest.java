package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

public final class ExampleReadWalkerWithVariantsIntegrationTest extends CommandLineProgramTest {

    private static final String TEST_DATA_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final String TEST_OUTPUT_DIRECTORY = exampleTestDir;

    @Test
    public void testExampleReadWalkerWithVariants() throws IOException {
        final IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + hg19MiniReference +
                " -I " + TEST_DATA_DIRECTORY + "reads_data_source_test1.bam" +
                " -V " + TEST_DATA_DIRECTORY + "example_variants_withSequenceDict.vcf" +
                " -O %s",
                Arrays.asList(TEST_OUTPUT_DIRECTORY + "expected_ExampleReadWalkerWithVariantsIntegrationTest_output.txt")
        );
        testSpec.executeTest("testExampleReadWalkerWithVariants", this);
    }

}
