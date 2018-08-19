package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class ExampleLocusWalkerIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_DATA_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final String TEST_OUTPUT_DIRECTORY = exampleTestDir;

    @Test
    public void testExampleLocusWalker() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
            " -L 1" +
                " -R " + hg19MiniReference +
                " -I " + TEST_DATA_DIRECTORY + "reads_data_source_test1.bam" +
                " -V " + TEST_DATA_DIRECTORY + "feature_data_source_test.vcf" +
                " -O %s",
            Arrays.asList(TEST_OUTPUT_DIRECTORY + "expected_ExampleLocusWalkerIntegrationTest_output.txt")
        );
        testSpec.executeTest("testExampleLocusWalker", this);
    }

}
