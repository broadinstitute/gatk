package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class ExampleReadSliderWalkerIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_DATA_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final String TEST_OUTPUT_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/tools/examples/";

    @Test
    public void testExampleReadSliderWalker() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
            " -L 1:200-1125" + // region with variants/reads
            " -L 4:15951-16000" + // region with reference bases
                " -R " + hg19MiniReference +
                " -I " + TEST_DATA_DIRECTORY + "reads_data_source_test1.bam" +
                " -V " + TEST_DATA_DIRECTORY + "feature_data_source_test.vcf" +
                " -O %s",
            Arrays.asList(TEST_OUTPUT_DIRECTORY + "expected_ExampleReadSliderWalkerIntegrationTest_output.txt")
        );
        testSpec.executeTest("testExampleReadSliderWalker", this);
    }

}
