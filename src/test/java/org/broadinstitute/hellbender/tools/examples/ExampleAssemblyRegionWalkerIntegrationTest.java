package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;

public class ExampleAssemblyRegionWalkerIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_OUTPUT_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/tools/examples/";

    @Test
    public void testExampleAssemblyRegionWalker() throws Exception {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + b37_reference_20_21 +
                " -I " + NA12878_20_21_WGS_bam +
                " -knownVariants " + dbsnp_138_b37_20_21_vcf +
                " -L 20:10000000-10050000 " +
                " -O %s",
                Arrays.asList(TEST_OUTPUT_DIRECTORY + "expected_ExampleAssemblyRegionWalkerIntegrationTest_output.txt")
        );

        testSpec.executeTest("testExampleAssemblyRegionWalker", this);
    }
}
