package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.util.Arrays;

public class ExampleAssemblyRegionWalkerIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_OUTPUT_DIRECTORY = exampleTestDir;

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
