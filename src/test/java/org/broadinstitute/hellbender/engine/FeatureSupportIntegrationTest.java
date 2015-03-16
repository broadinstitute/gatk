package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

public class FeatureSupportIntegrationTest extends CommandLineProgramTest {
    private static final String FEATURE_INTEGRATION_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";

    @Override
    public String getTestedClassName() {
        return "PrintReadsWithVariants";
    }

    @Test
    public void testFeatureSupportUsingVCF() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + hg19MiniReference +
                " -I " + FEATURE_INTEGRATION_TEST_DIRECTORY + "reads_data_source_test1.bam" +
                " -V " + "TestFeatures:" + FEATURE_INTEGRATION_TEST_DIRECTORY + "feature_data_source_test.vcf" +
                " --groupVariantsBySource" +
                " -O %s",
                Arrays.asList(FEATURE_INTEGRATION_TEST_DIRECTORY + "expected_testFeatureSupportUsingVCF_output.txt")
        );
        testSpec.executeTest("testFeatureSupportUsingVCF", this);
    }
}
