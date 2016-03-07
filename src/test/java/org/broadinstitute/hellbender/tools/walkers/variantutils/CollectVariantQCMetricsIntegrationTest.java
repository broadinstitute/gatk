package org.broadinstitute.hellbender.tools.walkers.variantutils;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

/**
 */
public class CollectVariantQCMetricsIntegrationTest extends CommandLineProgramTest {
    @Test
    public void test1() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " --variant " + getToolTestDataDir() + "CEUtrioTest.vcf" +
                        " -O %s",
                Arrays.asList(getToolTestDataDir() + "expected.CEUtrioTest.variantQC.table"));
        spec.executeTest("test1", this);
    }

    @Test(enabled = false)//disable for now until we are done with all metrics
    public void testLarge() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " --variant " + largeFileTestDir + "1000G.phase3.broad.withGenotypes.chr20.10100000.vcf" +
                        " -O %s",
                Arrays.asList(getToolTestDataDir() + "expected.CEUtrioTest.variantQC.table"));
        spec.executeTest("test1", this);
    }
}
