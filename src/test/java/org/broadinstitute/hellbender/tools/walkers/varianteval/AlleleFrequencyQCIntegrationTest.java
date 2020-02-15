package org.broadinstitute.hellbender.tools.walkers.varianteval;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

public class AlleleFrequencyQCIntegrationTest extends CommandLineProgramTest {

    private String evalVcf = getToolTestDataDir() + "af.na12878_array.vcf";
    private String comparisonVcf = getToolTestDataDir() + "af.thousand_genomes.10sites.vcf";


    private String getExpectedFile(String testName) {
        return getToolTestDataDir() + "expected/" + testName + ".expected.txt";
    }


    @Test(groups = {"R"})
    public void testAlleleFrequencyIntegrationTest() throws IOException {
        String name = "testAFQCIntegration";

        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + b37Reference +
                        " --eval " + evalVcf +
                        " --comp " + comparisonVcf +
                        " -eval:thousand_genomes " + comparisonVcf +
                        " -L " + comparisonVcf +
                        " -O %s"
                , Arrays.asList(getExpectedFile(name)));

        spec.executeTest(name, this);
    }

}

