package org.broadinstitute.hellbender.tools.walkers.sv;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.util.Arrays;

public class PairedEndAndSplitReadEvidenceCollectionIntegrationTest extends CommandLineProgramTest {

    private final String TEST_OUTPUT_DIRECTORY = getToolTestDataDir().toLowerCase() + "/";
    public static final String pesrTestDir = toolsTestDir + "walkers/sv/pesr";

    @Test
    public void testPESRCollection() throws Exception {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                "-I " + NA12878_20_21_WGS_bam + " --sample-name NA12878 -p %s -s %s",
                Arrays.asList(pesrTestDir + "/NA12878.disc.txt.gz", pesrTestDir + "/NA12878.split.txt.gz"));
        spec.executeTest("base PESR collection", this);
    }

}