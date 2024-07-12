package org.broadinstitute.hellbender.tools.walkers.conversion;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.annotations.Test;

import java.io.File;

public class GtfToBedIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir);
    private static final File GENCODE_PIK3CA_GTF = new File(TEST_SUB_DIR, "/funcotator/small_pik3ca_dbsnp_ds/gencode_pik3ca/hg19/gencode.v19.PIK3CA.gtf");

    @Test
    public void testDoesNotCrash() {
        final File outputFile = createTempFile("gtf_to_bed_test", ".bed");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .add("G", GENCODE_PIK3CA_GTF)
                .addOutput(outputFile);
        runCommandLine(argsBuilder);
    }
}
