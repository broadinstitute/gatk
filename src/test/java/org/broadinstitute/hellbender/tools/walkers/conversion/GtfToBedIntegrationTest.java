package org.broadinstitute.hellbender.tools.walkers.conversion;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.annotations.Test;

import java.io.File;

public class GtfToBedIntegrationTest extends CommandLineProgramTest {
    //private static final File TEST_SUB_DIR = new File(toolsTestDir);
    //private static final File GENCODE_PIK3CA_GTF = new File(TEST_SUB_DIR, "/funcotator/small_pik3ca_dbsnp_ds/gencode_pik3ca/hg19/gencode.v19.PIK3CA.gtf");
    private static final File input = new File("/Users/shahsana/TestGatk/gencode.v38.chr_patch_hapl_scaff.annotation.gtf");
    @Test
    public void testDoesNotCrash() {
        //final File outputFile = createTempFile("gtf_to_bed_test", ".bed");
        File outputFile = new File("/Users/shahsana/TestGatk/output.bed");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                //.add("G", GENCODE_PIK3CA_GTF)
                .add("G", input)
                .add(GtfToBed.SORT_BY_TRANSCRIPT_LONG_NAME, true)
                .addOutput(outputFile);
        runCommandLine(argsBuilder);
    }
}
