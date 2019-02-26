package org.broadinstitute.hellbender.tools.walkers.bqsr;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;

public final class BaseRecalibratorIntegrationTest extends AbstractBaseRecalibratorIntegrationTest {
    @Test
    public void testBQSRFailWithIncompatibleSequenceDictionaries() throws IOException {
        final String resourceDir =  getTestDataDir() + "/" + "BQSR" + "/";

        final String bam_chr20 = resourceDir + WGS_B37_CH20_1M_1M1K_BAM;
        final String dbSNPb37_chr17 =  getResourceDir() + "dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf";

        final String  NO_ARGS = "";
        final BQSRTest params = new BQSRTest(b37_reference_20_21, bam_chr20, dbSNPb37_chr17, NO_ARGS, resourceDir + "expected.txt");
        IntegrationTestSpec spec = new IntegrationTestSpec(
                params.getCommandLine(),
                1,
                UserException.IncompatibleSequenceDictionaries.class);
        spec.executeTest("testBQSRFailWithIncompatibleSequenceDictionaries", this);
    }
}
