package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.annotations.Test;

import java.io.IOException;

/**
 * Created by sam on 2/8/18.
 */
public class VariantTranchesFromInfoKeyIntegrationTest extends CommandLineProgramTest {
    final private static String inputVCF = largeFileTestDir + "VQSR/g94982_b37_chr20_1m_10m.vcf.gz";
    final private static String snpTruthVCF = largeFileTestDir + "VQSR/giab_chr20_1m_10m.vcf.gz";
    final private static String indelTruthVCF = largeFileTestDir + "VQSR/giab_chr20_1m_10m.vcf.gz";
    final private static String outputVCF = largeFileTestDir + "VQSR/expected/variant_tranches_python_expected.vcf";

    @Test(groups = {"python"})
    public void testTranches() throws IOException {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, outputVCF)
                .addArgument("snp-truth-vcf", snpTruthVCF)
                .addArgument("indel-truth-vcf", indelTruthVCF)
                .addArgument("tranche", "99.0")
                .addArgument("tranche", "95.0")
                .addArgument("max-sites", "2000")
                .addArgument("info-key", "VQSLOD");

        runCommandLine(argsBuilder);
    }
}
