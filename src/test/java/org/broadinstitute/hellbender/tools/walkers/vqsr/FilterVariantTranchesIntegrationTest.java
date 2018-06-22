package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;

import org.testng.annotations.Test;

import java.util.Arrays;
import java.io.IOException;

public class FilterVariantTranchesIntegrationTest  extends CommandLineProgramTest {


    /**
     * Run the tool on a small test VCF.
     */
    @Test
    public void testTrancheFiltering() throws IOException {
        final String trancheVCF = largeFileTestDir + "VQSR/g94982_20_1m_10m_python_2dcnn.vcf.gz";
        final String indelTruthVCF = largeFileTestDir + "VQSR/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.sites.20.1M-10M.vcf";
        final String snpTruthVCF = largeFileTestDir + "VQSR/Omni25_sites_1525_samples.b37.20.1M-10M.vcf";

        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, trancheVCF)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s")
                .addArgument("resource", snpTruthVCF)
                .addArgument("resource", indelTruthVCF)
                .addArgument("snp-tranche", "99.5")
                .addArgument("indel-tranche", "99.0")
                .addArgument("invalidate-previous-filters", "true")
                .addArgument("info-key", "MIX_SMALL_2D_W_DROPOUT")
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");


        final IntegrationTestSpec spec = new IntegrationTestSpec(argsBuilder.toString(),
                    Arrays.asList(largeFileTestDir + "VQSR/expected/g94982_20_1m_10m_tranched_99_99.5.vcf"));
        spec.executeTest("testTrancheFiltering", this);
    }

    @Test
    public void testTrancheFilteringDuplicate() throws IOException {
        final String trancheVCF = largeFileTestDir + "VQSR/g94982_20_1m_10m_python_2dcnn.vcf.gz";
        final String indelTruthVCF = largeFileTestDir + "VQSR/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.sites.20.1M-10M.vcf";
        final String snpTruthVCF = largeFileTestDir + "VQSR/Omni25_sites_1525_samples.b37.20.1M-10M.vcf";

        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, trancheVCF)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s")
                .addArgument("resource", snpTruthVCF)
                .addArgument("resource", indelTruthVCF)
                .addArgument("snp-tranche", "99.5")
                .addArgument("indel-tranche", "99.0")
                .addArgument("snp-tranche", "99.5")
                .addArgument("indel-tranche", "99.0")
                .addArgument("invalidate-previous-filters", "true")
                .addArgument("info-key", "MIX_SMALL_2D_W_DROPOUT")
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        final IntegrationTestSpec spec = new IntegrationTestSpec(argsBuilder.toString(),
                Arrays.asList(largeFileTestDir + "VQSR/expected/g94982_20_1m_10m_tranched_99_99.5.vcf"));
        spec.executeTest("testTrancheFilteringDuplicate", this);
    }

    @Test
    public void testTrancheFilteringTranches() throws IOException {
        final String trancheVCF = largeFileTestDir + "VQSR/g94982_20_1m_10m_python_2dcnn.vcf.gz";
        final String indelTruthVCF = largeFileTestDir + "VQSR/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.sites.20.1M-10M.vcf";
        final String snpTruthVCF = largeFileTestDir + "VQSR/Omni25_sites_1525_samples.b37.20.1M-10M.vcf";

        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, trancheVCF)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s")
                .addArgument("resource", snpTruthVCF)
                .addArgument("resource", indelTruthVCF)
                .addArgument("snp-tranche", "95.0")
                .addArgument("snp-tranche", "99.0")
                .addArgument("snp-tranche", "99.9")
                .addArgument("indel-tranche", "95.0")
                .addArgument("indel-tranche", "99.0")
                .addArgument("indel-tranche", "99.9")
                .addArgument("invalidate-previous-filters", "true")
                .addArgument("info-key", "MIX_SMALL_2D_W_DROPOUT")
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        final IntegrationTestSpec spec = new IntegrationTestSpec(argsBuilder.toString(),
                Arrays.asList(largeFileTestDir + "VQSR/expected/g94982_20_1m_10m_tranched_95_99_99.9.vcf"));
        spec.executeTest("testTrancheFilteringTranchesOrder", this);
    }

    @Test
    public void testTrancheFilteringTranchesOrder() throws IOException {
        final String trancheVCF = largeFileTestDir + "VQSR/g94982_20_1m_10m_python_2dcnn.vcf.gz";
        final String indelTruthVCF = largeFileTestDir + "VQSR/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.sites.20.1M-10M.vcf";
        final String snpTruthVCF = largeFileTestDir + "VQSR/Omni25_sites_1525_samples.b37.20.1M-10M.vcf";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, trancheVCF)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s")
                .addArgument("resource", snpTruthVCF)
                .addArgument("resource", indelTruthVCF)
                .addArgument("snp-tranche", "99.9")
                .addArgument("snp-tranche", "95.0")
                .addArgument("snp-tranche", "99.0")
                .addArgument("indel-tranche", "99.9")
                .addArgument("indel-tranche", "95.0")
                .addArgument("indel-tranche", "99.0")
                .addArgument("invalidate-previous-filters", "true")
                .addArgument("info-key", "MIX_SMALL_2D_W_DROPOUT")
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        final IntegrationTestSpec spec = new IntegrationTestSpec(argsBuilder.toString(),
                Arrays.asList(largeFileTestDir + "VQSR/expected/g94982_20_1m_10m_tranched_95_99_99.9.vcf"));
        spec.executeTest("testTrancheFilteringTranchesOrder", this);
    }

    @Test
    public void testTrancheFilteringWithOldFilters() throws IOException {
        final String trancheVCF = largeFileTestDir + "VQSR/g94982_20_1m_10m_python_2dcnn.vcf.gz";
        final String indelTruthVCF = largeFileTestDir + "VQSR/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.sites.20.1M-10M.vcf";
        final String snpTruthVCF = largeFileTestDir + "VQSR/Omni25_sites_1525_samples.b37.20.1M-10M.vcf";

        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, trancheVCF)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s")
                .addArgument("resource", snpTruthVCF)
                .addArgument("resource", indelTruthVCF)
                .addArgument("snp-tranche", "99.0")
                .addArgument("indel-tranche", "99.0")
                .addArgument("info-key", "MIX_SMALL_2D_W_DROPOUT")
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");


        final IntegrationTestSpec spec = new IntegrationTestSpec(argsBuilder.toString(),
                Arrays.asList(largeFileTestDir + "VQSR/expected/g94982_20_1m_10m_tranched_99_with_old_filters.vcf"));
        spec.executeTest("testTrancheFiltering", this);
    }

}
