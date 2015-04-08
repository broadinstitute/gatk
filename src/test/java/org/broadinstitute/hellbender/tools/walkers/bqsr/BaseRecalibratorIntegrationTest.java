package org.broadinstitute.hellbender.tools.walkers.bqsr;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

public class BaseRecalibratorIntegrationTest extends CommandLineProgramTest{

    private static class BQSRTest {
        final String reference;
        final String bam;
        final String knownSites;
        final String args;
        final String expectedFileName;

        private BQSRTest(String reference, String bam, String knownSites, String args, String expectedFileName) {
            this.reference = reference;
            this.bam = bam;
            this.knownSites = knownSites;
            this.args = args;
            this.expectedFileName = expectedFileName;
        }

        public String getCommandLine() {
            return  " -R " + reference +
                    " -I " + bam +
                    " " + args +
                    (knownSites.isEmpty() ? "": " -knownSites " + knownSites) +
                    " --RECAL_TABLE_FILE %s" +
                    " -sortAllCols";
        }

        @Override
        public String toString() {
            return String.format("BQSR(bam='%s', args='%s')", bam, args);
        }
    }

    @DataProvider(name = "BQSRTest")
    public Object[][] createBQSRTestData() {
        final String resourceDir = getTestDataDir() + "/" + "BQSR" + "/";

        final String hg18Reference = publicTestDir + "human_g1k_v37.chr17_1Mb.fasta";
        final String b36Reference = resourceDir + "human_b36_both.chr1_1k.fasta";
        final String HiSeqBam = resourceDir + "NA12878.chr17_69k_70k.dictFix.bam";
        final String dbSNPb37 =  resourceDir + "dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf";
        final String origQualsBam = resourceDir + "originalQuals.1kg.chr1.1-1K.1RG.dictFix.bam";
        final String dbSNPb36 = resourceDir + "dbsnp_132.b36.excluding_sites_after_129.chr1_1k.vcf";

        final String moreSites = resourceDir + "bqsr.fakeSitesForTesting.b37.chr17.vcf"; //for testing 2 input files

        return new Object[][]{
                {new BQSRTest(hg18Reference, HiSeqBam, dbSNPb37, "", resourceDir + "expected.NA12878.chr17_69k_70k.txt")},
                {new BQSRTest(hg18Reference, HiSeqBam, dbSNPb37, "-knownSites " + moreSites, resourceDir + "expected.NA12878.chr17_69k_70k.2inputs.txt")},
                {new BQSRTest(hg18Reference, HiSeqBam, dbSNPb37, "--no_standard_covs -cov ContextCovariate", resourceDir + "expected.NA12878.chr17_69k_70k.ContextCovariate.txt")},
                {new BQSRTest(hg18Reference, HiSeqBam, dbSNPb37, "--no_standard_covs -cov CycleCovariate", resourceDir + "expected.NA12878.chr17_69k_70k.CycleCovariate.txt")},
                {new BQSRTest(hg18Reference, HiSeqBam, dbSNPb37, "--indels_context_size 4", resourceDir + "expected.NA12878.chr17_69k_70k.indels_context_size4.txt")},
                {new BQSRTest(hg18Reference, HiSeqBam, dbSNPb37, "--low_quality_tail 5", resourceDir + "expected.NA12878.chr17_69k_70k.low_quality_tail5.txt")},
                {new BQSRTest(hg18Reference, HiSeqBam, dbSNPb37, "--quantizing_levels 6", resourceDir + "expected.NA12878.chr17_69k_70k.quantizing_levels6.txt")},
                {new BQSRTest(hg18Reference, HiSeqBam, dbSNPb37, "--mismatches_context_size 4", resourceDir + "expected.NA12878.chr17_69k_70k.mismatches_context_size4.txt")},
                {new BQSRTest(b36Reference, origQualsBam, dbSNPb36, "-OQ", resourceDir + "expected.originalQuals.1kg.chr1.1-1K.1RG.dictFix.OQ.txt")},
        };
    }
    @Test(dataProvider = "BQSRTest")
    public void testBQSR(BQSRTest params) throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                params.getCommandLine(),
                Arrays.asList(params.expectedFileName));
        spec.executeTest("testBQSR-" + params.args, this);
    }


    @Test
    public void testBQSRFailWithoutDBSNP() throws IOException {
        final String resourceDir =  getTestDataDir() + "/" + "BQSR" + "/";

        final String hg18Reference = publicTestDir + "human_g1k_v37.chr17_1Mb.fasta";
        final String HiSeqBam = resourceDir + "NA12878.chr17_69k_70k.dictFix.bam";

        final String  NO_DBSNP = "";
        final String  NO_ARGS = "";
        final BQSRTest params = new BQSRTest(hg18Reference, HiSeqBam, NO_DBSNP, NO_ARGS, resourceDir + "expected.NA12878.chr17_69k_70k.txt");
        IntegrationTestSpec spec = new IntegrationTestSpec(
                params.getCommandLine(),
                1,
                UserException.CommandLineException.class);
        spec.executeTest("testBQSRFailWithoutDBSNP", this);
    }

    @Test
    public void testBQSRFailWithIncompatibleReference() throws IOException {
        final String resourceDir =  getTestDataDir() + "/" + "BQSR" + "/";

        final String HiSeqBam_Hg18 = resourceDir + "HiSeq.1mb.1RG.2k_lines.bam";

        final String  NO_ARGS = "";
        final BQSRTest params = new BQSRTest(hg19MiniReference, HiSeqBam_Hg18, hg19_chr1_1M_dbSNP, NO_ARGS, resourceDir + "expected.txt");
        IntegrationTestSpec spec = new IntegrationTestSpec(
                params.getCommandLine(),
                1,
                UserException.MissingContigInSequenceDictionary.class);
        spec.executeTest("testBQSRFailWithIncompatibleReference", this);
    }
}
