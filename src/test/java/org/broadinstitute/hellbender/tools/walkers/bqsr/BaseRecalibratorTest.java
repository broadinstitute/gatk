package org.broadinstitute.hellbender.tools.walkers.bqsr;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.WalkerTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

public class BaseRecalibratorTest extends CommandLineProgramTest{

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
                    " -knownSites " + knownSites +
                    " --allow_potentially_misencoded_quality_scores" +  // TODO -- remove me when we get new SOLiD bams
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

        return new Object[][]{
                {new BQSRTest(hg18Reference, HiSeqBam, dbSNPb37, "", resourceDir + "expected.NA12878.chr17_69k_70k.txt")},
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
        WalkerTestSpec spec = new WalkerTestSpec(
                params.getCommandLine(),
                Arrays.asList(params.expectedFileName));
        spec.executeTest("testBQSR-" + params.args, this);
    }

    @Test(enabled = false)
    public void testBQSRFailWithoutDBSNP() throws IOException {
        final String resourceDir = getTestDataDir() + getCommandLineProgramName();

        final String hg18Reference = publicTestDir + "human_g1k_v37.chr17_1Mb.fasta";
        final String HiSeqBam = resourceDir + "NA12878.chr17_69k_70k.dictFix.bam";

        final BQSRTest params = new BQSRTest(hg18Reference, HiSeqBam, "", "", resourceDir + "expected.NA12878.chr17_69k_70k.txt");
        WalkerTestSpec spec = new WalkerTestSpec(
                params.getCommandLine(),
                1,
                UserException.CommandLineException.class);
        spec.executeTest("testBQSRFailWithoutDBSNP", this);
    }
}
