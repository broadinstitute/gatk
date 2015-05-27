package org.broadinstitute.hellbender.tools.walkers.bqsr;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.dev.tools.walkers.bqsr.BaseRecalibratorDataflow;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.ApplyBQSR;
import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

/**
 *
 * To run the cloud tests, you need:
 * (1) a client-secrets.json in the project's base folder
 * (2) set the TEST_PROJECT and TEST_STAGING_GCS_FOLDER environment variables
 * (3) set the HELLBENDER_TEST_INPUTS environment variable to a GCP folder that has the test inputs this needs
 *     (from src/test/resources/org/broadinstitute/hellbender/tools/BQSR/, see test code
 *      to see which specific files are used)
 *     e.g. gs://MYBUCKET/TESTINPUTS/
 */
public final class BaseRecalibratorDataflowIntegrationTest extends CommandLineProgramTest{

    private final static String CLOUD_INPUTS = System.getenv("HELLBENDER_TEST_INPUTS");
    private final static String TEST_PROJECT = System.getenv("TEST_PROJECT");
    private final static String TEST_STAGING_GCS_FOLDER = System.getenv("TEST_STAGING_GCS_FOLDER");

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
                    " --secret client-secrets.json" +
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

    private String getResourceDir(){
        return getTestDataDir() + "/" + "BQSR" + "/";
    }

    @DataProvider(name = "BQSRTest")
    public Object[][] createBQSRTestData() {
        final String hg18Reference = publicTestDir + "human_g1k_v37.chr17_1Mb.fasta";
        final String b36Reference = getResourceDir() + "human_b36_both.chr1_1k.fasta";
        final String HiSeqBam = getResourceDir() + "NA12878.chr17_69k_70k.dictFix.bam";
        final String dbSNPb37 =  getResourceDir() + "dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf";
        final String origQualsBam = getResourceDir() + "originalQuals.1kg.chr1.1-1K.1RG.dictFix.bam";
        final String dbSNPb36 = getResourceDir() + "dbsnp_132.b36.excluding_sites_after_129.chr1_1k.vcf";

        final String moreSites = getResourceDir() + "bqsr.fakeSitesForTesting.b37.chr17.vcf"; //for testing 2 input files


        return new Object[][]{
                // local files and computation
                {new BQSRTest(hg18Reference, HiSeqBam, dbSNPb37, "", getResourceDir() + "expected.NA12878.chr17_69k_70k.txt")},
                {new BQSRTest(hg18Reference, HiSeqBam, dbSNPb37, "-knownSites " + moreSites, getResourceDir() + "expected.NA12878.chr17_69k_70k.2inputs.txt")},
                {new BQSRTest(hg18Reference, HiSeqBam, dbSNPb37, "--indels_context_size 4", getResourceDir() + "expected.NA12878.chr17_69k_70k.indels_context_size4.txt")},
                {new BQSRTest(hg18Reference, HiSeqBam, dbSNPb37, "--low_quality_tail 5", getResourceDir() + "expected.NA12878.chr17_69k_70k.low_quality_tail5.txt")},
                {new BQSRTest(hg18Reference, HiSeqBam, dbSNPb37, "--quantizing_levels 6", getResourceDir() + "expected.NA12878.chr17_69k_70k.quantizing_levels6.txt")},
                {new BQSRTest(hg18Reference, HiSeqBam, dbSNPb37, "--mismatches_context_size 4", getResourceDir() + "expected.NA12878.chr17_69k_70k.mismatches_context_size4.txt")},
                {new BQSRTest(b36Reference, origQualsBam, dbSNPb36, "-OQ", getResourceDir() + "expected.originalQuals.1kg.chr1.1-1K.1RG.dictFix.OQ.txt")},
        };
    }

    @DataProvider(name = "BQSRTestCloud")
    public Object[][] createBQSRTestDataCloud() {
        final String cloudArgs = "--runner BLOCKING --project " + TEST_PROJECT + " --staging " + TEST_STAGING_GCS_FOLDER;
        final String hg18Reference = publicTestDir + "human_g1k_v37.chr17_1Mb.fasta";
        final String hg18ReferenceCloud = CLOUD_INPUTS + "human_g1k_v37.chr17_1Mb.fasta";
        final String HiSeqBam = getResourceDir() + "NA12878.chr17_69k_70k.dictFix.bam";
        final String HiSeqBamCloud = CLOUD_INPUTS + "NA12878.chr17_69k_70k.dictFix.bam";
        final String dbSNPb37 =  getResourceDir() + "dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf";

        final String moreSites = getResourceDir() + "bqsr.fakeSitesForTesting.b37.chr17.vcf"; //for testing 2 input files


        return new Object[][]{
                // reference in cloud, computation local.
                {new BQSRTest(hg18ReferenceCloud, HiSeqBam, dbSNPb37, "-knownSites " + moreSites, getResourceDir() + "expected.NA12878.chr17_69k_70k.2inputs.txt")},
                // input in cloud, computation local.
                {new BQSRTest(hg18Reference, HiSeqBamCloud, dbSNPb37, "", getResourceDir() + "expected.NA12878.chr17_69k_70k.txt")},
                // reference in cloud, compute in cloud.
                {new BQSRTest(hg18ReferenceCloud, HiSeqBam, dbSNPb37, cloudArgs + " -knownSites " + moreSites, getResourceDir() + "expected.NA12878.chr17_69k_70k.2inputs.txt")},
                // reference and input in cloud, computation in cloud.
                {new BQSRTest(hg18ReferenceCloud, HiSeqBamCloud, dbSNPb37, cloudArgs, getResourceDir() + "expected.NA12878.chr17_69k_70k.txt")},
        };
    }


    @Test(dataProvider = "BQSRTest")
    public void testBQSRLocal(BQSRTest params) throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                params.getCommandLine(),
                Arrays.asList(params.expectedFileName));
        spec.executeTest("testBQSR-" + params.args, this);
    }

    @Test(dataProvider = "BQSRTestCloud", groups = {"cloud"})
    public void testBQSRCloud(BQSRTest params) throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                params.getCommandLine(),
                Arrays.asList(params.expectedFileName));
        spec.executeTest("testBQSR-" + params.args, this);
    }


    @Test(description = "This is to test https://github.com/broadinstitute/hellbender/issues/322", groups = {"cloud"})
    public void testPlottingWorkflow() throws IOException {
        final String cloudArgs = "--secret client-secrets.json ";
        final String resourceDir = getTestDataDir() + "/" + "BQSR" + "/";
        final String hg18Reference = publicTestDir + "human_g1k_v37.chr17_1Mb.fasta";
        final String dbSNPb37 =  getResourceDir() + "dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf";
        final String HiSeqBam = getResourceDir() + "NA12878.chr17_69k_70k.dictFix.bam";

        final File actualHiSeqBam_recalibrated = createTempFile("actual.NA12878.chr17_69k_70k.dictFix.recalibrated", ".bam");

        final String tablePre = createTempFile("gatk4.pre.cols", ".table").getAbsolutePath();
        final String argPre = cloudArgs + "-R " + hg18Reference + " --knownSites " + dbSNPb37 + " -I " + HiSeqBam + " -RECAL_TABLE_FILE " + tablePre + " --sort_by_all_columns true";
        new BaseRecalibratorDataflow().instanceMain(Utils.escapeExpressions(argPre));

        final String argApply = "-I " + HiSeqBam + " --bqsr_recal_file " + tablePre+ "  -O " + actualHiSeqBam_recalibrated.getAbsolutePath();
        new ApplyBQSR().instanceMain(Utils.escapeExpressions(argApply));

        final File actualTablePost = createTempFile("gatk4.post.cols", ".table");
        final String argsPost = cloudArgs
                + " -R " + hg18Reference + " --knownSites " + dbSNPb37 + " -I " + actualHiSeqBam_recalibrated.getAbsolutePath()
                + " -RECAL_TABLE_FILE " + actualTablePost.getAbsolutePath() + " --sort_by_all_columns true";
        new BaseRecalibratorDataflow().instanceMain(Utils.escapeExpressions(argsPost));

        final File expectedHiSeqBam_recalibrated = new File(resourceDir + "expected.NA12878.chr17_69k_70k.dictFix.recalibrated.bam");

        //this fails, disable for now: https://github.com/broadinstitute/hellbender/issues/419
        //IntegrationTestSpec.compareBamFiles(actualHiSeqBam_recalibrated, expectedHiSeqBam_recalibrated);

        final File expectedTablePost = new File(getResourceDir() + "expected.NA12878.chr17_69k_70k.postRecalibrated.txt");
        IntegrationTestSpec.compareTextFiles(actualTablePost, expectedTablePost);
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
