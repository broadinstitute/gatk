package org.broadinstitute.hellbender.tools.walkers.bqsr;

import htsjdk.samtools.ValidationStringency;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.dev.tools.walkers.bqsr.BaseRecalibratorDataflow;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.ApplyBQSR;
import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.read.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

/**
 *
 * To run the cloud tests, you need the following environment variables:
 * HELLBENDER_TEST_PROJECT - the short name of your Google Cloud Platform project
 * HELLBENDER_TEST_APIKEY - the API key associated with that project
 * HELLBENDER_TEST_STAGING - a Google Cloud Storage folder to hold temporary files, e.g. gs://MYBUCKET/staging/
 * HELLBENDER_TEST_INPUTS - a Google Cloud Storage path (ending in "/") that contains this folder:
 *                          org/broadinstitute/hellbender/tools/dataflow/BaseRecalibratorDataflow/
 *                          with the BaseRecalibratorDataflow test inputs.
 *                          The files must be shared publicly (i.e. world-readable).
 *                          The specific files used for this test are:
 *                          - human_g1k_v37.chr17_1Mb.*
 *                          - NA12878.chr17_69k_70k.dictFix.*
 */
public final class BaseRecalibratorDataflowIntegrationTest extends CommandLineProgramTest {

    private final static String THIS_TEST_FOLDER = "org/broadinstitute/hellbender/tools/BQSR/";

    // this is a hack so we can run the stuff from a Jar file without too much hassle
    public static void main(String[] args) throws Exception {
        new BaseRecalibratorDataflowIntegrationTest().runAllTests();
    }

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

    private String getResourceDir(){
        return getTestDataDir() + "/" + "BQSR" + "/";
    }

    private String getCloudInputs() {
        return getDataflowTestInputPath() + THIS_TEST_FOLDER;
    }


    public void runAllTests() throws Exception {
        System.out.println("testBQSRLocal");
        for(Object[] in : createBQSRTestData()) {
            BQSRTest testCase = (BQSRTest)in[0];
            testBQSRLocal(testCase);
        }
        System.out.println("testBQSRFailWithoutDBSNP");
        testBQSRFailWithoutDBSNP();
        System.out.println("testBQSRFailWithIncompatibleReference");
        testBQSRFailWithIncompatibleReference();
        System.out.println("testBQSRBucket");
        for(Object[] in : createBQSRTestDataBucket()) {
            BQSRTest testCase = (BQSRTest)in[0];
            testBQSRBucket(testCase);
        }
        System.out.println("testBQSRCloud");
        for(Object[] in : createBQSRTestDataCloud()) {
            BQSRTest testCase = (BQSRTest)in[0];
            testBQSRCloud(testCase);
        }
        System.out.println("Tests passed.");
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

    @DataProvider(name = "BQSRTestBucket")
    public Object[][] createBQSRTestDataBucket() {
        final String hg18Reference = publicTestDir + "human_g1k_v37.chr17_1Mb.fasta";
        final String hg18ReferenceCloud = getCloudInputs() + "human_g1k_v37.chr17_1Mb.fasta";
        final String HiSeqBam = getResourceDir() + "NA12878.chr17_69k_70k.dictFix.bam";
        final String HiSeqBamCloud = getCloudInputs() + "NA12878.chr17_69k_70k.dictFix.bam";
        final String dbSNPb37 =  getResourceDir() + "dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf";

        final String moreSites = getResourceDir() + "bqsr.fakeSitesForTesting.b37.chr17.vcf"; //for testing 2 input files


        return new Object[][]{
                // reference in cloud, computation local.
                {new BQSRTest(hg18ReferenceCloud, HiSeqBam, dbSNPb37, "-knownSites " + moreSites, getResourceDir() + "expected.NA12878.chr17_69k_70k.2inputs.txt")},
                // input in cloud, computation local.
                {new BQSRTest(hg18Reference, HiSeqBamCloud, dbSNPb37, "", getResourceDir() + "expected.NA12878.chr17_69k_70k.txt")},
        };
    }

    @DataProvider(name = "BQSRTestCloud")
    public Object[][] createBQSRTestDataCloud() {
        final String cloudArgs = "--runner BLOCKING --apiKey " + getDataflowTestApiKey() + " --project " + getDataflowTestProject() + " --staging " + getDataflowTestStaging();
        final String hg18ReferenceCloud = getCloudInputs() + "human_g1k_v37.chr17_1Mb.fasta";
        final String HiSeqBam = getResourceDir() + "NA12878.chr17_69k_70k.dictFix.bam";
        final String HiSeqBamCloud = getCloudInputs() + "NA12878.chr17_69k_70k.dictFix.bam";
        final String dbSNPb37 =  getResourceDir() + "dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf";

        final String moreSites = getResourceDir() + "bqsr.fakeSitesForTesting.b37.chr17.vcf"; //for testing 2 input files


        return new Object[][]{
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


    @Test(dataProvider = "BQSRTestBucket", groups = {"bucket"})
    public void testBQSRBucket(BQSRTest params) throws IOException {
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
        final String cloudArgs = "--apiKey " + getDataflowTestApiKey() + " ";
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

        SamAssertionUtils.assertSamsEqual(actualHiSeqBam_recalibrated, expectedHiSeqBam_recalibrated, ValidationStringency.LENIENT);

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
                UserException.IncompatibleSequenceDictionaries.class);
        spec.executeTest("testBQSRFailWithIncompatibleReference", this);
    }
}
