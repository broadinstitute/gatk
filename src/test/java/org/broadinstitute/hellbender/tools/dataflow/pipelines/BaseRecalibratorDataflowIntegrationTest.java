package org.broadinstitute.hellbender.tools.dataflow.pipelines;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefAPISource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.ApplyBQSR;
import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
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
        final String referenceSetID;
        final String bam;
        final String knownSites;
        final String args;
        final String expectedFileName;

        private BQSRTest(String referenceSetID, String bam, String knownSites, String args, String expectedFileName) {
            this.referenceSetID = referenceSetID;
            this.bam = bam;
            this.knownSites = knownSites;
            this.args = args;
            this.expectedFileName = expectedFileName;
        }

        public String getCommandLine() {
            return  " -R " + RefAPISource.URL_PREFIX + referenceSetID +
                    " -I " + bam +
                    " " + args +
                    (knownSites.isEmpty() ? "": " -knownSites " + knownSites) +
                    " -O %s" +
                    " -sortAllCols" +
                    " --apiKey " + getDataflowTestApiKey();
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
        // we need an API key to get to the reference
        final String apiArgs = " --project " + getDataflowTestProject() + " ";
        final String localResources =  getResourceDir();
        final String hg19Ref = "EMWV_ZfLxrDY-wE";
        final String GRCh37Ref = RefAPISource.GRCH37_REF_ID;
        final String HiSeqBam = localResources + "CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam";
        final String dbSNPb37 =  getResourceDir() + "dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf";
        final String moreSites = getResourceDir() + "bqsr.fakeSitesForTesting.b37.chr17.vcf"; //for testing 2 input files

        return new Object[][]{
            // local computation and files (except for the reference)
            {new BQSRTest(GRCh37Ref, HiSeqBam, dbSNPb37, apiArgs + "", localResources + "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.recal.txt")},
            {new BQSRTest(GRCh37Ref, HiSeqBam, dbSNPb37, apiArgs + "-knownSites " + moreSites, localResources + "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.2inputs.recal.txt")},
            {new BQSRTest(GRCh37Ref, HiSeqBam, dbSNPb37, apiArgs + "--indels_context_size 4",  localResources + "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.indels_context_size_4.recal.txt")},
            {new BQSRTest(GRCh37Ref, HiSeqBam, dbSNPb37, apiArgs + "--low_quality_tail 5",     localResources + "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.low_quality_tail_5.recal.txt")},
            {new BQSRTest(GRCh37Ref, HiSeqBam, dbSNPb37, apiArgs + "--quantizing_levels 6",    localResources + "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.quantizing_levels_6.recal.txt")},
            {new BQSRTest(GRCh37Ref, HiSeqBam, dbSNPb37, apiArgs + "--mismatches_context_size 4", localResources + "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.mismatches_context_size_4.recal.txt")},
            //// //{new BQSRTest(b36Reference, origQualsBam, dbSNPb36, "-OQ", getResourceDir() + "expected.originalQuals.1kg.chr1.1-1K.1RG.dictFix.OQ.txt")},
        };
    }

    @DataProvider(name = "BQSRTestBucket")
    public Object[][] createBQSRTestDataBucket() {
        final String apiArgs = " --project " + getDataflowTestProject();
        final String GRCh37Ref = RefAPISource.GRCH37_REF_ID;
        final String localResources =  getResourceDir();
        final String HiSeqBamCloud = getCloudInputs() + "CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam";
        final String dbSNPb37 =  getResourceDir() + "dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf";

        return new Object[][]{
            // input in cloud, computation local.
            {new BQSRTest(GRCh37Ref, HiSeqBamCloud, dbSNPb37, apiArgs + "", localResources + "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.recal.txt")},
        };
    }

    @DataProvider(name = "BQSRTestCloud")
    public Object[][] createBQSRTestDataCloud() {
        final String cloudArgs = "--runner BLOCKING --project " + getDataflowTestProject() + " --staging " + getDataflowTestStaging();
        final String GRCh37Ref = RefAPISource.GRCH37_REF_ID;
        final String localResources =  getResourceDir();
        final String HiSeqBamCloud = getCloudInputs() + "CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam";
        final String dbSNPb37 =  getResourceDir() + "dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf";

        return new Object[][]{
            // reference and input in cloud, computation in cloud.
            {new BQSRTest(GRCh37Ref, HiSeqBamCloud, dbSNPb37, cloudArgs, localResources + "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.recal.txt")},
        };
    }

    // "local", but we're still getting the reference from the cloud.
    @Test(dataProvider = "BQSRTest", groups = {"cloud"})
    public void testBQSRLocal(BQSRTest params) throws IOException {
        ArgumentsBuilder ab = new ArgumentsBuilder().add(params.getCommandLine());
        addDataflowRunnerArgs(ab);
        IntegrationTestSpec spec = new IntegrationTestSpec(
                ab.getString(),
                Arrays.asList(params.expectedFileName));
        spec.executeTest("testBQSR-" + params.args, this);
    }


    @Test(dataProvider = "BQSRTestBucket", groups = {"bucket"})
    public void testBQSRBucket(BQSRTest params) throws IOException {
        ArgumentsBuilder ab = new ArgumentsBuilder().add(params.getCommandLine());
        addDataflowRunnerArgs(ab);
        IntegrationTestSpec spec = new IntegrationTestSpec(
                ab.getString(),
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

    // TODO: enable this once we figure out how to read bams without requiring them to be indexed.
    @Test(description = "This is to test https://github.com/broadinstitute/hellbender/issues/322", groups = {"cloud"}, enabled = false)
    public void testPlottingWorkflow() throws IOException {
        final String resourceDir = getTestDataDir() + "/" + "BQSR" + "/";
        final String GRCh37Ref = RefAPISource.GRCH37_REF_ID; // that's the "full" version
        final String dbSNPb37 =  getResourceDir() + "dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf";
        final String HiSeqBam = getResourceDir() + "CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam";

        final File actualHiSeqBam_recalibrated = createTempFile("actual.NA12878.chr17_69k_70k.dictFix.recalibrated", ".bam");

        final String tablePre = createTempFile("gatk4.pre.cols", ".table").getAbsolutePath();
        final String argPre = " -apiref " + RefAPISource.URL_PREFIX + GRCh37Ref + " -knownSites " + dbSNPb37 + " -I " + HiSeqBam
                + " -O " + tablePre + " --sort_by_all_columns true";
        new BaseRecalibratorDataflow().instanceMain(Utils.escapeExpressions(argPre));

        final String argApply = "-I " + HiSeqBam + " --bqsr_recal_file " + tablePre+ "  -O " + actualHiSeqBam_recalibrated.getAbsolutePath();
        new ApplyBQSR().instanceMain(Utils.escapeExpressions(argApply));

        final File actualTablePost = createTempFile("gatk4.post.cols", ".table");
        final String argsPost = " -apiref " + RefAPISource.URL_PREFIX + GRCh37Ref + " -knownSites " + dbSNPb37 + " -I " + actualHiSeqBam_recalibrated.getAbsolutePath()
                + " -O " + actualTablePost.getAbsolutePath() + " --sort_by_all_columns true";
        // currently fails with:
        // org.broadinstitute.hellbender.exceptions.UserException: A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
        // Please index all input files:
        // samtools index /tmp/actual.NA12878.chr17_69k_70k.dictFix.recalibrated360293023111112542.bam
        new BaseRecalibratorDataflow().instanceMain(Utils.escapeExpressions(argsPost));

        final File expectedHiSeqBam_recalibrated = new File(resourceDir + "expected.NA12878.chr17_69k_70k.dictFix.recalibrated.bam");

        SamAssertionUtils.assertSamsEqual(actualHiSeqBam_recalibrated, expectedHiSeqBam_recalibrated, ValidationStringency.LENIENT);

        final File expectedTablePost = new File(getResourceDir() + "expected.NA12878.chr17_69k_70k.postRecalibrated.txt");
        // this fails because "Cannot compare coordinate-sorted SAM files because sequence dictionaries differ."
        IntegrationTestSpec.assertEqualTextFiles(actualTablePost, expectedTablePost);
    }

    @Test(groups = {"cloud"})
    public void testBQSRFailWithoutDBSNP() throws IOException {
        final String resourceDir =  getTestDataDir() + "/" + "BQSR" + "/";
        final String apiArgs = " --project " + getDataflowTestProject() + " ";
        final String localResources =  getResourceDir();

        final String GRCh37Ref = RefAPISource.GRCH37_REF_ID; // that's the "full" version
        final String HiSeqBam = resourceDir + "NA12878.chr17_69k_70k.dictFix.bam";

        final String  NO_DBSNP = "";
        final BQSRTest params = new BQSRTest(GRCh37Ref, HiSeqBam, NO_DBSNP , apiArgs, localResources + "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.recal.txt");
        IntegrationTestSpec spec = new IntegrationTestSpec(
                params.getCommandLine(),
                1,
                UserException.CommandLineException.class);
        spec.executeTest("testBQSRFailWithoutDBSNP", this);
    }

    @Test(groups = {"cloud"})
    public void testBQSRFailWithIncompatibleReference() throws IOException {
        final String resourceDir =  getTestDataDir() + "/" + "BQSR" + "/";
        final String apiArgs = " --project " + getDataflowTestProject() + " ";
        final String localResources =  getResourceDir();

        final String hg19Ref = "EMWV_ZfLxrDY-wE";
        final String HiSeqBam = resourceDir + "NA12878.chr17_69k_70k.dictFix.bam";

        final String dbSNPb37 =  getResourceDir() + "dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf";
        final BQSRTest params = new BQSRTest(hg19Ref, HiSeqBam, dbSNPb37 , apiArgs, localResources + "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.recal.txt");
        IntegrationTestSpec spec = new IntegrationTestSpec(
                params.getCommandLine(),
                1,
                UserException.IncompatibleSequenceDictionaries.class);
        spec.executeTest("testBQSRFailWithIncompatibleReference", this);
    }
}
