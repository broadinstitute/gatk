package org.broadinstitute.hellbender.tools.spark;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.bqsr.BQSRTestData;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

public class BaseRecalibratorSparkShardedIntegrationTest extends CommandLineProgramTest {

    private final static String THIS_TEST_FOLDER = toolsTestDir + "BQSR/";

    private static class BQSRTest {
        final String referenceURL;
        final String bam;
        final String knownSites;
        final String args;
        final String expectedFileName;

        private BQSRTest(String referenceURL, String bam, String knownSites, String args, String expectedFileName) {
            this.referenceURL = referenceURL;
            this.bam = bam;
            this.knownSites = knownSites;
            this.args = args;
            this.expectedFileName = expectedFileName;
        }

        public String getCommandLine() {
            return  " -R " + referenceURL +
                    " -I " + bam +
                    " " + args +
                    (knownSites.isEmpty() ? "": " --known-sites " + knownSites) +
                    " -O %s";
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
        return getGCPTestInputPath() + THIS_TEST_FOLDER;
    }

    @DataProvider(name = "BQSRTest")
    public Object[][] createBQSRTestData() {
        final String localResources =  getResourceDir();
        final String GRCh37Ref = GCS_b37_CHR20_21_REFERENCE;
        final String HiSeqBam_chr20 = localResources + WGS_B37_CH20_1M_1M1K_BAM;
        final String dbSNPb37_chr20 = localResources + DBSNP_138_B37_CH20_1M_1M1K_VCF;
        final String GRCh37RefLocal = b37_reference_20_21;
        final String more20Sites = getResourceDir() + "bqsr.fakeSitesForTesting.b37.chr20.vcf"; //for testing 2 input files

        return new Object[][]{
                // local computation and files (except for the reference)
                {new BQSRTest(GRCh37Ref, HiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq ", localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_RECAL)},
                {new BQSRTest(GRCh37Ref, HiSeqBam_chr20, dbSNPb37_chr20, " ", localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_NOINDEL_NOBAQ_RECAL)},
                {new BQSRTest(GRCh37Ref, HiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +"--indels-context-size 4",  localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_INDELS_CONTEXT_SIZE_4_RECAL)},
                {new BQSRTest(GRCh37Ref, HiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +"--low-quality-tail 5",     localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_LOW_QUALITY_TAIL_5_RECAL)},
                {new BQSRTest(GRCh37Ref, HiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +"--quantizing-levels 6",    localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_QUANTIZING_LEVELS_6_RECAL)},
                {new BQSRTest(GRCh37Ref, HiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +"--mismatches-context-size 4", localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_MISMATCHES_CONTEXT_SIZE_4_RECAL)},

                //multiple known sites; expected output generated with GATK4 walker BQSR
                {new BQSRTest(GRCh37Ref, HiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +"--known-sites " + more20Sites, localResources + "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.2inputs.recal.txt")},

                //// //{new BQSRTest(b36Reference, origQualsBam, dbSNPb36, "-OQ", getResourceDir() + "expected.originalQuals.1kg.chr1.1-1K.1RG.dictFix.OQ.txt")},

                // local reference
                {new BQSRTest(GRCh37RefLocal, HiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq ", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_RECAL)},
                {new BQSRTest(GRCh37RefLocal, HiSeqBam_chr20, dbSNPb37_chr20, " ", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_NOINDEL_NOBAQ_RECAL)},
        };
    }

    @DataProvider(name = "BQSRTestBucket")
    public Object[][] createBQSRTestDataBucket() {
        final String GRCh37Ref = GCS_b37_CHR20_21_REFERENCE;
        final String localResources =  getResourceDir();
        final String HiSeqBamCloud_chr20 = getCloudInputs() + WGS_B37_CH20_1M_1M1K_BAM;
        final String dbSNPb37_chr20 = localResources + DBSNP_138_B37_CH20_1M_1M1K_VCF;

        return new Object[][]{
                // input in cloud, computation local.
                {new BQSRTest(GRCh37Ref, HiSeqBamCloud_chr20, dbSNPb37_chr20, "-indels --enable-baq ", localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_RECAL)},
                {new BQSRTest(GRCh37Ref, HiSeqBamCloud_chr20, dbSNPb37_chr20, " ", localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_NOINDEL_NOBAQ_RECAL)},
        };
    }


    // "local", but we're still getting the reference from the cloud.
    @Test(dataProvider = "BQSRTest", groups = {"spark", "bucket"}, enabled = false) //FIXME: disabled because it fails. https://github.com/broadinstitute/gatk/issues/1119
    public void testBQSRLocal(BQSRTest params) throws IOException {
        ArgumentsBuilder ab = new ArgumentsBuilder().add(params.getCommandLine());
        IntegrationTestSpec spec = new IntegrationTestSpec(
                ab.getString(),
                Arrays.asList(params.expectedFileName));
        spec.executeTest("testBQSR-" + params.args, this);
    }

    // this one actually passes, but it takes too long (>10min)! Adding -L speeds it up but then the output's wrong (?)
    @Test(dataProvider = "BQSRTestBucket", groups = {"spark", "bucket"}, enabled = false)
    public void testBQSRBucket(BQSRTest params) throws IOException {
        ArgumentsBuilder ab = new ArgumentsBuilder().add(params.getCommandLine());
        IntegrationTestSpec spec = new IntegrationTestSpec(
                ab.getString(),
                Arrays.asList(params.expectedFileName));
        spec.executeTest("testBQSR-" + params.args, this);
    }

    // TODO: We need to update the expected output files for this test, then it can be re-enabled.
    @Test(description = "This is to test https://github.com/broadinstitute/hellbender/issues/322", groups = {"spark", "bucket"}, enabled = false)
    public void testPlottingWorkflow() throws IOException {
        final String resourceDir = getTestDataDir() + "/" + "BQSR" + "/";
        final String dbSNPb37 =  getResourceDir() + "dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf";
        final String HiSeqBam = getResourceDir() + WGS_B37_CH20_1M_1M1K_BAM;

        final File actualHiSeqBam_recalibrated = createTempFile("actual.NA12878.chr17_69k_70k.dictFix.recalibrated", ".bam");

        final String tablePre = createTempFile("gatk4.pre.cols", ".table").getAbsolutePath();
        final String argPre = " -R " +GCS_b37_CHR20_21_REFERENCE + "-indels --enable-baq " +" --known-sites " + dbSNPb37 + " -I " + HiSeqBam
                + " -O " + tablePre;
        new BaseRecalibratorSpark().instanceMain(Utils.escapeExpressions(argPre));

        final String argApply = "-I " + HiSeqBam + " --bqsr-recal-file " + tablePre + " -O " + actualHiSeqBam_recalibrated.getAbsolutePath();
        new ApplyBQSRSpark().instanceMain(Utils.escapeExpressions(argApply));

        final File actualTablePost = createTempFile("gatk4.post.cols", ".table");
        final String argsPost = " -R " + GCS_b37_CHR20_21_REFERENCE + "-indels --enable-baq " +" --known-sites " + dbSNPb37 + " -I " + actualHiSeqBam_recalibrated.getAbsolutePath()
                + " -O " + actualTablePost.getAbsolutePath();
        new BaseRecalibratorSpark().instanceMain(Utils.escapeExpressions(argsPost));

        final File expectedHiSeqBam_recalibrated = new File(resourceDir + "expected.NA12878.chr17_69k_70k.dictFix.recalibrated.DIQ.bam");

        SamAssertionUtils.assertSamsEqual(actualHiSeqBam_recalibrated, expectedHiSeqBam_recalibrated, ValidationStringency.LENIENT);

        final File expectedTablePost = new File(getResourceDir() + "expected.NA12878.chr17_69k_70k.postRecalibrated.txt");
        IntegrationTestSpec.assertEqualTextFiles(actualTablePost, expectedTablePost);
    }

    @Test(groups = {"spark", "bucket"})
    public void testBQSRFailWithoutDBSNP() throws IOException {
        final String resourceDir =  getTestDataDir() + "/" + "BQSR" + "/";
        final String localResources =  getResourceDir();

        final String GRCh37Ref = GCS_b37_CHR20_21_REFERENCE;
        final String HiSeqBam = resourceDir + "NA12878.chr17_69k_70k.dictFix.bam";

        final String  NO_DBSNP = "";
        final BQSRTest params = new BQSRTest(GRCh37Ref, HiSeqBam, NO_DBSNP, "-indels --enable-baq ", localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_RECAL);
        IntegrationTestSpec spec = new IntegrationTestSpec(
                params.getCommandLine(),
                1,
                CommandLineException.class);
        spec.executeTest("testBQSRFailWithoutDBSNP", this);
    }

    @Test(groups = {"spark", "bucket"})
    public void testBQSRFailWithIncompatibleReference() throws IOException {
        final String resourceDir =  getTestDataDir() + "/" + "BQSR" + "/";
        final String localResources =  getResourceDir();

        final String hg19Ref = hg19MicroReference;
        final String HiSeqBam = resourceDir + "NA12878.chr17_69k_70k.dictFix.bam";

        final String dbSNPb37 =  getResourceDir() + "dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf";
        final BQSRTest params = new BQSRTest(hg19Ref, HiSeqBam, dbSNPb37, "-indels --enable-baq ", localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_RECAL);
        IntegrationTestSpec spec = new IntegrationTestSpec(
                params.getCommandLine(),
                1,
                UserException.IncompatibleSequenceDictionaries.class);
        spec.executeTest("testBQSRFailWithIncompatibleReference", this);
    }
}
