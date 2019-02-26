package org.broadinstitute.hellbender.tools.spark;

import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.walkers.bqsr.AbstractBaseRecalibratorIntegrationTest;
import org.broadinstitute.hellbender.tools.walkers.bqsr.BQSRTestData;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

@Test(groups = "spark")
public final class BaseRecalibratorSparkIntegrationTest extends AbstractBaseRecalibratorIntegrationTest {

    private final static String THIS_TEST_FOLDER = toolsTestDir + "BQSR/";

    private String getCloudInputs() {
        return getGCPTestInputPath() + THIS_TEST_FOLDER;
    }

    //This data provider is for tests that use reference (but not BAM) files stored in buckets
    @DataProvider(name = "BQSRCloudTest")
    public Object[][] createBQSRCloudTestData() {
        final String localResources =  getResourceDir();

        final String GRCh37RefCloud = GCS_b37_CHR20_21_REFERENCE;
        final String hiSeqBam_chr20 = localResources + WGS_B37_CH20_1M_1M1K_BAM;
        final String hiSeqBam_1read = localResources + "overlappingRead.bam";
        final String dbSNPb37_chr20 = localResources + DBSNP_138_B37_CH20_1M_1M1K_VCF;
        final String dbSNPb37_chr2021 = dbsnp_138_b37_20_21_vcf;

        return new Object[][]{
                //Note: recal tables were created using GATK3.4 with one change from 2.87 to 2.88 in expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.recal.txt
                // The reason is that GATK4 uses a multiplier in summing doubles in RecalDatum.
                // See MathUtilsUniTest.testAddDoubles for a demonstration how that can change the results.
                // See RecalDatum for explanation of why the multiplier is needed.

                // local input/computation, using a gs:// reference
                {new BQSRTest(GRCh37RefCloud, hiSeqBam_1read, dbSNPb37_chr2021, "-indels --enable-baq", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1READ_RECAL)},
                {new BQSRTest(GRCh37RefCloud, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq", localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_RECAL)},
                {new BQSRTest(GRCh37RefCloud, hiSeqBam_chr20, dbSNPb37_chr20, "", localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_NOINDEL_NOBAQ_RECAL)},
                {new BQSRTest(GRCh37RefCloud, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq --indels-context-size 4",  localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_INDELS_CONTEXT_SIZE_4_RECAL)},
                {new BQSRTest(GRCh37RefCloud, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq --low-quality-tail 5",     localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_LOW_QUALITY_TAIL_5_RECAL)},
                {new BQSRTest(GRCh37RefCloud, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq --quantizing-levels 6",    localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_QUANTIZING_LEVELS_6_RECAL)},
                {new BQSRTest(GRCh37RefCloud, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq --mismatches-context-size 4", localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_MISMATCHES_CONTEXT_SIZE_4_RECAL)},
        };
    }

    @Test(dataProvider = "BQSRCloudTest", groups = {"bucket", "spark"})
    public void testBQSRSparkCloud(final BQSRTest params) throws IOException {
        ArgumentsBuilder ab = new ArgumentsBuilder().add(params.getCommandLine());
        IntegrationTestSpec spec = new IntegrationTestSpec(
                ab.getString(),
                Arrays.asList(params.expectedFileName));
        spec.executeTest("testBQSRSparkCloud-" + params.args, this);
    }

    //This data provider is for tests that use BAM files stored in buckets
    @DataProvider(name = "BQSRTestBucket")
    public Object[][] createBQSRTestDataBucket() {
        final String GRCh37RefCloud = GCS_b37_CHR20_21_REFERENCE;
        final String chr2021Reference2bit = GCS_b37_CHR20_21_REFERENCE_2BIT;
        final String localResources = getResourceDir();
        final String HiSeqBamCloud_chr20 = getCloudInputs() + WGS_B37_CH20_1M_1M1K_BAM;
        final String dbSNPb37_chr20 = localResources + DBSNP_138_B37_CH20_1M_1M1K_VCF;

        return new Object[][]{
                // input in cloud, computation local.
                {new BQSRTest(GRCh37RefCloud, HiSeqBamCloud_chr20, dbSNPb37_chr20, "-indels --enable-baq "+" --join-strategy SHUFFLE", localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_RECAL)},
                {new BQSRTest(GRCh37RefCloud, HiSeqBamCloud_chr20, dbSNPb37_chr20, " --join-strategy SHUFFLE", localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_NOINDEL_NOBAQ_RECAL)},
                {new BQSRTest(chr2021Reference2bit, HiSeqBamCloud_chr20, dbSNPb37_chr20, "-indels --enable-baq "+" --join-strategy BROADCAST", localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_RECAL)},
                {new BQSRTest(chr2021Reference2bit, HiSeqBamCloud_chr20, dbSNPb37_chr20, " --join-strategy BROADCAST", localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_NOINDEL_NOBAQ_RECAL)},
        };
    }

    // TODO: re-enable once ReadsSparkSource natively supports files in GCS buckets
    @Test(dataProvider = "BQSRTestBucket", groups = {"spark", "bucket"}, enabled = false)
    public void testBQSRBucket(BQSRTest params) throws IOException {
        ArgumentsBuilder ab = new ArgumentsBuilder().add(params.getCommandLine());
        IntegrationTestSpec spec = new IntegrationTestSpec(
                ab.getString(),
                Arrays.asList(params.expectedFileName));
        spec.executeTest("testBQSRBucket-" + params.args, this);
    }
}
