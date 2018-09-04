package org.broadinstitute.hellbender.tools.spark;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.bqsr.BQSRTestData;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

public final class BaseRecalibratorSparkIntegrationTest extends CommandLineProgramTest {

    private final static String THIS_TEST_FOLDER = toolsTestDir + "BQSR/";

    private static final class BQSRTest {
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

        final String GRCh37Ref2bit_chr2021 = b37_2bit_reference_20_21;
        final String GRCh37Ref_chr2021 = b37_reference_20_21;
        final String hiSeqBam_chr20 = localResources + WGS_B37_CH20_1M_1M1K_BAM;
        final String hiSeqBam_1read = localResources + "overlappingRead.bam";
        final String dbSNPb37_chr20 = localResources + DBSNP_138_B37_CH20_1M_1M1K_VCF;
        final String dbSNPb37_chr2021 = dbsnp_138_b37_20_21_vcf;

        final String hg19Chr171Mb = publicTestDir + "human_g1k_v37.chr17_1Mb.fasta";
        final String hg19Chr171Mb_2bit = publicTestDir + "human_g1k_v37.chr17_1Mb.2bit";
        final String HiSeqBam_chr17 = localResources + "NA12878.chr17_69k_70k.dictFix.bam";
        final String dbSNPb37_chr17 =  localResources + "dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf";
        final String more17Sites = localResources + "bqsr.fakeSitesForTesting.b37.chr17.vcf";

        final String hiSeqBam_20_21_100000 = localResources + "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.10m-10m100.bam";
        final String hiSeqCram_20_21_100000 = localResources + "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.10m-10m100.cram";
        final String more20Sites = localResources + "dbsnp_138.b37.20.10m-10m100.vcf"; //for testing 2 input files
        final String more21Sites = localResources + "dbsnp_138.b37.21.10m-10m100.vcf"; //for testing 2 input files

        return new Object[][]{
                //Note: recal tables were created using GATK3.4 with one change from 2.87 to 2.88 in expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.recal.txt
                // The reason is that GATK4 uses a multiplier in summing doubles in RecalDatum.
                // See MathUtilsUniTest.testAddDoubles for a demonstration how that can change the results.
                // See RecalDatum for explanation of why the multiplier is needed.

                // local input/computation/reference, SHUFFLE
                {new BQSRTest(GRCh37Ref_chr2021, hiSeqBam_1read, dbSNPb37_chr2021, "-indels --enable-baq " + "--join-strategy SHUFFLE", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1READ_RECAL)},
                {new BQSRTest(GRCh37Ref_chr2021, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +"--join-strategy SHUFFLE", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_RECAL)},
                {new BQSRTest(GRCh37Ref_chr2021, hiSeqBam_chr20, dbSNPb37_chr20, "--join-strategy SHUFFLE", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_NOINDEL_NOBAQ_RECAL)},
                {new BQSRTest(GRCh37Ref_chr2021, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +"--join-strategy SHUFFLE --indels-context-size 4", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_INDELS_CONTEXT_SIZE_4_RECAL)},
                {new BQSRTest(GRCh37Ref_chr2021, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +"--join-strategy SHUFFLE --low-quality-tail 5", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_LOW_QUALITY_TAIL_5_RECAL)},
                {new BQSRTest(GRCh37Ref_chr2021, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +"--join-strategy SHUFFLE --quantizing-levels 6", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_QUANTIZING_LEVELS_6_RECAL)},
                {new BQSRTest(GRCh37Ref_chr2021, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +"--join-strategy SHUFFLE --mismatches-context-size 4", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_MISMATCHES_CONTEXT_SIZE_4_RECAL)},
                // multiple known sites; output generated with GATK4 walker
                {new BQSRTest(GRCh37Ref_chr2021 , hiSeqBam_20_21_100000, more20Sites, "-indels --enable-baq " +" --join-strategy SHUFFLE --known-sites " + more21Sites, getResourceDir() + "expected.CEUTrio.HiSeq.WGS.b37.ch20.ch21.10m-10m100.recal.txt")},
                {new BQSRTest(GRCh37Ref_chr2021 , hiSeqCram_20_21_100000, more20Sites, "-indels --enable-baq " +" --join-strategy SHUFFLE --known-sites " + more21Sites, getResourceDir() + "expected.CEUTrio.HiSeq.WGS.b37.ch20.ch21.10m-10m100.recal.txt")},
                // multiple known sites  with SHUFFLE; entire test case shared with walker version
                {new BQSRTest(hg19Chr171Mb, HiSeqBam_chr17, dbSNPb37_chr17, "-indels --enable-baq " +" --join-strategy SHUFFLE --known-sites " + more17Sites, getResourceDir() + "expected.NA12878.chr17_69k_70k.2inputs.txt")},

                // multiple known sites  with BROADCAST; entire test case shared with walker version
                {new BQSRTest(hg19Chr171Mb_2bit, HiSeqBam_chr17, dbSNPb37_chr17, "-indels --enable-baq " +" --join-strategy BROADCAST --known-sites " + more17Sites, getResourceDir() + "expected.NA12878.chr17_69k_70k.2inputs.txt")},

                // local input/computation, 2Bit Reference, BROADCAST
                {new BQSRTest(GRCh37Ref2bit_chr2021, hiSeqBam_1read, dbSNPb37_chr2021, "-indels --enable-baq " +"--join-strategy BROADCAST", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1READ_RECAL)},
                {new BQSRTest(GRCh37Ref2bit_chr2021, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +"--join-strategy BROADCAST", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_RECAL)},
                {new BQSRTest(GRCh37Ref2bit_chr2021, hiSeqBam_chr20, dbSNPb37_chr20, "--join-strategy BROADCAST", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_NOINDEL_NOBAQ_RECAL)},
                {new BQSRTest(GRCh37Ref2bit_chr2021, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +"--join-strategy BROADCAST --indels-context-size 4", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_INDELS_CONTEXT_SIZE_4_RECAL)},
                {new BQSRTest(GRCh37Ref2bit_chr2021, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +"--join-strategy BROADCAST --low-quality-tail 5", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_LOW_QUALITY_TAIL_5_RECAL)},
                {new BQSRTest(GRCh37Ref2bit_chr2021, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +"--join-strategy BROADCAST --quantizing-levels 6", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_QUANTIZING_LEVELS_6_RECAL)},
                {new BQSRTest(GRCh37Ref2bit_chr2021, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +"--join-strategy BROADCAST --mismatches-context-size 4", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_MISMATCHES_CONTEXT_SIZE_4_RECAL)},
                // multiple known sites with 2bit BROADCAST; same output used for multiple known sites SHUFFLE test above
                {new BQSRTest(GRCh37Ref2bit_chr2021, hiSeqBam_20_21_100000, more20Sites, "-indels --enable-baq " +" --known-sites " + more21Sites, getResourceDir() + "expected.CEUTrio.HiSeq.WGS.b37.ch20.ch21.10m-10m100.recal.txt")},
                // Can't use 2 bit reference with a CRAM file: https://github.com/broadinstitute/gatk/issues/1443
                //{new BQSRTest(GRCh37Ref2bit_chr2021, hiSeqCram_20_21_100000, more20Sites, " -known-sites " + more21Sites, getResourceDir() + "expected.CEUTrio.HiSeq.WGS.b37.ch20.ch21.1m-1m100.recal.txt")},

                //// //{new BQSRTest(b36Reference, origQualsBam, dbSNPb36, "-OQ", getResourceDir() + "expected.originalQuals.1kg.chr1.1-1K.1RG.dictFix.OQ.txt")},

                // multiple known sites  with OVERLAPS_PARTITIONER; entire test case shared with walker version
                {new BQSRTest(hg19Chr171Mb_2bit, HiSeqBam_chr17, dbSNPb37_chr17, "-indels --enable-baq " +" --join-strategy OVERLAPS_PARTITIONER --known-sites " + more17Sites, getResourceDir() + "expected.NA12878.chr17_69k_70k.2inputs.txt")},

                // local input/computation, 2Bit Reference, OVERLAPS_PARTITIONER
                {new BQSRTest(GRCh37Ref2bit_chr2021, hiSeqBam_1read, dbSNPb37_chr2021, "-indels --enable-baq " +"--join-strategy OVERLAPS_PARTITIONER", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1READ_RECAL)},
                {new BQSRTest(GRCh37Ref2bit_chr2021, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +"--join-strategy OVERLAPS_PARTITIONER", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_RECAL)},
                {new BQSRTest(GRCh37Ref2bit_chr2021, hiSeqBam_chr20, dbSNPb37_chr20, "--join-strategy OVERLAPS_PARTITIONER", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_NOINDEL_NOBAQ_RECAL)},
                {new BQSRTest(GRCh37Ref2bit_chr2021, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +"--join-strategy OVERLAPS_PARTITIONER --indels-context-size 4", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_INDELS_CONTEXT_SIZE_4_RECAL)},
                {new BQSRTest(GRCh37Ref2bit_chr2021, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +"--join-strategy OVERLAPS_PARTITIONER --low-quality-tail 5", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_LOW_QUALITY_TAIL_5_RECAL)},
                {new BQSRTest(GRCh37Ref2bit_chr2021, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +"--join-strategy OVERLAPS_PARTITIONER --quantizing-levels 6", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_QUANTIZING_LEVELS_6_RECAL)},
                {new BQSRTest(GRCh37Ref2bit_chr2021, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +"--join-strategy OVERLAPS_PARTITIONER --mismatches-context-size 4", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_MISMATCHES_CONTEXT_SIZE_4_RECAL)},
                // multiple known sites with 2bit OVERLAPS_PARTITIONER; same output used for multiple known sites SHUFFLE test above
                {new BQSRTest(GRCh37Ref2bit_chr2021, hiSeqBam_20_21_100000, more20Sites, "-indels --enable-baq " +" --join-strategy OVERLAPS_PARTITIONER --known-sites " + more21Sites, getResourceDir() + "expected.CEUTrio.HiSeq.WGS.b37.ch20.ch21.10m-10m100.recal.txt")},
        };
    }

    @Test(dataProvider = "BQSRTest", groups = "spark")
    public void testBQSRSpark(BQSRTest params) throws IOException {
        ArgumentsBuilder ab = new ArgumentsBuilder().add(params.getCommandLine());
        IntegrationTestSpec spec = new IntegrationTestSpec(
                ab.getString(),
                Arrays.asList(params.expectedFileName));
        spec.executeTest("testBQSRSpark-" + params.args, this);
    }

    //This data provider is for tests that use reference (but not BAM) files stored in buckets
    @DataProvider(name = "BQSRCloudTest")
    public Object[][] createBQSRCloudTestData() {
        final String localResources =  getResourceDir();

        final String GRCh37RefCloud = GCS_b37_CHR20_21_REFERENCE;
        final String chr2021Reference2bit = GCS_b37_CHR20_21_REFERENCE_2BIT;
        final String hiSeqBam_chr20 = localResources + WGS_B37_CH20_1M_1M1K_BAM;
        final String hiSeqBam_1read = localResources + "overlappingRead.bam";
        final String dbSNPb37_chr20 = localResources + DBSNP_138_B37_CH20_1M_1M1K_VCF;
        final String dbSNPb37_chr2021 = dbsnp_138_b37_20_21_vcf;

        return new Object[][]{
                //Note: recal tables were created using GATK3.4 with one change from 2.87 to 2.88 in expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.recal.txt
                // The reason is that GATK4 uses a multiplier in summing doubles in RecalDatum.
                // See MathUtilsUniTest.testAddDoubles for a demonstration how that can change the results.
                // See RecalDatum for explanation of why the multiplier is needed.

                // local input/computation, using a gs:// reference , SHUFFLE
                {new BQSRTest(GRCh37RefCloud, hiSeqBam_1read, dbSNPb37_chr2021, "-indels --enable-baq " +" --join-strategy SHUFFLE", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1READ_RECAL)},
                {new BQSRTest(GRCh37RefCloud, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +" --join-strategy SHUFFLE", localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_RECAL)},
                {new BQSRTest(GRCh37RefCloud, hiSeqBam_chr20, dbSNPb37_chr20, " --join-strategy SHUFFLE", localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_NOINDEL_NOBAQ_RECAL)},
                {new BQSRTest(GRCh37RefCloud, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +" --join-strategy SHUFFLE --indels-context-size 4",  localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_INDELS_CONTEXT_SIZE_4_RECAL)},
                {new BQSRTest(GRCh37RefCloud, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +" --join-strategy SHUFFLE --low-quality-tail 5",     localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_LOW_QUALITY_TAIL_5_RECAL)},
                {new BQSRTest(GRCh37RefCloud, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +" --join-strategy SHUFFLE --quantizing-levels 6",    localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_QUANTIZING_LEVELS_6_RECAL)},
                {new BQSRTest(GRCh37RefCloud, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +" --join-strategy SHUFFLE --mismatches-context-size 4", localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_MISMATCHES_CONTEXT_SIZE_4_RECAL)},

                // local input/computation, using a 2bit reference file in a GCS bucket, BROADCAST
                {new BQSRTest(chr2021Reference2bit, hiSeqBam_1read, dbSNPb37_chr2021, "-indels --enable-baq " +" --join-strategy BROADCAST", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1READ_RECAL)},
                {new BQSRTest(chr2021Reference2bit, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +" --join-strategy BROADCAST", localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_RECAL)},
                {new BQSRTest(chr2021Reference2bit, hiSeqBam_chr20, dbSNPb37_chr20, " --join-strategy BROADCAST", localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_NOINDEL_NOBAQ_RECAL)},
                {new BQSRTest(chr2021Reference2bit, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +" --join-strategy BROADCAST --indels-context-size 4",  localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_INDELS_CONTEXT_SIZE_4_RECAL)},
                {new BQSRTest(chr2021Reference2bit, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +" --join-strategy BROADCAST --low-quality-tail 5",     localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_LOW_QUALITY_TAIL_5_RECAL)},
                {new BQSRTest(chr2021Reference2bit, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +" --join-strategy BROADCAST --quantizing-levels 6",    localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_QUANTIZING_LEVELS_6_RECAL)},
                {new BQSRTest(chr2021Reference2bit, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +" --join-strategy BROADCAST --mismatches-context-size 4", localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_MISMATCHES_CONTEXT_SIZE_4_RECAL)},
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

    @Test(groups = "spark")
    public void testBlowUpOnBroadcastIncompatibleReference() throws IOException {
        //this should blow up because broadcast requires a 2bit reference
        final String hiSeqBam_chr20 = getResourceDir() + WGS_B37_CH20_1M_1M1K_BAM;
        final String dbSNPb37_chr20 = getResourceDir() + DBSNP_138_B37_CH20_1M_1M1K_VCF;

        BQSRTest params = new BQSRTest(b37_reference_20_21, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +"--join-strategy BROADCAST", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_RECAL);

        ArgumentsBuilder ab = new ArgumentsBuilder().add(params.getCommandLine());
        IntegrationTestSpec spec = new IntegrationTestSpec(
                ab.getString(),
                1,
                UserException.Require2BitReferenceForBroadcast.class);
        spec.executeTest("testBQSR-" + params.args, this);
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

    // TODO: This test is disabled because a new expected result needs to be created.
    @Test(description = "This is to test https://github.com/broadinstitute/hellbender/issues/322", groups = {"cloud", "spark"}, enabled = false)
    public void testPlottingWorkflow() throws IOException {
        final String resourceDir = getTestDataDir() + "/" + "BQSR" + "/";
        final String chr2021Reference2bit = GCS_b37_CHR20_21_REFERENCE_2BIT;
        final String dbSNPb37_chr2021 = resourceDir + DBSNP_138_B37_CH20_1M_1M1K_VCF;
        final String HiSeqBam_chr20 = getResourceDir() + WGS_B37_CH20_1M_1M1K_BAM;

        final File actualHiSeqBam_recalibrated = createTempFile("actual.recalibrated", ".bam");

        final String tablePre = createTempFile("gatk4.pre.cols", ".table").getAbsolutePath();
        final String argPre = " -R " + chr2021Reference2bit + "-indels --enable-baq " +" --known-sites " + dbSNPb37_chr2021 + " -I " + HiSeqBam_chr20
                + " -O " + tablePre;
        new BaseRecalibratorSpark().instanceMain(Utils.escapeExpressions(argPre));

        final String argApply = "-I " + HiSeqBam_chr20 + " --bqsr-recal-file " + tablePre + " -O " + actualHiSeqBam_recalibrated.getAbsolutePath();
        new ApplyBQSRSpark().instanceMain(Utils.escapeExpressions(argApply));

        final File actualTablePost = createTempFile("gatk4.post.cols", ".table");
        final String argsPost = " -R " + chr2021Reference2bit + "-indels --enable-baq " +" --known-sites " + dbSNPb37_chr2021 + " -I " + actualHiSeqBam_recalibrated.getAbsolutePath()
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

        final String chr2021Reference2bit = GCS_b37_CHR20_21_REFERENCE_2BIT;
        final String HiSeqBam_chr17 = resourceDir + "NA12878.chr17_69k_70k.dictFix.bam";

        final String  NO_DBSNP = "";
        final BQSRTest params = new BQSRTest(chr2021Reference2bit, HiSeqBam_chr17, NO_DBSNP, "", localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_RECAL);
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

        final String HiSeqBam_chr17 = resourceDir + "NA12878.chr17_69k_70k.dictFix.bam";

        final String dbSNPb37_chr2021 = resourceDir + DBSNP_138_B37_CH20_1M_1M1K_VCF;
        final BQSRTest params = new BQSRTest(GATKBaseTest.hg19MiniReference, HiSeqBam_chr17, dbSNPb37_chr2021, "", localResources + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_RECAL);
        IntegrationTestSpec spec = new IntegrationTestSpec(
                params.getCommandLine(),
                1,
                UserException.IncompatibleSequenceDictionaries.class);
        spec.executeTest("testBQSRFailWithIncompatibleReference", this);
    }
}
