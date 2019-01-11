package org.broadinstitute.hellbender.tools.walkers.bqsr;

import htsjdk.samtools.ValidationStringency;

import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

public final class BaseRecalibratorIntegrationTest extends CommandLineProgramTest{

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

    @DataProvider(name = "BQSRTest")
    public Object[][] createBQSRTestData() {
        final String hg18Reference = publicTestDir + "human_g1k_v37.chr17_1Mb.fasta";
        final String b36Reference = getResourceDir() + "human_b36_both.chr1_1k.fasta";
        final String hiSeqBam_chr20 = getResourceDir() + WGS_B37_CH20_1M_1M1K_BAM;
        final String hiSeqBam_1read = getResourceDir() + "overlappingRead.bam";
        final String hiSeqBam_readNithNoRefBases = getResourceDir() + "NA12878.oq.read_consumes_zero_ref_bases.chr20.bam";
        final String HiSeqBam_chr17 = getResourceDir() + "NA12878.chr17_69k_70k.dictFix.bam";
        final String HiSeqCram_chr17 = getResourceDir() + "NA12878.chr17_69k_70k.dictFix.cram";
        final String trickyBam_chr20 = getResourceDir() + "CEUTrio.HiSeq.WGS.b37.ch20.4379150-4379157.bam";
        final String dbSNPb37_chr17 =  getResourceDir() + "dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf";
        final String dbSNPb37_chr20 = getResourceDir() + DBSNP_138_B37_CH20_1M_1M1K_VCF;
        final String origQualsBam_chr1 = getResourceDir() + "originalQuals.1kg.chr1.1-1K.1RG.dictFix.bam";
        final String dbSNPb36_chr1 = getResourceDir() + "dbsnp_132.b36.excluding_sites_after_129.chr1_1k.vcf";
        final String GRCh37Ref_chr2021 = "src/test/resources/large/human_g1k_v37.20.21.fasta";
        final String more17Sites = getResourceDir() + "bqsr.fakeSitesForTesting.b37.chr17.vcf"; //for testing 2 input files

        final String hiSeqBam_20_21_100000 = getResourceDir() + "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.10m-10m100.bam";
        final String more20Sites = getResourceDir() + "dbsnp_138.b37.20.10m-10m100.vcf"; //for testing 2 input files
        final String more21Sites = getResourceDir() + "dbsnp_138.b37.21.10m-10m100.vcf"; //for testing 2 input files
        return new Object[][]{
                //Note: recal tables were created using GATK3.4 with one change from 2.87 to 2.88 in expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.recal.txt
                // The reason is that GATK4 uses a multiplier in summing doubles in RecalDatum.
                // See MathUtilsUniTest.testAddDoubles for a demonstration how that can change the results.
                // See RecalDatum for explanation of why the multiplier is needed.

                {new BQSRTest(GRCh37Ref_chr2021, hiSeqBam_chr20, dbSNPb37_chr20, "", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_NOINDEL_NOBAQ_RECAL)},
                {new BQSRTest(GRCh37Ref_chr2021, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_RECAL)},
                {new BQSRTest(GRCh37Ref_chr2021, hiSeqBam_1read, dbsnp_138_b37_20_21_vcf, "-indels --enable-baq", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1READ_RECAL)},
                {new BQSRTest(GRCh37Ref_chr2021, hiSeqBam_readNithNoRefBases, dbsnp_138_b37_20_21_vcf, "-indels --enable-baq", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1READ_NOREFBASES_RECAL)},

                {new BQSRTest(GRCh37Ref_chr2021, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " + "--indels-context-size 4",  getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_INDELS_CONTEXT_SIZE_4_RECAL)},
                {new BQSRTest(GRCh37Ref_chr2021, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " + "--low-quality-tail 5",     getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_LOW_QUALITY_TAIL_5_RECAL)},
                {new BQSRTest(GRCh37Ref_chr2021, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +"--quantizing-levels 6",    getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_QUANTIZING_LEVELS_6_RECAL)},
                {new BQSRTest(GRCh37Ref_chr2021, hiSeqBam_chr20, dbSNPb37_chr20, "-indels --enable-baq " +"--mismatches-context-size 4", getResourceDir() + BQSRTestData.EXPECTED_WGS_B37_CH20_1M_1M1K_MISMATCHES_CONTEXT_SIZE_4_RECAL)},

                {new BQSRTest(hg18Reference, HiSeqCram_chr17, dbSNPb37_chr17, "-indels --enable-baq ", getResourceDir() + "expected.NA12878.chr17_69k_70k.txt")},
                {new BQSRTest(hg18Reference, HiSeqBam_chr17, dbSNPb37_chr17, "-indels --enable-baq ", getResourceDir() + "expected.NA12878.chr17_69k_70k.txt")},
                {new BQSRTest(GRCh37Ref_chr2021, trickyBam_chr20, dbSNPb37_chr20, "-indels --enable-baq ", getResourceDir() + "expected.CEUTrio.HiSeq.WGS.b37.ch20.4379150-4379157.recal.txt")},
                {new BQSRTest(hg18Reference, HiSeqBam_chr17, dbSNPb37_chr17, "-indels --enable-baq " +"--known-sites " + more17Sites, getResourceDir() + "expected.NA12878.chr17_69k_70k.2inputs.txt")},
                {new BQSRTest(hg18Reference, HiSeqBam_chr17, dbSNPb37_chr17, "-indels --enable-baq " +"--indels-context-size 4", getResourceDir() + "expected.NA12878.chr17_69k_70k.indels_context_size4.txt")},
                {new BQSRTest(hg18Reference, HiSeqBam_chr17, dbSNPb37_chr17, "-indels --enable-baq " +"--low-quality-tail 5", getResourceDir() + "expected.NA12878.chr17_69k_70k.low_quality_tail5.txt")},
                {new BQSRTest(hg18Reference, HiSeqBam_chr17, dbSNPb37_chr17, "-indels --enable-baq " +"--quantizing-levels 6", getResourceDir() + "expected.NA12878.chr17_69k_70k.quantizing_levels6.txt")},
                {new BQSRTest(hg18Reference, HiSeqBam_chr17, dbSNPb37_chr17, "-indels --enable-baq " +"--mismatches-context-size 4", getResourceDir() + "expected.NA12878.chr17_69k_70k.mismatches_context_size4.txt")},
                {new BQSRTest(b36Reference, origQualsBam_chr1, dbSNPb36_chr1, "-indels --enable-baq " +"-OQ", getResourceDir() + "expected.originalQuals.1kg.chr1.1-1K.1RG.dictFix.OQ.txt")},
        };
    }
    @Test(dataProvider = "BQSRTest")
    public void testBQSR(BQSRTest params) throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                params.getCommandLine(),
                Arrays.asList(params.expectedFileName));
        spec.executeTest("testBQSR-" + params.args, this);
    }

    @Test(description = "This is to test https://github.com/broadinstitute/hellbender/issues/322")
    public void testPlottingWorkflow() throws IOException {
        final String resourceDir = getTestDataDir() + "/" + "BQSR" + "/";
        final String hg18Reference = publicTestDir + "human_g1k_v37.chr17_1Mb.fasta";
        final String dbSNPb37_chr17 =  getResourceDir() + "dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf";
        final String HiSeqBam_chr17 = getResourceDir() + "NA12878.chr17_69k_70k.dictFix.bam";

        final File actualHiSeqBam_recalibrated_chr17 = createTempFile("actual.recalibrated", ".bam");

        final String tablePre = createTempFile("gatk4.pre.cols", ".table").getAbsolutePath();
        final String argPre = "-R " + hg18Reference + " -indels --enable-baq " + " --known-sites " + dbSNPb37_chr17 + " -I " + HiSeqBam_chr17 + " -O " + tablePre + " ";
        new BaseRecalibrator().instanceMain(Utils.escapeExpressions(argPre));

        final String argApply = "-I " + HiSeqBam_chr17 + " --" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME +" " + tablePre + " -O " + actualHiSeqBam_recalibrated_chr17.getAbsolutePath();
        new ApplyBQSR().instanceMain(Utils.escapeExpressions(argApply));

        final File actualTablePost = createTempFile("gatk4.post.cols", ".table");
        final String argsPost = "-R " + hg18Reference + " -indels --enable-baq " + " --known-sites " + dbSNPb37_chr17 + " -I " + actualHiSeqBam_recalibrated_chr17.getAbsolutePath() + " -O " + actualTablePost.getAbsolutePath() + " ";
        new BaseRecalibrator().instanceMain(Utils.escapeExpressions(argsPost));

        final File expectedHiSeqBam_recalibrated_chr17 = new File(resourceDir + "expected.NA12878.chr17_69k_70k.dictFix.recalibrated.DIQ.bam");

        SamAssertionUtils.assertSamsEqual(actualHiSeqBam_recalibrated_chr17, expectedHiSeqBam_recalibrated_chr17, ValidationStringency.LENIENT);

        final File expectedTablePost = new File(getResourceDir() + "expected.NA12878.chr17_69k_70k.postRecalibrated.txt");
        IntegrationTestSpec.assertEqualTextFiles(actualTablePost, expectedTablePost);
    }

    @Test
    public void testBQSRFailWithoutDBSNP() throws IOException {
        final String resourceDir =  getTestDataDir() + "/" + "BQSR" + "/";

        final String hg18Reference = publicTestDir + "human_g1k_v37.chr17_1Mb.fasta";
        final String HiSeqBam_chr17 = resourceDir + "NA12878.chr17_69k_70k.dictFix.bam";

        final String  NO_DBSNP = "";
        final String  NO_ARGS = "";
        final BQSRTest params = new BQSRTest(hg18Reference, HiSeqBam_chr17, NO_DBSNP, NO_ARGS, resourceDir + "expected.NA12878.chr17_69k_70k.txt");
        IntegrationTestSpec spec = new IntegrationTestSpec(
                params.getCommandLine(),
                1,
                CommandLineException.class);
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

    @Test
    public void testBQSRFailWithIncompatibleSequenceDictionaries() throws IOException {
        final String resourceDir =  getTestDataDir() + "/" + "BQSR" + "/";

        final String bam_chr20 = resourceDir + WGS_B37_CH20_1M_1M1K_BAM;
        final String dbSNPb37_chr17 =  getResourceDir() + "dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf";

        final String  NO_ARGS = "";
        final BQSRTest params = new BQSRTest(b37_reference_20_21, bam_chr20, dbSNPb37_chr17, NO_ARGS, resourceDir + "expected.txt");
        IntegrationTestSpec spec = new IntegrationTestSpec(
                params.getCommandLine(),
                1,
                UserException.IncompatibleSequenceDictionaries.class);
        spec.executeTest("testBQSRFailWithIncompatibleSequenceDictionaries", this);
    }

}
