package org.broadinstitute.hellbender.tools.spark.pipelines;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

public class ReadsPipelineSparkIntegrationTest extends CommandLineProgramTest {

    private static final class PipelineTest {
        final String referenceURL;
        final String bam;
        final String knownSites;
        final String args;
        final String expectedFileName;

        private PipelineTest(String referenceURL, String bam, String knownSites, String args, String expectedFileName) {
            this.referenceURL = referenceURL;
            this.bam = bam;
            this.knownSites = knownSites;
            this.args = args;
            this.expectedFileName = expectedFileName;
        }

        public String getCommandLineNoApiKey() {
            return  " -R " + referenceURL +
                    " -I " + bam +
                    " " + args +
                    (knownSites.isEmpty() ? "": " --knownSites " + knownSites) +
                    " -O %s";
        }

        public String getCommandLine() {
            return  getCommandLineNoApiKey() +
                    " --apiKey " + getGCPTestApiKey();
        }

        @Override
        public String toString() {
            return String.format("ReadsPipeline(bam='%s', args='%s')", bam, args);
        }
    }

    private String getResourceDir(){
        return getTestDataDir() + "/" + "BQSR" + "/";
    }

    @DataProvider(name = "ReadsPipeline")
    public Object[][] createReadsPipelineSparkTestData() {
        final String GRCh37Ref_2021 = b37_reference_20_21;
        final String GRCh37Ref2bit_chr2021 = b37_2bit_reference_20_21;
        final String hiSeqBam_chr20 = getResourceDir() + "CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.noMD.noBQSR.bam";
        final String dbSNPb37_20 = getResourceDir() + DBSNP_138_B37_CH20_1M_1M1K_VCF;

        final String EXPECTED_BAM = "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.noMD.noBQSR.md.bqsr.bam";

        return new Object[][]{
                // input local, computation local.
                //Note: these output files were created by running Picard 1.130 and GATK3.46
                {new PipelineTest(GRCh37Ref2bit_chr2021, hiSeqBam_chr20, dbSNPb37_20, "", getResourceDir() + EXPECTED_BAM)},
        };
    }

    @Test(dataProvider = "ReadsPipeline")
    public void testReadsPipelineSpark(PipelineTest params) throws IOException {
        ArgumentsBuilder ab = new ArgumentsBuilder().add(params.getCommandLineNoApiKey());
        IntegrationTestSpec spec = new IntegrationTestSpec(
                ab.getString(),
                Arrays.asList(params.expectedFileName));
        spec.setCompareBamFilesSorted(true);
        spec.setValidationStringency(ValidationStringency.SILENT);
        spec.executeTest("testReadsPipeline-" + params.args, this);
    }
}
