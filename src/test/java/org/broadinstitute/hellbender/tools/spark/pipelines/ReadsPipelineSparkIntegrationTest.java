package org.broadinstitute.hellbender.tools.spark.pipelines;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerIntegrationTest;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.stream.Stream;

public class ReadsPipelineSparkIntegrationTest extends CommandLineProgramTest {

    private static final class PipelineTest {
        final String referenceURL;
        final String bam;
        final String outputExtension;
        final String knownSites;
        final String args;
        final String expectedFileName;

        private PipelineTest(String referenceURL, String bam, String outputExtension, String knownSites, String args, String expectedFileName) {
            this.referenceURL = referenceURL;
            this.bam = bam;
            this.outputExtension = outputExtension;
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
        final String GRCh37Ref2bit_chr2021 = b37_2bit_reference_20_21;
        final String hiSeqBam_chr20 = getResourceDir() + "CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.noMD.noBQSR.bam";
        final String dbSNPb37_20 = getResourceDir() + DBSNP_138_B37_CH20_1M_1M1K_VCF;
        final String more20Sites = getResourceDir() + "dbsnp_138.b37.20.10m-10m100.vcf"; //for testing 2 input files

        final String expectedMultipleKnownSites = "expected.MultiSite.reads.pipeline.vcf";

        return new Object[][]{
                // input local, computation local.
                // Output generated with GATK4 using
                // ./gatk-launch HaplotypeCaller \
                //        -I src/test/resources/org/broadinstitute/hellbender/tools/BQSR/expected.MultiSite.reads.pipeline.bam \
                //        -R src/test/resources/large/human_g1k_v37.20.21.fasta \
                //        -O src/test/resources/org/broadinstitute/hellbender/tools/BQSR/expected.MultiSite.reads.pipeline.vcf \
                //        -pairHMM AVX_LOGLESS_CACHING \
                //        -stand_call_conf 30.0
                {new PipelineTest(GRCh37Ref2bit_chr2021, hiSeqBam_chr20, ".vcf", dbSNPb37_20, "-indelBQSR -enableBAQ -pairHMM AVX_LOGLESS_CACHING -stand_call_conf 30.0 " +"--joinStrategy BROADCAST --knownSites " + more20Sites, getResourceDir() + expectedMultipleKnownSites)},
        };
    }

    @Test(dataProvider = "ReadsPipeline", groups = "spark")
    public void testReadsPipelineSpark(PipelineTest params) throws IOException {
        File outFile = BaseTest.createTempFile("readSparkPipelineTest", params.outputExtension);
        final ArrayList<String> args = new ArrayList<>();

        args.add("-I");
        args.add(new File(params.bam).getAbsolutePath());
        args.add("-O");
        args.add(outFile.getAbsolutePath());

        File referenceFile = null;
        if (params.referenceURL != null) {
            referenceFile = new File(params.referenceURL);
            args.add("-R");
            args.add(referenceFile.getAbsolutePath());
        }
        args.add("-indelBQSR");
        args.add("-enableBAQ");
        args.add("--knownSites");
        args.add(params.knownSites);
        if (params.args != null) {
            Stream.of(params.args.trim().split(" ")).forEach(args::add);
        }

        runCommandLine(args);

        final double concordance = HaplotypeCallerIntegrationTest.calculateConcordance(outFile, new File(params.expectedFileName));
        Assert.assertTrue(concordance >= 0.99, "Concordance with GATK 4 in VCF mode is < 99% (" +  concordance + ")");
    }
}
