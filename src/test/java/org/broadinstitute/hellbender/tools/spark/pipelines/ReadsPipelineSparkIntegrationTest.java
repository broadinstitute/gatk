package org.broadinstitute.hellbender.tools.spark.pipelines;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MarkDuplicatesSparkArgumentCollection;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceTwoBitSource;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerIntegrationTest;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
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
        final String expectedBamFileName;
        final String expectedVcfFileName;

        private PipelineTest(String referenceURL, String bam, String outputExtension, String knownSites, String args, String expectedBamFileName, String expectedVcfFileName) {
            this.referenceURL = referenceURL;
            this.bam = bam;
            this.outputExtension = outputExtension;
            this.knownSites = knownSites;
            this.args = args;
            this.expectedBamFileName = expectedBamFileName;
            this.expectedVcfFileName = expectedVcfFileName;
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
        final String GRCh37Ref_2021 = GATKBaseTest.b37_reference_20_21;
        final String GRCh37Ref2bit_chr2021 = GATKBaseTest.b37_2bit_reference_20_21;
        final String GRCh37Ref_2021_img = GATKBaseTest.b37_reference_20_21_img;
        final String hiSeqBam_chr20 = getResourceDir() + "CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.noMD.noBQSR.bam";
        final String hiSeqBam_chr20_queryNameSorted = getResourceDir() + "CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.noMD.noBQSR.queryNameSorted.bam";
        final String hiSeqCram_chr20 = getResourceDir() + "CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.noMD.noBQSR.cram";
        final String unalignedBam = largeFileTestDir + "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.tiny.unaligned.bam";
        final String dbSNPb37_20 = getResourceDir() + DBSNP_138_B37_CH20_1M_1M1K_VCF;
        final String more20Sites = getResourceDir() + "dbsnp_138.b37.20.10m-10m100.vcf"; //for testing 2 input files

        final String expectedSingleKnownSites = "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.noMD.noBQSR.md.bqsr.bam";
        final String expectedSingleKnownSitesVcf = "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.noMD.noBQSR.md.bqsr.vcf";
        final String expectedMultipleKnownSites = "expected.MultiSite.reads.pipeline.bam";
        final String expectedMultipleKnownSitesCram = "expected.MultiSite.reads.pipeline.cram";
        final String expectedMultipleKnownSitesVcf = "expected.MultiSite.reads.pipeline.vcf";
        final String expectedMultipleKnownSitesFromUnalignedVcf = "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.tiny.vcf";

        return new Object[][]{
                // input local, computation local.
                // Note: these output files were created by running Picard 1.130 and GATK3.46.
                // VCF files were generated using GATK4 with a command like:
                // gatk HaplotypeCaller \
                //        -I src/test/resources/org/broadinstitute/hellbender/tools/BQSR/expected.MultiSite.reads.pipeline.bam \
                //        -R src/test/resources/large/human_g1k_v37.20.21.fasta \
                //        -O src/test/resources/org/broadinstitute/hellbender/tools/BQSR/expected.MultiSite.reads.pipeline.vcf \
                //        -pairHMM AVX_LOGLESS_CACHING \
                //        -stand_call_conf 30.0
                {new PipelineTest(GRCh37Ref2bit_chr2021, hiSeqBam_chr20, ".bam", dbSNPb37_20, "--join-strategy BROADCAST", null, getResourceDir() + expectedSingleKnownSitesVcf)}, // don't write intermediate BAM
                {new PipelineTest(GRCh37Ref2bit_chr2021, hiSeqBam_chr20, ".bam", dbSNPb37_20, "--join-strategy BROADCAST", getResourceDir() + expectedSingleKnownSites, getResourceDir() + expectedSingleKnownSitesVcf)},
                {new PipelineTest(GRCh37Ref2bit_chr2021, hiSeqBam_chr20, ".bam", dbSNPb37_20, "--join-strategy OVERLAPS_PARTITIONER --read-shard-padding 1000", getResourceDir() + expectedSingleKnownSites, getResourceDir() + expectedSingleKnownSitesVcf)},
                {new PipelineTest(GRCh37Ref2bit_chr2021, hiSeqBam_chr20_queryNameSorted, ".bam", dbSNPb37_20, "--join-strategy BROADCAST", getResourceDir() + expectedMultipleKnownSites, getResourceDir() + expectedSingleKnownSitesVcf)},

                // Output generated with GATK4
                // CRAM test fails since can't use 2bit with CRAM, can only use 2bit with HC.
//                {new PipelineTest(GRCh37Ref2bit_chr2021, hiSeqCram_chr20, ".cram", dbSNPb37_20, "--joinStrategy BROADCAST --knownSites " + more20Sites, getResourceDir() + expectedMultipleKnownSitesCram, getResourceDir() + expectedMultipleKnownSitesVcf)},
                {new PipelineTest(GRCh37Ref2bit_chr2021, hiSeqBam_chr20, ".bam", dbSNPb37_20, "--join-strategy BROADCAST --known-sites " + more20Sites, getResourceDir() + expectedMultipleKnownSites, getResourceDir() + expectedMultipleKnownSitesVcf)},
                {new PipelineTest(GRCh37Ref2bit_chr2021, hiSeqBam_chr20, ".bam", dbSNPb37_20, "--join-strategy OVERLAPS_PARTITIONER --read-shard-padding 1000 --known-sites " + more20Sites, getResourceDir() + expectedMultipleKnownSites, getResourceDir() + expectedMultipleKnownSitesVcf)},

                // BWA-MEM
                {new PipelineTest(GRCh37Ref2bit_chr2021, unalignedBam, ".bam", dbSNPb37_20, "--align --bwa-mem-index-image " + GRCh37Ref_2021_img + " --join-strategy BROADCAST --known-sites " + more20Sites, null, largeFileTestDir + expectedMultipleKnownSitesFromUnalignedVcf)},
        };
    }

    @Test(dataProvider = "ReadsPipeline", groups = "spark")
    public void testReadsPipelineSpark(PipelineTest params) throws IOException {
        File outFile = BaseTest.createTempFile("readSparkPipelineTest", ".vcf");
        File outFileBam = BaseTest.createTempFile("readSparkPipelineTest", params.outputExtension);
        final ArrayList<String> args = new ArrayList<>();

        args.add("-I");
        args.add(new File(params.bam).getAbsolutePath());
        args.add("-O");
        args.add(outFile.getAbsolutePath());
        if (params.expectedBamFileName != null) {
            args.add("--output-bam");
            args.add(outFileBam.getAbsolutePath());
        }

        File referenceFile = null;
        if (params.referenceURL != null) {
            referenceFile = new File(params.referenceURL);
            args.add("-R");
            args.add(referenceFile.getAbsolutePath());
        }
        args.add("--"+ MarkDuplicatesSparkArgumentCollection.DO_NOT_MARK_UNMAPPED_MATES_LONG_NAME);
        args.add("-indels");
        args.add("--enable-baq");
        args.add("-pairHMM");
        args.add("AVX_LOGLESS_CACHING");
        args.add("-stand-call-conf");
        args.add("30.0");
        args.add("--known-sites");
        args.add(params.knownSites);
        if (params.args != null) {
            Stream.of(params.args.trim().split(" ")).forEach(args::add);
        }

        runCommandLine(args);

        if (params.expectedBamFileName != null) {
            if (referenceFile != null && !ReferenceTwoBitSource.isTwoBit(referenceFile.getName())) { // htsjdk can't handle 2bit reference files
                SamAssertionUtils.assertEqualBamFiles(outFileBam, new File(params.expectedBamFileName), referenceFile, true, ValidationStringency.SILENT);
            } else {
                SamAssertionUtils.assertEqualBamFiles(outFileBam, new File(params.expectedBamFileName), true, ValidationStringency.SILENT);
            }
        }

        final double concordance = HaplotypeCallerIntegrationTest.calculateConcordance(outFile, new File(params.expectedVcfFileName));
        Assert.assertTrue(concordance >= 0.99, "Concordance with GATK 4 in VCF mode is < 99% (" +  concordance + ")");
    }
}
