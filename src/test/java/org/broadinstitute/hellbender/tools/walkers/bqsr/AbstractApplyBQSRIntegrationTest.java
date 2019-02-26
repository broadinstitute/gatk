package org.broadinstitute.hellbender.tools.walkers.bqsr;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Stream;

public abstract class AbstractApplyBQSRIntegrationTest extends CommandLineProgramTest {
    public static class ABQSRTest {
        public final String bam;
        public final String reference;
        public final String outputExtension;
        public final String args[];
        public final String expectedFile;

        public ABQSRTest(String bam, String reference, String outputExtension, String args[], String expectedFile) {
            this.bam= bam;
            this.reference = reference;
            this.outputExtension = outputExtension;
            this.args = args;
            this.expectedFile = expectedFile;
        }

        @Override
        public String toString() {
            return String.format("ApplyBQSR(args='%s')", args == null ? "" : StringUtils.join(args));
        }
    }

    protected final String resourceDir = getTestDataDir() + "/" + "BQSR" + "/";
    protected final String hg18Reference = publicTestDir + "human_g1k_v37.chr17_1Mb.fasta";
    protected final String hiSeqBam = resourceDir + "HiSeq.1mb.1RG.2k_lines.alternate.bam";
    protected final String hiSeqCram = resourceDir + "HiSeq.1mb.1RG.2k_lines.alternate.cram";
    protected final String hiSeqBamAligned = resourceDir + "HiSeq.1mb.1RG.2k_lines.alternate_allaligned.bam";
    protected final String hiSeqCramAligned = resourceDir + "HiSeq.1mb.1RG.2k_lines.alternate_allaligned.cram";

    protected String getDevNull() throws IOException {
        return "/dev/null";
    }

    @DataProvider(name = "ApplyBQSRTest")
    public Object[][] createABQSRTestData() {
        List<Object[]> tests = new ArrayList<>();

        //Note: these outputs were created using GATK3
        tests.add(new Object[]{new ABQSRTest(hiSeqBam, null, ".bam", null, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBam, null, ".bam", new String[] {"-OQ"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.OQ.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBam, null, ".bam", new String[] {"--quantize-quals", "-1"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.qq-1.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBam, null, ".bam", new String[] {"--quantize-quals", "6"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.qq6.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBam, null, ".bam", new String[] {"--static-quantized-quals", "10", "--static-quantized-quals", "20", "--static-quantized-quals", "30"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.SQQ102030.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBam, null, ".bam", new String[] {"--static-quantized-quals", "10", "--static-quantized-quals", "20", "--static-quantized-quals", "30", "--round-down-quantized"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.SQQ102030RDQ.bam")});

        tests.add(new Object[]{new ABQSRTest(hiSeqBamAligned, null, ".bam", null, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBamAligned, null, ".bam", new String[] {"-OQ"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.OQ.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBamAligned, null, ".bam", new String[] {"--quantize-quals", "-1"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.qq-1.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBamAligned, null, ".bam", new String[] {"--quantize-quals", "6"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.qq6.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBamAligned, null, ".bam", new String[] {"--static-quantized-quals", "10", "--static-quantized-quals", "20", "--static-quantized-quals", "30"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.SQQ102030.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBamAligned, null, ".bam", new String[] {"--static-quantized-quals", "10", "--static-quantized-quals", "20", "--static-quantized-quals", "30", "--round-down-quantized"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.SQQ102030RDQ.bam")});

        //CRAM - input and output crams generated by direct conversion of the corresponding BAM test files with samtools 1.3
        tests.add(new Object[]{new ABQSRTest(hiSeqCram, hg18Reference, ".cram", new String[] {"--" + StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME, "true"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.cram")});
        tests.add(new Object[]{new ABQSRTest(hiSeqCramAligned, hg18Reference, ".cram", new String[] {"--quantize-quals", "6", "--" + StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME, "true"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.qq6.cram")});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ApplyBQSRTest")
    public void testApplyBQSRFile(ABQSRTest params) throws IOException {
        File outFile = GATKBaseTest.createTempFile("applyBQSRTest", params.outputExtension);
        final ArrayList<String> args = new ArrayList<>();
        File refFile = null;

        args.add("-I");
        args.add(new File(params.bam).getAbsolutePath());
        args.add("--" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME);
        args.add(new File(resourceDir + "HiSeq.20mb.1RG.table.gz").getAbsolutePath());
        args.add("-O");
        args.add(outFile.getAbsolutePath());
        if (params.reference != null) {
            refFile = new File(params.reference);
            args.add("-R");
            args.add(refFile.getAbsolutePath());
        }
        if (params.args != null) {
            Stream.of(params.args).forEach(arg -> args.add(arg));
        }

        runCommandLine(args);

        SamAssertionUtils.assertSamsEqual(outFile, new File(params.expectedFile), refFile);
    }

    @Test
    public void testOverfiltering() throws IOException {
        final File zeroRefBasesReadBam = new File(resourceDir, "NA12878.oq.read_consumes_zero_ref_bases.bam");
        final File outFile = GATKBaseTest.createTempFile("testReadThatConsumesNoReferenceBases", ".bam");
        final String[] args = new String[] {
                "--input", zeroRefBasesReadBam.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME, resourceDir + "NA12878.oq.gatk4.recal.gz",
                "--use-original-qualities",
                "--output", outFile.getAbsolutePath()
        };
        runCommandLine(args);
        //The expected output is actually the same as inputs for this read
        SamAssertionUtils.assertSamsEqual(outFile, zeroRefBasesReadBam);
    }

    @Test
    public void testMissingReadGroup() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -I " + hiSeqBamAligned +
                        " --" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME + " " + resourceDir + "HiSeq.20mb.1RG.table.missingRG.gz" +
                        " -O " + getDevNull(), 0,
                IllegalStateException.class);
        spec.executeTest("testMissingReadGroup", this);
    }

    @Test
    public void testemptyBqsrRecalFile() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -I " + hiSeqBamAligned +
                        " --" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME + " " + createTempFile("emptyBqsrRecal", "").toString() +
                        " -O " + getDevNull(), 0,
                UserException.class);
        spec.executeTest("testemptyBqsrRecalFile", this);
    }

    @Test
    public void testPRNoFailWithHighMaxCycle() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -I " + hiSeqBamAligned +
                        " --" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME + " " + resourceDir + "HiSeq.1mb.1RG.highMaxCycle.table.gz" +
                        " -O " + getDevNull(),
                Arrays.<String>asList());
        spec.executeTest("testPRNoFailWithHighMaxCycle", this);      //this just checks that the tool does not blow up
    }


    @Test
    public void testHelp() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -I " + hiSeqBamAligned +
                        " --help --" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME + " " + resourceDir + "HiSeq.1mb.1RG.highMaxCycle.table.gz" +
                        " -O " + getDevNull(),
                Arrays.<String>asList());
        spec.executeTest("testHelp", this);      //this just checks that the tool does not blow up
    }

    @Test
    public void testPRFailWithLowMaxCycle() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -I " + hiSeqBamAligned +
                        " --"  + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME+ " " + resourceDir + "HiSeq.1mb.1RG.lowMaxCycle.table.gz" +
                        " -O " + getDevNull(),
                0,
                UserException.class);
        spec.executeTest("testPRFailWithLowMaxCycle", this);
    }

    @Test
    public void testPRWithConflictingArguments_qqAndSQQ() throws IOException {
        // --quantize-quals and --static-quantized-quals shouldn't be able to be run in the same command
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " -I " + hiSeqBam +
                        " --" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME + " " + resourceDir + "HiSeq.20mb.1RG.table.gz" +
                        " --static-quantized-quals 9 --quantize-quals 4 " +
                        " -O " + getDevNull(),
                0,
                CommandLineException.class);
        spec.executeTest("testPRWithConflictingArguments_qqAndSQQ", this);
    }

    @Test
    public void testPRWithConflictingArguments_qqAndRDQ() throws IOException {
        // --quantize-quals and --static-quantized-quals shouldn't be able to be run in the same command
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " -I " + hiSeqBam +
                        " --" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME + " "  + resourceDir + "HiSeq.20mb.1RG.table.gz" +
                        " --round-down-quantized --quantize-quals 4 " +
                        " -O " + getDevNull(),
                0,
                CommandLineException.class);
        spec.executeTest("testPRWithConflictingArguments_qqAndSQQ", this);
    }
}
