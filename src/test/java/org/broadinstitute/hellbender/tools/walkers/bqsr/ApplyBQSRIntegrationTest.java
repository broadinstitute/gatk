package org.broadinstitute.hellbender.tools.walkers.bqsr;

import com.google.common.jimfs.Configuration;
import com.google.common.jimfs.Jimfs;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReaderFactory;
import java.nio.file.FileSystem;
import java.nio.file.Path;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.ApplyBQSRArgumentCollection;
import org.broadinstitute.hellbender.tools.ApplyBQSRUniqueArgumentCollection;
import org.broadinstitute.hellbender.transformers.BQSRReadTransformer;
import org.broadinstitute.hellbender.utils.collections.NestedIntegerArray;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.recalibration.RecalDatum;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationReport;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationTables;
import org.broadinstitute.hellbender.utils.recalibration.covariates.BQSRCovariateList;
import org.broadinstitute.hellbender.utils.recalibration.covariates.QualityScoreCovariate;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Stream;

public final class ApplyBQSRIntegrationTest extends CommandLineProgramTest {
    private static class ApplyBQSRTestData {
        final String bam;
        final String reference;
        final String outputExtension;
        final String args[];
        final String expectedFile;

        private ApplyBQSRTestData(String bam, String reference, String outputExtension, String args[], String expectedFile) {
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

    @Override
    public String getTestedClassName() {
        return ApplyBQSR.class.getSimpleName();
    }

    private final String resourceDir = getTestDataDir() + "/" + "BQSR" + "/";
    private final String hg18Reference = publicTestDir + "human_g1k_v37.chr17_1Mb.fasta";
    private final String hiSeqBam = resourceDir + "HiSeq.1mb.1RG.2k_lines.alternate.bam";
    private final String hiSeqCram = resourceDir + "HiSeq.1mb.1RG.2k_lines.alternate.cram";
    private final String hiSeqBamAligned = resourceDir + "HiSeq.1mb.1RG.2k_lines.alternate_allaligned.bam";
    private final String hiSeqCramAligned = resourceDir + "HiSeq.1mb.1RG.2k_lines.alternate_allaligned.cram";

    @DataProvider(name = "ApplyBQSRTest")
    public Object[][] createABQSRTestData() {
        List<Object[]> tests = new ArrayList<>();

        //Note: these outputs were created using GATK3
        tests.add(new Object[]{new ApplyBQSRTestData(hiSeqBam, null, ".bam", null, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.bam")});
        tests.add(new Object[]{new ApplyBQSRTestData(hiSeqBam, null, ".bam", new String[] {"-OQ"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.OQ.bam")});
        tests.add(new Object[]{new ApplyBQSRTestData(hiSeqBam, null, ".bam", new String[] {"--quantize-quals", "-1"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.qq-1.bam")});
        tests.add(new Object[]{new ApplyBQSRTestData(hiSeqBam, null, ".bam", new String[] {"--quantize-quals", "6"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.qq6.bam")});
        tests.add(new Object[]{new ApplyBQSRTestData(hiSeqBam, null, ".bam", new String[] {"--static-quantized-quals", "10", "--static-quantized-quals", "20", "--static-quantized-quals", "30"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.SQQ102030.bam")});
        tests.add(new Object[]{new ApplyBQSRTestData(hiSeqBam, null, ".bam", new String[] {"--static-quantized-quals", "10", "--static-quantized-quals", "20", "--static-quantized-quals", "30", "--round-down-quantized"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.SQQ102030RDQ.bam")});

        tests.add(new Object[]{new ApplyBQSRTestData(hiSeqBamAligned, null, ".bam", null, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.bam")});
        tests.add(new Object[]{new ApplyBQSRTestData(hiSeqBamAligned, null, ".bam", new String[] {"-OQ"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.OQ.bam")});
        tests.add(new Object[]{new ApplyBQSRTestData(hiSeqBamAligned, null, ".bam", new String[] {"--quantize-quals", "-1"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.qq-1.bam")});
        tests.add(new Object[]{new ApplyBQSRTestData(hiSeqBamAligned, null, ".bam", new String[] {"--quantize-quals", "6"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.qq6.bam")});
        tests.add(new Object[]{new ApplyBQSRTestData(hiSeqBamAligned, null, ".bam", new String[] {"--static-quantized-quals", "10", "--static-quantized-quals", "20", "--static-quantized-quals", "30"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.SQQ102030.bam")});
        tests.add(new Object[]{new ApplyBQSRTestData(hiSeqBamAligned, null, ".bam", new String[] {"--static-quantized-quals", "10", "--static-quantized-quals", "20", "--static-quantized-quals", "30", "--round-down-quantized"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.SQQ102030RDQ.bam")});

        //CRAM - input and output crams generated by direct conversion of the corresponding BAM test files with samtools 1.3
        tests.add(new Object[]{new ApplyBQSRTestData(hiSeqCram, hg18Reference, ".cram", new String[] {"--" + StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME, "true"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.cram")});
        tests.add(new Object[]{new ApplyBQSRTestData(hiSeqCramAligned, hg18Reference, ".cram", new String[] {"--quantize-quals", "6", "--" + StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME, "true"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.qq6.cram")});

        return tests.toArray(new Object[][]{});
    }

    @DataProvider(name = "MiniApplyBQSRTest")
    public Object[][] createMiniABQSRTestData() {
        List<Object[]> tests = new ArrayList<>();

        //Note: these outputs were created using GATK3
        tests.add(new Object[]{new ApplyBQSRTestData(hiSeqBam, null, ".bam", null, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.bam")});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ApplyBQSRTest")
    public void testApplyBQSRFile(ApplyBQSRTestData params) throws IOException {
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
            if (params.args != null) {
                Stream.of(params.args).forEach(arg -> args.add(arg));
            }

            runCommandLine(args);

            SamAssertionUtils.assertSamsEqual(outFile, new File(params.expectedFile), refFile);
        }
    }

    @Test(dataProvider = "MiniApplyBQSRTest")
    public void testApplyBQSRPath(ApplyBQSRTestData params) throws IOException {
        try (FileSystem jimfs = Jimfs.newFileSystem(Configuration.unix())) {
            final Path outPath = jimfs.getPath("applyBQSRTest"+params.outputExtension);

            final ArrayList<String> args = new ArrayList<>();
            Path refPath = null;

            args.add("-I");
            args.add(new File(params.bam).getAbsolutePath());
            args.add("--" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME);
            args.add(new File(resourceDir + "HiSeq.20mb.1RG.table.gz").getAbsolutePath());
            args.add("-O"); args.add(outPath.toUri().toString());
            if (params.reference != null) {
                File refFile = new File(params.reference);
                args.add("-R"); args.add(refFile.getAbsolutePath());
                refPath = refFile.toPath();
            }
            if (params.args != null) {
                Stream.of(params.args).forEach(arg -> args.add(arg));
            }

            runCommandLine(args);

            SamAssertionUtils.assertSamsEqual(outPath, new File(params.expectedFile).toPath(), refPath);
        }
    }

    @Test(dataProvider = "ApplyBQSRTest", groups={"bucket"})
    public void testApplyBQSRCloud(ApplyBQSRTestData params) throws IOException {
        // getTempFilePath also deletes the file on exit.
        final String outString = BucketUtils.getTempFilePath(getGCPTestStaging() + "tmp/testApplyBQSRCloud",  params.outputExtension);
        final Path outPath = BucketUtils.getPathOnGcs(outString);
        final ArrayList<String> args = new ArrayList<>();
        Path refPath = null;

        args.add("-I");
        args.add(new File(params.bam).getAbsolutePath());
        args.add("--" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME);
        args.add(new File(resourceDir + "HiSeq.20mb.1RG.table.gz").getAbsolutePath());
        args.add("-O");
        args.add(outString);
        if (params.reference != null) {
            File refFile = new File(params.reference);
            args.add("-R");
            args.add(refFile.getAbsolutePath());
            refPath = refFile.toPath();
        }
        if (params.args != null) {
            Stream.of(params.args).forEach(arg -> args.add(arg));
        }

        runCommandLine(args);

        SamAssertionUtils.assertSamsEqual(outPath, new File(params.expectedFile).toPath(), refPath);
    }

    @Test
    public void testMissingReadGroup() {
        final String inputBam = resourceDir + WGS_B37_CH20_1M_1M1K_BAM;
        final String knownSites = resourceDir + DBSNP_138_B37_CH20_1M_1M1K_VCF;

        // Base Recalibrator //
        final ArgumentsBuilder argsForGATKRecalibrator = new ArgumentsBuilder();
        final File recalTableOutput = createTempFile();
        final String readGroupToFilterOut = "20FUKAAXX100202.1";
        argsForGATKRecalibrator.addInput(inputBam);
        argsForGATKRecalibrator.addOutput(recalTableOutput);
        argsForGATKRecalibrator.add("read-filter", "ReadGroupBlackListReadFilter");
        argsForGATKRecalibrator.add(ReadFilterArgumentDefinitions.READ_GROUP_BLACK_LIST_LONG_NAME, "PU:" + readGroupToFilterOut);
        argsForGATKRecalibrator.add(BaseRecalibrator.KNOWN_SITES_ARG_FULL_NAME,  knownSites);
        argsForGATKRecalibrator.addReference(GCS_b37_CHR20_21_REFERENCE);
        argsForGATKRecalibrator.add(ApplyBQSRArgumentCollection.USE_ORIGINAL_QUALITIES_LONG_NAME, true);
        runCommandLine(argsForGATKRecalibrator, BaseRecalibrator.class.getSimpleName());

        // Apply BQSR //
        final ArgumentsBuilder argsForApplyBQSR = new ArgumentsBuilder();
        final File recalibratedBam = createTempFile("missingReadGroupTest", ".bam");
        argsForApplyBQSR.addInput(inputBam);
        argsForApplyBQSR.add(StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME, recalTableOutput);
        argsForApplyBQSR.addOutput(recalibratedBam);
        argsForApplyBQSR.add(ApplyBQSRArgumentCollection.ALLOW_MISSING_READ_GROUPS_LONG_NAME, true);
        // per the warp pipeline
        argsForApplyBQSR.add(ApplyBQSRUniqueArgumentCollection.STATIC_QUANTIZED_QUALS_LONG_NAME, 10);
        argsForApplyBQSR.add(ApplyBQSRUniqueArgumentCollection.STATIC_QUANTIZED_QUALS_LONG_NAME, 20);
        argsForApplyBQSR.add(ApplyBQSRUniqueArgumentCollection.STATIC_QUANTIZED_QUALS_LONG_NAME, 30);
        argsForApplyBQSR.add(ApplyBQSRUniqueArgumentCollection.STATIC_QUANTIZED_QUALS_LONG_NAME, 40);
        argsForApplyBQSR.add(ApplyBQSRArgumentCollection.USE_ORIGINAL_QUALITIES_LONG_NAME, true);

        runCommandLine(argsForApplyBQSR, ApplyBQSR.class.getSimpleName());

        // Validation //
        final ReadsPathDataSource originalReads = new ReadsPathDataSource(Path.of(inputBam));
        final Iterator<GATKRead> originalReadsIterator = originalReads.iterator();

        final ReadsPathDataSource recalibratedReads = new ReadsPathDataSource(recalibratedBam.toPath());
        final Iterator<GATKRead> recalibratedReadsIterator = recalibratedReads.iterator();

        // Counts the number of times (new qual) != (old qual) to detect the pathological case where none of the bases is recalibrated.
        int numDifferingQualBases = 0;

        while (recalibratedReadsIterator.hasNext()) {
            final GATKRead originalRead = originalReadsIterator.next();
            final GATKRead recalibratedRead = recalibratedReadsIterator.next();

            Assert.assertEquals(originalRead.getReadGroup(), recalibratedRead.getReadGroup());
            final SAMFileHeader originalBamHeader = originalReads.getHeader();

            final byte[] newQuals = recalibratedRead.getBaseQualities();
            final byte[] oldQuals = ReadUtils.getOriginalBaseQualities(originalRead);

            if (ReadUtils.getPlatformUnit(originalRead, originalBamHeader).equals(readGroupToFilterOut)) {
                // These are the read groups that are not in the recal table.
                final List<Byte> possibleQuals = Arrays.asList((byte) 2, (byte) 6, (byte) 10, (byte) 20, (byte) 30, (byte) 40);
                for (int i = 0; i < originalRead.getLength(); i++){
                    final byte newQual = newQuals[i];
                    final byte oldQual = oldQuals[i];
                    final int diff = Math.abs(newQual - oldQual);
                    // When the read group is missing we simply round to the closest bin (static bin in this particular test).
                    // But rounding is done in probability space, so let's set the allowable difference to be 10.
                    Assert.assertTrue(diff <= 10);
                    Assert.assertTrue(possibleQuals.contains(newQual), "newQual = " + newQual);
                    if (newQual != oldQual){
                        numDifferingQualBases++;
                    }
                }
            }
        }

        Assert.assertTrue(numDifferingQualBases > 0, "No recalibration was done");
        Assert.assertFalse(originalReadsIterator.hasNext(), "the original and recalibrated bam must have the same number of reads");
    }

    @Test
    public void testemptyBqsrRecalFile() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -I " + hiSeqBamAligned +
                        " --" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME + " " + createTempFile("emptyBqsrRecal", "").toString() +
                        " -O /dev/null", 0,
                UserException.class);
        spec.executeTest("testemptyBqsrRecalFile", this);
    }

    @Test
    public void testPRNoFailWithHighMaxCycle() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                        " -I " + hiSeqBamAligned +
                        " --" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME + " " + resourceDir + "HiSeq.1mb.1RG.highMaxCycle.table.gz" +
                        " -O /dev/null",
                Arrays.<String>asList());
        spec.executeTest("testPRNoFailWithHighMaxCycle", this);      //this just checks that the tool does not blow up
    }


    @Test
    public void testHelp() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -I " + hiSeqBamAligned +
                        " --help --" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME + " " + resourceDir + "HiSeq.1mb.1RG.highMaxCycle.table.gz" +
                        " -O /dev/null",
                Arrays.<String>asList());
        spec.executeTest("testHelp", this);      //this just checks that the tool does not blow up
    }

    @Test
    public void testPRFailWithLowMaxCycle() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                        " -I " + hiSeqBamAligned +
                        " --"  + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME+ " " + resourceDir + "HiSeq.1mb.1RG.lowMaxCycle.table.gz" +
                        " -O /dev/null",
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
                        " -O /dev/null",
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
                        " -O /dev/null",
                0,
                CommandLineException.class);
        spec.executeTest("testPRWithConflictingArguments_qqAndSQQ", this);
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
    public void testAddingPG() throws IOException {
        final File inFile = new File(resourceDir, "NA12878.oq.read_consumes_zero_ref_bases.bam");
        final File outFile = GATKBaseTest.createTempFile("testAddingPG", ".bam");
        final String[] args = new String[] {
                "--input", inFile.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME, resourceDir + "NA12878.oq.gatk4.recal.gz",
                "--use-original-qualities",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_SAM_PROGRAM_RECORD,
                "--output", outFile.getAbsolutePath()
        };
        runCommandLine(args);

        //The expected output is actually the same as inputs for this read (this ignores the PGs header)
        SamAssertionUtils.assertSamsEqual(outFile, inFile);

        //input has no GATK ApplyBQSR in headers
        Assert.assertNull(SamReaderFactory.makeDefault().open(inFile).getFileHeader().getProgramRecord("GATK ApplyBQSR"));

        //output has a GATK ApplyBQSR in headers
        Assert.assertNotNull(SamReaderFactory.makeDefault().open(outFile).getFileHeader().getProgramRecord("GATK ApplyBQSR"));
    }

    @Test
    public void testCustomCovariates() {
        final File inputBam = new File(resourceDir + WGS_B37_CH20_1M_1M1K_BAM);
        final File knownSites = new File(resourceDir + DBSNP_138_B37_CH20_1M_1M1K_VCF);

        //*** BaseRecalibrator with custom covariates ***//
        final File recalTableOutput = createTempFile("recal_table", ".txt");
        final ArgumentsBuilder recalibratorArgs = getArgumentsForRecalibration(inputBam, recalTableOutput,
                knownSites, Arrays.asList("RepeatLengthCovariate"));
        runCommandLine(recalibratorArgs, BaseRecalibrator.class.getSimpleName());

        //*** Apply BQSR with customCovariates ***//
        final File recalibratedBam = createTempFile("customCovariateTest", ".bam");
        final ArgumentsBuilder argsForApplyBQSR = getArgumentsForApplyBQSR(inputBam, recalibratedBam, recalTableOutput);
        runCommandLine(argsForApplyBQSR, ApplyBQSR.class.getSimpleName());

        //*** Baseline BaseRecalibrator ***//
        final File baselineRecalTableOutput = createTempFile("baseline_recal_table", ".txt");
        final ArgumentsBuilder recalibratorArgsBaseline = getArgumentsForRecalibration(inputBam, baselineRecalTableOutput,
                knownSites, Collections.<String>emptyList());
        runCommandLine(recalibratorArgsBaseline, BaseRecalibrator.class.getSimpleName());

        //*** Baseline ApplyBQSR ***//
        final File baselineRecalibratedBam = createTempFile("baseline_customCovariateTest", ".bam");
        final ArgumentsBuilder baselineArgsForApplyBQSR = getArgumentsForApplyBQSR(inputBam, baselineRecalibratedBam, baselineRecalTableOutput);
        runCommandLine(baselineArgsForApplyBQSR, ApplyBQSR.class.getSimpleName());

        //*** Validation ***//
        // Check that...the output is the same for the rest of the other covariates when this special covariates are not used.
        final RecalibrationReport evalRecalReport = new RecalibrationReport(recalTableOutput);
        final RecalibrationReport baselineRecalReport = new RecalibrationReport(baselineRecalTableOutput);

        // worry about porting this to recalibrator test class later
        final RecalibrationTables evalTables = evalRecalReport.getRecalibrationTables();
        final NestedIntegerArray<RecalDatum> evalReadGroupTable =  evalTables.getReadGroupTable();
        final NestedIntegerArray<RecalDatum> evalQualityScoreTable =  evalTables.getQualityScoreTable();
        final NestedIntegerArray<RecalDatum> evalContextTable =  evalTables.getTable(BQSRCovariateList.CONTEXT_COVARIATE_DEFAULT_INDEX);
        final NestedIntegerArray<RecalDatum> evalCycleTable =  evalTables.getTable(BQSRCovariateList.CYCLE_COVARIATE_DEFAULT_INDEX);
        final NestedIntegerArray<RecalDatum> evalRepeatLengthTable =  evalTables.getTable(BQSRCovariateList.CYCLE_COVARIATE_DEFAULT_INDEX + 1);


        final RecalibrationTables baselineTables = baselineRecalReport.getRecalibrationTables();
        final NestedIntegerArray<RecalDatum> baselineReadGroupTable =  baselineTables.getReadGroupTable();
        final NestedIntegerArray<RecalDatum> baselineQualityScoreTable =  baselineTables.getQualityScoreTable();
        // tsato: The dimension for context is something like 1000, this has got to be so wasteful.
        final NestedIntegerArray<RecalDatum> baselineContextTable =  baselineTables.getTable(BQSRCovariateList.CONTEXT_COVARIATE_DEFAULT_INDEX);
        final NestedIntegerArray<RecalDatum> baselineCycleTable =  baselineTables.getTable(BQSRCovariateList.CYCLE_COVARIATE_DEFAULT_INDEX);
        int d = 3; // I forget --- are cycle and context tables separate? But printed as the same in the file?

        // Read groups shouldn't have changed...
        int numReadGroups = new ReadsPathDataSource(inputBam.toPath()).getHeader().getReadGroups().size();
        long readGroupCount = 0;
        for (final NestedIntegerArray.Leaf<RecalDatum> leaf : baselineReadGroupTable.getAllLeaves()) {
            int[] keys = leaf.keys;
            RecalDatum baselineReadGroupDatum = leaf.value;
            RecalDatum evalReadGroupDatum = evalReadGroupTable.get2Keys(keys[0], BQSRReadTransformer.BASE_SUBSTITUTION_INDEX); // Ah, this is probably the right thing to do. 0 indexes the event --- get snp
            Assert.assertEquals(evalReadGroupDatum, baselineReadGroupDatum);
            readGroupCount += baselineReadGroupDatum.getNumObservations();
        }

        long reportedQualityCount = 0;
        for (final NestedIntegerArray.Leaf<RecalDatum> leaf : baselineQualityScoreTable.getAllLeaves()) {
            int[] keys = leaf.keys;
            RecalDatum baselineQualityScoreDatum = leaf.value;
            RecalDatum evalQualityScoreDatum = evalQualityScoreTable.get3Keys(keys[0], keys[1], BQSRReadTransformer.BASE_SUBSTITUTION_INDEX); // Ah, this is probably the right thing to do. 0 indexes the event --- get snp
            Assert.assertEquals(evalQualityScoreDatum, baselineQualityScoreDatum);
            reportedQualityCount += baselineQualityScoreDatum.getNumObservations();
        }



        // Check that the context covariate did not change after adding a custom covariate
        long contextCount = 0L;
        for (final NestedIntegerArray.Leaf<RecalDatum> leaf : baselineContextTable.getAllLeaves()) {
            int[] keys = leaf.keys;
            final RecalDatum baselineContextDatum = leaf.value;
            final RecalDatum evalContextDatum = evalContextTable.get4Keys(keys[0], keys[1], keys[2], BQSRReadTransformer.BASE_SUBSTITUTION_INDEX);
            Assert.assertEquals(evalContextDatum, baselineContextDatum);
            contextCount += baselineContextDatum.getNumObservations();
        }

        // Ditto cycle covariate
        long cycleCount = 0L;
        for (final NestedIntegerArray.Leaf<RecalDatum> leaf : baselineCycleTable.getAllLeaves()) {
            int[] keys = leaf.keys;
            final RecalDatum baselineCycleDatum = leaf.value;
            final RecalDatum evalCycleDatum = evalCycleTable.get4Keys(keys[0], keys[1], keys[2], BQSRReadTransformer.BASE_SUBSTITUTION_INDEX);
            Assert.assertEquals(evalCycleDatum, baselineCycleDatum);
            cycleCount += baselineCycleDatum.getNumObservations();
        }

        // Some basics checks on the repeat length covariates
        long repeatLengthCount = 0;
        for (final NestedIntegerArray.Leaf<RecalDatum> leaf : evalRepeatLengthTable.getAllLeaves()) {
            RecalDatum datum = leaf.value;
            repeatLengthCount += datum.getNumObservations();
        }

        // TODO: contextCount isn't the same as the rest. Investigate.
        Assert.assertEquals(repeatLengthCount, readGroupCount);
    }

    private ArgumentsBuilder getArgumentsForRecalibration(final File inputBam, final File outputBam,
                                                          final File knownSites, final List<String> customCovariates){
        final ArgumentsBuilder result = new ArgumentsBuilder();
        result.addInput(inputBam);
        result.addOutput(outputBam);
        result.add(BaseRecalibrator.KNOWN_SITES_ARG_FULL_NAME, knownSites);
        result.addReference(GCS_b37_CHR20_21_REFERENCE);
        result.add(ApplyBQSRArgumentCollection.USE_ORIGINAL_QUALITIES_LONG_NAME, true);
        for (String customCovariate : customCovariates){
            result.add(RecalibrationArgumentCollection.COVARIATES_LONG_NAME, customCovariate);
        }

        return result;
    }

    private ArgumentsBuilder getArgumentsForApplyBQSR(final File inputBam, final File outputBam, final File recalTable){
        final ArgumentsBuilder result = new ArgumentsBuilder();
        result.addInput(inputBam);
        result.add(StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME, recalTable);
        result.addOutput(outputBam);
        result.add(ApplyBQSRArgumentCollection.ALLOW_MISSING_READ_GROUPS_LONG_NAME, true);
        // per the warp pipeline
        result.add(ApplyBQSRUniqueArgumentCollection.STATIC_QUANTIZED_QUALS_LONG_NAME, 10);
        result.add(ApplyBQSRUniqueArgumentCollection.STATIC_QUANTIZED_QUALS_LONG_NAME, 20);
        result.add(ApplyBQSRUniqueArgumentCollection.STATIC_QUANTIZED_QUALS_LONG_NAME, 30);
        result.add(ApplyBQSRUniqueArgumentCollection.STATIC_QUANTIZED_QUALS_LONG_NAME, 40);
        result.add(ApplyBQSRArgumentCollection.USE_ORIGINAL_QUALITIES_LONG_NAME, true);
        return result;
    }
}
