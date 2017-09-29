package org.broadinstitute.hellbender.tools.spark.pipelines;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.filters.ReadLengthReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadNameReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public final class PrintReadsSparkIntegrationTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = getTestDataDir();

    @Override
    public String getTestedClassName() {
        return PrintReadsSpark.class.getSimpleName();
    }

    @DataProvider(name="testingData")
    public Object[][] testingData() {
        return new String[][]{
                {"print_reads.sorted.sam", ".sam", null},
                {"print_reads.sorted.sam", ".bam", null},
                {"print_reads.sorted.bam", ".sam", null},
                {"print_reads.sorted.bam", ".bam", null},
                {"print_reads.sorted.sam", ".cram", "print_reads.fasta"},
                {"print_reads.sorted.bam", ".cram", "print_reads.fasta"},
                {"print_reads.sorted.cram", ".sam", "print_reads.fasta"},
                {"print_reads.sorted.cram", ".bam", "print_reads.fasta"},
                {"print_reads.sorted.cram", ".cram", "print_reads.fasta"}
        };
    }

    @DataProvider
    public Object[][] gcsTestingData() {
        return new Object[][]{
            {"org/broadinstitute/hellbender/engine/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.10000000-10000020.with.unmapped.bam",
                ".bam", false, null},
            {"org/broadinstitute/hellbender/engine/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.10000000-10000020.with.unmapped.bam",
                ".bam", true, null},
        };
    }

    /**
     * Test the Spark code locally, including GCS access.
     *
     * For this to work, the settings in src/main/resources/core-site.xml must be correct,
     * and the project name and credential file it points to must be present.
     */
    @Test(dataProvider = "gcsTestingData", groups = "bucket")
    public void testGCSInputsAndOutputs(final String gcsInput, final String outputExtension,
        final boolean outputToGCS, final File expectedOutput) {
        final String gcsInputPath = getGCPTestInputPath() + gcsInput;
        final String outputPrefix = outputToGCS ? getGCPTestStaging() : "testGCSInputsAndOutputs";
        final String outputPath = BucketUtils.getTempFilePath(outputPrefix, outputExtension);

        final ArgumentsBuilder argBuilder = new ArgumentsBuilder();
        argBuilder.addArgument("input", gcsInputPath)
            .addArgument("output", outputPath);
        runCommandLine(argBuilder);
    }

    @Test(dataProvider="testingData", groups="spark")
    public void testFileToFile(String fileIn, String extOut, String reference) throws Exception {
        final File outFile = BaseTest.createTempFile(fileIn + ".", extOut);
        outFile.deleteOnExit();
        final File originalFile = new File(TEST_DATA_DIR, fileIn);
        final File refFile;
        final String[] args;
        if (reference == null) {
            refFile = null;
            args = new String[]{
                "--input", originalFile.getAbsolutePath(),
                "--output", outFile.getAbsolutePath(),
            };
        } else {
            refFile = new File(TEST_DATA_DIR, reference);
            args = new String[]{
                    "--input", originalFile.getAbsolutePath(),
                    "--output", outFile.getAbsolutePath(),
                    "-R", refFile.getAbsolutePath()
            };
        }
        runCommandLine(args);

        SamAssertionUtils.assertSamsEqual(outFile, originalFile, refFile);
    }

    @Test(groups = "spark")
    public void testCoordinateSorted() throws Exception {
        final File inBam = new File(getTestDataDir(), "print_reads.sorted.bam");
        final File outBam = BaseTest.createTempFile("print_reads_spark", ".bam");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(inBam.getCanonicalPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outBam.getCanonicalPath());

        this.runCommandLine(args.getArgsArray());

        SamAssertionUtils.assertSamsEqual(outBam, inBam);
    }

    @Test(groups = "spark")
    public void testCoordinateSortedInRegion() throws Exception {
        final File inBam = new File(getTestDataDir(), "print_reads.sorted.bam");
        final File expectedBam = new File(getTestDataDir(), "print_reads.sorted.chr1_1.bam");
        final File outBam = BaseTest.createTempFile("print_reads_spark", ".bam");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(inBam.getCanonicalPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outBam.getCanonicalPath());
        args.add("-L chr7:1-100 -XL chr7:2-100");

        this.runCommandLine(args.getArgsArray());

        SamAssertionUtils.assertSamsEqual(outBam, expectedBam);
    }

    @Test(expectedExceptions = UserException.IncompatibleSequenceDictionaries.class, groups="spark")
    public void testSequenceDictionaryValidation() throws Exception {
        final File inCram = new File(getTestDataDir(), "print_reads.sorted.cram");
        final File inRef = new File(getTestDataDir(), "print_reads.chr1only.fasta");
        final File outBam = BaseTest.createTempFile("print_reads_spark", ".bam");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(inCram.getCanonicalPath());
        args.add("-R");
        args.add(inRef.getCanonicalPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outBam.getCanonicalPath());

        this.runCommandLine(args.getArgsArray());
    }

    @DataProvider(name="testFileToFile_queryNameSorted")
    public Object[][] testFileToFile_queryNameSorted() {
        return new String[][]{
                {"print_reads.sorted.queryname.sam", ".sam", null},
                {"print_reads.sorted.queryname.sam", ".bam", null},
                {"print_reads.sorted.queryname.bam", ".sam", null},
                {"print_reads.sorted.queryname.bam", ".bam", null},
                {"print_reads.sorted.queryname.cram", ".sam", "print_reads.fasta"},
                {"print_reads.sorted.queryname.cram", ".bam", "print_reads.fasta"},
                {"print_reads.sorted.queryname.cram", ".cram", "print_reads.fasta"},
                {"print_reads.sorted.queryname.sam", ".cram", "print_reads.fasta"},
                {"print_reads.sorted.queryname.bam", ".cram", "print_reads.fasta"}
        };
    }

    @Test(dataProvider="testFileToFile_queryNameSorted", expectedExceptions = UserException.class, groups="spark")
    public void testFileToFile_queryNameSorted(String fileIn, String extOut, String reference) throws Exception {
        final File outFile = BaseTest.createTempFile(fileIn + ".", extOut);
        outFile.deleteOnExit();
        final File originalFile = new File(TEST_DATA_DIR, fileIn);
        final File refFile;
        final String[] args;
        if (reference == null) {
            refFile = null;
            args = new String[]{
                    "--input", originalFile.getAbsolutePath(),
                    "--output", outFile.getAbsolutePath(),
            };
        } else {
            refFile = new File(TEST_DATA_DIR, reference);
            args = new String[]{
                    "--input", originalFile.getAbsolutePath(),
                    "--output", outFile.getAbsolutePath(),
                    "-R", refFile.getAbsolutePath()
            };
        }
        runCommandLine(args);
        SamAssertionUtils.assertSamsEqual(outFile, originalFile, refFile);
    }

    @Test(expectedExceptions = UserException.class, groups = "spark")
    public void testNameSorted() throws Exception {
        final File inBam = new File(getTestDataDir(), "print_reads.bam");
        final File outBam = BaseTest.createTempFile("print_reads", ".bam");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(inBam.getCanonicalPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outBam.getCanonicalPath());

        this.runCommandLine(args.getArgsArray());

        SamAssertionUtils.assertSamsEqual(outBam, inBam);
    }

    /**
     * Test that PrintReadsSpark is correctly applying the WellformedReadFilter
     */
    @Test(groups = "spark")
    public void testReadFiltering() throws IOException {
        final File samWithOneMalformedRead = new File(getTestDataDir(), "print_reads_one_malformed_read.sam");
        final File outBam = BaseTest.createTempFile("print_reads_testReadFiltering", ".bam");

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(samWithOneMalformedRead.getCanonicalPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outBam.getCanonicalPath());

        runCommandLine(args.getArgsArray());
        SamAssertionUtils.assertSamsEqual(outBam, new File(getTestDataDir(), "expected.print_reads_one_malformed_read.bam"));
    }

    @Test(groups = "spark")
    public void testReadFiltering_asIntegrationTestSpec() throws IOException {
        final File samWithOneMalformedRead = new File(getTestDataDir(), "print_reads_one_malformed_read.sam");

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " --" + StandardArgumentDefinitions.INPUT_LONG_NAME + " " + samWithOneMalformedRead.getCanonicalPath() +
                " --" + StandardArgumentDefinitions.OUTPUT_LONG_NAME + " " + "%s",
                Arrays.asList(new File(getTestDataDir(), "expected.print_reads_one_malformed_read.bam").getCanonicalPath())
        );
        spec.executeTest("testReadFiltering_asIntegrationTestSpec", this);
    }

    @Test(expectedExceptions = UserException.MissingReference.class)
    public void testNonExistentReference() throws Exception {
        final File inCram = new File(TEST_DATA_DIR, "print_reads.sorted.cram");
        final File outCram = BaseTest.createTempFile("print_reads_bad_reference", ".cram");

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(inCram.getCanonicalPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outCram.getCanonicalPath());
        args.add("-R");
        args.add(BaseTest.getSafeNonExistentFile("Nonexistent.fasta").getCanonicalPath());

        runCommandLine(args.getArgsArray());
    }

    @DataProvider(name = "UnmappedReadInclusionTestData")
    public Object[][] unmappedReadInclusionTestData() {
        // This bam has mapped reads from various contigs, plus a few unmapped reads with no mapped mate
        final File unmappedBam = new File(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1_with_unmapped.bam");

        // This is a snippet of the CEUTrio.HiSeq.WGS.b37.NA12878 bam from large, with mapped reads
        // from chromosome 20 (with one mapped read having an unmapped mate), plus several unmapped
        // reads with no mapped mate.
        final File ceuSnippet = new File(publicTestDir + "org/broadinstitute/hellbender/engine/CEUTrio.HiSeq.WGS.b37.NA12878.snippet_with_unmapped.bam");
        final File ceuSnippetCram = new File(publicTestDir + "org/broadinstitute/hellbender/engine/CEUTrio.HiSeq.WGS.b37.NA12878.snippet_with_unmapped.cram");

        return new Object[][] {
                { unmappedBam, null, Arrays.asList("unmapped"), Arrays.asList("u1", "u2", "u3", "u4", "u5") },
                // The same interval as above in an intervals file
                { unmappedBam, null, Arrays.asList(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1_unmapped.intervals"), Arrays.asList("u1", "u2", "u3", "u4", "u5") },
                { unmappedBam, null, Arrays.asList("1:200-300", "unmapped"), Arrays.asList("a", "b", "c", "u1", "u2", "u3", "u4", "u5") },
                { unmappedBam, null, Arrays.asList("1:200-300", "4:700-701", "unmapped"), Arrays.asList("a", "b", "c", "k", "u1", "u2", "u3", "u4", "u5") },
                // The same intervals as above in an intervals file
                { unmappedBam, null, Arrays.asList(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1_unmapped2.intervals"), Arrays.asList("a", "b", "c", "k", "u1", "u2", "u3", "u4", "u5") },
                { ceuSnippet, null, Arrays.asList("unmapped"), Arrays.asList("g", "h", "h", "i", "i") },
                { ceuSnippet, null, Arrays.asList("20:10000009-10000011"), Arrays.asList("a", "b", "c", "d", "e") },
                { ceuSnippet, null, Arrays.asList("20:10000009-10000011", "unmapped"), Arrays.asList("a", "b", "c", "d", "e", "g", "h", "h", "i", "i") },
                { ceuSnippet, null, Arrays.asList("20:10000009-10000013", "unmapped"), Arrays.asList("a", "b", "c", "d", "e", "f", "f", "g", "h", "h", "i", "i") },
                { ceuSnippetCram, b37_reference_20_21, Arrays.asList("unmapped"), Arrays.asList("g", "h", "h", "i", "i") },
                { ceuSnippetCram, b37_reference_20_21, Arrays.asList("20:10000009-10000011", "unmapped"), Arrays.asList("a", "b", "c", "d", "e", "g", "h", "h", "i", "i") },
                { ceuSnippetCram, b37_reference_20_21, Arrays.asList("20:10000009-10000013", "unmapped"), Arrays.asList("a", "b", "c", "d", "e", "f", "f", "g", "h", "h", "i", "i") }
        };
    }

    @Test(dataProvider = "UnmappedReadInclusionTestData")
    public void testUnmappedReadInclusion( final File input, final String reference, final List<String> intervalStrings, final List<String> expectedReadNames ) {
        final File outFile = createTempFile("testUnmappedReadInclusion", ".bam");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("-I"); args.add(input.getAbsolutePath());
        args.add("-O"); args.add(outFile.getAbsolutePath());
        for ( final String intervalString : intervalStrings ) {
            args.add("-L"); args.add(intervalString);
        }
        if ( reference != null ) {
            args.add("-R"); args.add(reference);
        }

        runCommandLine(args);

        try ( final ReadsDataSource outputReadsSource = new ReadsDataSource(outFile.toPath()) ) {
            final List<GATKRead> actualReads = new ArrayList<>();
            for ( final GATKRead read : outputReadsSource ) {
                actualReads.add(read);
            }

            if (actualReads.size() != expectedReadNames.size()) {
                System.out.println("actual: " + actualReads);
                System.out.println("expectedReadNames: " + expectedReadNames);
            }
            Assert.assertEquals(actualReads.size(), expectedReadNames.size(), "Wrong number of reads output");

            for ( int readNumber = 0; readNumber < actualReads.size(); ++readNumber ) {
                Assert.assertEquals(actualReads.get(readNumber).getName(), expectedReadNames.get(readNumber), "Unexpected read name");
            }
        }
    }

    @DataProvider(name="readFilterTestData")
    public Object[][] testReadFilterData() {
        return new Object[][]{
                {"print_reads_one_malformed_read.sam", null, ".sam", Collections.emptyList(), 7},
                {"print_reads_one_malformed_read.sam", null, ".sam", Arrays.asList("--disableToolDefaultReadFilters"), 8},
                {"print_reads_one_malformed_read.sam", null, ".sam",
                        Arrays.asList("--disableReadFilter", "WellformedReadFilter"), 8},
                {"print_reads.sorted.sam", null, ".sam",
                        Arrays.asList(
                                "--readFilter", ReadNameReadFilter.class.getSimpleName(),
                                "--readName", "both_reads_align_clip_adapter"),
                        2},
                {"print_reads.sorted.sam", null, ".sam",
                        Arrays.asList(
                                "--RF", ReadLengthReadFilter.class.getSimpleName(),
                                "--minReadLength", "100",
                                "--maxReadLength", "200"),
                        8},
                {"print_reads.sorted.sam", null, ".sam",
                        Arrays.asList(
                                "--RF", ReadLengthReadFilter.class.getSimpleName(),
                                "--minReadLength", "1",
                                "--maxReadLength", "10"),
                        0},
                {"print_reads.sorted.sam", null, ".sam",
                        Arrays.asList(
                                "--readFilter", ReadNameReadFilter.class.getSimpleName(),
                                "--readName", "both_reads_align_clip_adapter",
                                "--RF", ReadLengthReadFilter.class.getSimpleName(),
                                "--minReadLength", "100",
                                "--maxReadLength", "101"),
                        2},
                {"print_reads.sorted.bam", null, ".sam", Arrays.asList("--disableToolDefaultReadFilters"), 8},
                {"print_reads.sorted.bam", null, ".sam",
                        Arrays.asList(
                                "--readFilter", ReadNameReadFilter.class.getSimpleName(),
                                "--readName", "both_reads_align_clip_adapter"),
                        2},
                {"print_reads.sorted.bam", null, ".sam",
                        Arrays.asList(
                                "--RF", ReadLengthReadFilter.class.getSimpleName(),
                                "--minReadLength", "100",
                                "--maxReadLength", "101"),
                        8},
                {"print_reads.sorted.bam", null, ".sam",
                        Arrays.asList(
                                "--readFilter", ReadNameReadFilter.class.getSimpleName(),
                                "--readName", "both_reads_align_clip_adapter",
                                "--RF", ReadLengthReadFilter.class.getSimpleName(),
                                "--minReadLength", "100",
                                "--maxReadLength", "101"),
                        2},
                {"print_reads.sorted.cram", "print_reads.fasta", ".sam",
                        Arrays.asList(
                                "--readFilter", ReadNameReadFilter.class.getSimpleName(),
                                "--readName", "both_reads_align_clip_adapter",
                                "--RF", ReadLengthReadFilter.class.getSimpleName(),
                                "--minReadLength", "100",
                                "--maxReadLength", "101"),
                        2},
        };
    }

    @Test(dataProvider = "readFilterTestData", groups = "spark")
    public void testReadFilters(
            final String input,
            final String reference,
            final String extOut,
            final List<String> inputArgs,
            final int expectedCount) throws IOException
    {
        final File outFile = createTempFile("testReadFilter", extOut);

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("-I"); args.add(new File(TEST_DATA_DIR, input).getAbsolutePath());
        args.add("-O"); args.add(outFile.getAbsolutePath());
        if ( reference != null ) {
            args.add("-R"); args.add(new File(TEST_DATA_DIR, reference).getAbsolutePath());
        }
        for (final String filter : inputArgs) {
            args.add(filter);
        }

        runCommandLine(args);


        SamReaderFactory factory = SamReaderFactory.makeDefault();
        if (reference != null) {
            factory = factory.referenceSequence(new File(TEST_DATA_DIR, reference));
        }
        int count = 0;
        try (final SamReader reader = factory.open(outFile)) {
            Iterator<SAMRecord> it = reader.iterator();
            while (it.hasNext()) {
                SAMRecord rec = it.next();
                count++;
            }
        }
        Assert.assertEquals(count, expectedCount);
    }

}
