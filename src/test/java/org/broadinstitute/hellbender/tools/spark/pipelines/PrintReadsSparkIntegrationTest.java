package org.broadinstitute.hellbender.tools.spark.pipelines;

import htsjdk.samtools.cram.build.CramIO;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

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

    @Test(dataProvider="testingData")
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

    @Test
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

    @Test
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

    @Test(expectedExceptions = UserException.IncompatibleSequenceDictionaries.class)
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

    @Test(dataProvider="testFileToFile_queryNameSorted", expectedExceptions = UserException.class)
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

    @Test(expectedExceptions = UserException.class)
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
    @Test
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

    @Test
    public void testReadFiltering_asIntegrationTestSpec() throws IOException {
        final File samWithOneMalformedRead = new File(getTestDataDir(), "print_reads_one_malformed_read.sam");

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " --" + StandardArgumentDefinitions.INPUT_LONG_NAME + " " + samWithOneMalformedRead.getCanonicalPath() +
                " --" + StandardArgumentDefinitions.OUTPUT_LONG_NAME + " " + "%s",
                Arrays.asList(new File(getTestDataDir(), "expected.print_reads_one_malformed_read.bam").getCanonicalPath())
        );
        spec.executeTest("testReadFiltering_asIntegrationTestSpec", this);
    }
}
