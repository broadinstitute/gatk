package org.broadinstitute.hellbender.tools;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public final class PrintReadsIntegrationTest extends CommandLineProgramTest{

    private static final File TEST_DATA_DIR = getTestDataDir();

    @Override
    public String getTestedClassName() {
        return PrintReads.class.getSimpleName();
    }

    @Test(dataProvider="testingData")
    public void testFileToFile(String fileIn, String extOut) throws Exception {
        String samFile= fileIn;
        final File outFile = BaseTest.createTempFile(samFile + ".", extOut);
        File ORIG_BAM = new File(TEST_DATA_DIR, samFile);
        final String[] args = new String[]{
                "--input" , ORIG_BAM.getAbsolutePath(),
                "--output", outFile.getAbsolutePath()
        };
        Assert.assertEquals(runCommandLine(args), null);
        SamAssertionUtils.assertSamsEqual(outFile, ORIG_BAM);
    }

    @Test(dataProvider="testingData")
    public void testFileToFileWithMD5(String fileIn, String extOut) throws Exception {
        String samFile= fileIn;
        final File outFile = BaseTest.createTempFile(samFile + ".", extOut);
        File ORIG_BAM = new File(TEST_DATA_DIR, samFile);
        final String[] args = new String[]{
                "--input" , ORIG_BAM.getAbsolutePath(),
                "--createOutputBamMD5", "true",
                "--output", outFile.getAbsolutePath()
        };
        Assert.assertEquals(runCommandLine(args), null);
        SamAssertionUtils.assertSamsEqual(outFile, ORIG_BAM);
        checkMD5asExpected(outFile);
    }

    @Test(dataProvider="testingDataCRAM")
    public void testFileToFileCRAMWithMD5(String fileIn, String extOut, String reference) throws Exception {
        String samFile= fileIn;
        final File outFile = BaseTest.createTempFile(samFile + ".", extOut);
        final File ORIG_BAM = new File(TEST_DATA_DIR, samFile);
        final File refFile = new File(TEST_DATA_DIR, reference);
        final String[] args = new String[]{
                "--input" , ORIG_BAM.getAbsolutePath(),
                "--createOutputBamMD5", "true",
                "--output", outFile.getAbsolutePath(),
                "-R", refFile.getAbsolutePath()
        };
        runCommandLine(args);

        //TODO: Make this comparison non-lenient in all cases: https://github.com/broadinstitute/gatk/issues/1087
        if (extOut.equals(".cram")){
            SamAssertionUtils.assertSamsEqualLenient(outFile, ORIG_BAM, refFile);
        } else {
            SamAssertionUtils.assertSamsEqual(outFile, ORIG_BAM, refFile);
        }
        checkMD5asExpected(outFile);
    }

    private void checkMD5asExpected(final File outFile) throws IOException {
        final File md5File = new File(outFile.getAbsolutePath() + ".md5");
        if (md5File.exists()) {
            md5File.deleteOnExit();
        }
        Assert.assertTrue(md5File.exists(), md5File + " does not exist");
        final String expectedMD5 = Utils.calculateFileMD5(outFile);
        final String actualMD5 = FileUtils.readFileToString(md5File);
        Assert.assertEquals(actualMD5, expectedMD5);
    }

    @Test(dataProvider="testingDataCRAM")
    public void testFileToFileCRAM(String fileIn, String extOut, String reference) throws Exception {
        String samFile= fileIn;
        final File outFile = BaseTest.createTempFile(samFile + ".", extOut);
        final File ORIG_BAM = new File(TEST_DATA_DIR, samFile);
        final File refFile = new File(TEST_DATA_DIR, reference);
        final String[] args = new String[]{
                "--input" , ORIG_BAM.getAbsolutePath(),
                "--output", outFile.getAbsolutePath(),
                "-R", refFile.getAbsolutePath()
        };
        runCommandLine(args);

        //TODO: Make this comparison non-lenient in all cases: https://github.com/broadinstitute/gatk/issues/1087
        if (extOut.equals(".cram")){
            SamAssertionUtils.assertSamsEqualLenient(outFile, ORIG_BAM, refFile);
        } else {
            SamAssertionUtils.assertSamsEqual(outFile, ORIG_BAM, refFile);
        }
    }

    @DataProvider(name="testingData")
    public Object[][] testingData() {
        return new String[][]{
                {"print_reads.sam", ".sam"},
                {"print_reads.sam", ".bam"},
                {"print_reads.bam", ".sam"},
                {"print_reads.bam", ".bam"},
                {"print_reads.sorted.queryname.sam", ".sam"},
                {"print_reads.sorted.queryname.sam", ".bam"},
                {"print_reads.sorted.queryname.bam", ".sam"},
                {"print_reads.sorted.queryname.bam", ".bam"},

        };
    }

    @DataProvider(name="testingDataCRAM")
    public Object[][] testingDataCRAM() {
        return new String[][]{
                {"print_reads.sorted.cram", ".sam", "print_reads.fasta"},
                {"print_reads.sorted.cram", ".bam", "print_reads.fasta"},
                {"print_reads.sorted.cram", ".cram", "print_reads.fasta"},
                {"print_reads.sorted.sam", ".cram", "print_reads.fasta"},
                {"print_reads.sorted.bam", ".cram", "print_reads.fasta"},
                {"print_reads.sorted.queryname.cram", ".sam", "print_reads.fasta"},
                {"print_reads.sorted.queryname.cram", ".bam", "print_reads.fasta"},

                //comparing queryname sorted cram files blows up: https://github.com/broadinstitute/gatk/issues/1271
//                {"print_reads.sorted.queryname.cram", ".cram", "print_reads.fasta"},
//                {"print_reads.sorted.queryname.sam", ".cram", "print_reads.fasta"},
//                {"print_reads.sorted.queryname.bam", ".cram", "print_reads.fasta"}
        };
    }

    @Test
    public void testReadThatConsumesNoReferenceBases() throws IOException {
        final File zeroRefBasesReadBam = new File(TEST_DATA_DIR, "read_consumes_zero_ref_bases.bam");
        final File outFile = BaseTest.createTempFile("testReadThatConsumesNoReferenceBases", ".bam");
        final String[] args = new String[] {
                "--input" , zeroRefBasesReadBam.getAbsolutePath(),
                "--output", outFile.getAbsolutePath()
        };
        // Make sure no exception is thrown given an input containing a read that consumes no reference bases
        runCommandLine(args);

        //Make sure we print the read, ie not lose it.
        SamAssertionUtils.assertSamsEqual(outFile, zeroRefBasesReadBam);
    }
}