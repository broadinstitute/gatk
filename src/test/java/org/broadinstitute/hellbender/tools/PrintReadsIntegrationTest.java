package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.read.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

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
        SamAssertionUtils.assertSamsEqual(ORIG_BAM, outFile);
    }

    @Test(dataProvider="testingDataCRAM")
    public void testFileToFileCRAM(String fileIn, String extOut, String reference) throws Exception {
        String samFile= fileIn;
        final File outFile = File.createTempFile(samFile + ".", extOut);
        outFile.deleteOnExit();
        final File ORIG_BAM = new File(TEST_DATA_DIR, samFile);
        final File refFile = new File(TEST_DATA_DIR, reference);
        final String[] args = new String[]{
                "--input" , ORIG_BAM.getAbsolutePath(),
                "--output", outFile.getAbsolutePath(),
                "-R", refFile.getAbsolutePath()
        };
        runCommandLine(args);
        SamAssertionUtils.assertSamsEqual(ORIG_BAM, outFile, refFile);
    }

    /**
     * This test is currently disabled because indexed CRAMs do not currently work.
     * TODO: Reenable this test once indexing with CRAM is fixed
     */
    @Test(dataProvider="testingDataIntervals")
    public void testFileToFileIntervals(String fileIn, String extOut, String reference, String intervals) throws Exception {
        String samFile= fileIn;
        final File outFile = File.createTempFile(samFile + ".", extOut);
        outFile.deleteOnExit();
        final File ORIG_BAM = new File(TEST_DATA_DIR, samFile);
        final File refFile = new File(TEST_DATA_DIR, reference);
        final File intervalFile = new File(TEST_DATA_DIR, intervals);
        final String[] args = new String[]{
                "--input" , ORIG_BAM.getAbsolutePath(),
                "--output", outFile.getAbsolutePath(),
                "-R", refFile.getAbsolutePath(),
                "-L", intervalFile.getAbsolutePath()
        };
        runCommandLine(args);
        SamAssertionUtils.assertSamsEqual(ORIG_BAM, outFile, refFile);
    }

    @DataProvider(name="testingData")
    public Object[][] testingData() {
        return new String[][]{
                {"print_reads.sam", ".sam"},
                {"print_reads.sam", ".bam"},
                {"print_reads.bam", ".sam"},
                {"print_reads.bam", ".bam"},
        };
    }

    @DataProvider(name="testingDataCRAM")
    public Object[][] testingDataCRAM() {
        return new String[][]{
                {"print_reads.sorted.cram", ".sam", "print_reads.fasta"},
                {"print_reads.sorted.cram", ".bam", "print_reads.fasta"},
                {"print_reads.sorted.cram", ".cram", "print_reads.fasta"},
                {"print_reads.sorted.sam", ".cram", "print_reads.fasta"},
                {"print_reads.sorted.bam", ".cram", "print_reads.fasta"}
        };
    }

    @DataProvider(name="testingDataIntervals")
    public Object[][] testingDataIntervals() {
        return new String[][]{
                {"print_reads.sorted.cram", ".sam", "print_reads.fasta", "print_reads.intervals"},
                {"print_reads.sorted.cram", ".bam", "print_reads.fasta", "print_reads.intervals"},
        };
    }

}