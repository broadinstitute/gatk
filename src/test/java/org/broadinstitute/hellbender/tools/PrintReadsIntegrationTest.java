package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
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

        //TODO: Make this comparison non-lenient in all cases: https://github.com/broadinstitute/gatk/issues/1087
        if (extOut.equals(".cram")){
            SamAssertionUtils.assertSamsEqualLenient(ORIG_BAM, outFile, refFile);
        } else {
            SamAssertionUtils.assertSamsEqual(ORIG_BAM, outFile, refFile);
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
    public void testReadThatConsumesNoReferenceBases() {
        final File zeroRefBasesReadBam = new File(TEST_DATA_DIR, "read_consumes_zero_ref_bases.bam");
        final File outFile = BaseTest.createTempFile("testReadThatConsumesNoReferenceBases", ".bam");
        final String[] args = new String[] {
                "--input" , zeroRefBasesReadBam.getAbsolutePath(),
                "--output", outFile.getAbsolutePath()
        };
        // Make sure no exception is thrown given an input containing a read that consumes no reference bases
        runCommandLine(args);
    }
}