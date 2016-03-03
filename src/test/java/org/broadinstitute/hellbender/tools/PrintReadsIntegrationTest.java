package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.cram.build.CramIO;
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
import java.util.ArrayList;

public final class PrintReadsIntegrationTest extends CommandLineProgramTest{

    private static final File TEST_DATA_DIR = getTestDataDir();

    @Override
    public String getTestedClassName() {
        return PrintReads.class.getSimpleName();
    }

    public void doFileToFile(String fileIn, String extOut, String reference, boolean testMD5) throws Exception {
        String samFile = fileIn;
        final File outFile = BaseTest.createTempFile(samFile + ".", extOut);
        final File ORIG_BAM = new File(TEST_DATA_DIR, samFile);

        final ArrayList<String> args = new ArrayList<>();
        args.add("--input"); args.add(ORIG_BAM.getAbsolutePath());
        args.add("--output"); args.add(outFile.getAbsolutePath());
        if (reference != null) {
            final File refFile = new File(TEST_DATA_DIR, reference);
            args.add("-R"); args.add(refFile.getAbsolutePath());
        };
        if (testMD5) {
            args.add("--createOutputBamMD5");
            args.add("true");
        }
        runCommandLine(args);

        if (fileIn.endsWith(CramIO.CRAM_FILE_EXTENSION) || extOut.equals(CramIO.CRAM_FILE_EXTENSION)) {
            final File refFile = new File(TEST_DATA_DIR, reference);
            if (!fileIn.endsWith(CramIO.CRAM_FILE_EXTENSION)) {
                // TODO: We still need lenient comparison due to https://github.com/samtools/htsjdk/issues/455
                // (bam->cram conversion strips out NM/MD tags)
                SamAssertionUtils.assertSamsEqualLenient(outFile, ORIG_BAM, refFile);
            }
            else {
                SamAssertionUtils.assertSamsEqual(outFile, ORIG_BAM, refFile);
            }
        } else {
            SamAssertionUtils.assertSamsEqual(outFile, ORIG_BAM);
        }

        if (testMD5) {
            checkMD5asExpected(outFile);
        }
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

    @Test(dataProvider="testingData")
    public void testFileToFile(String fileIn, String extOut, String reference) throws Exception {
        doFileToFile(fileIn, extOut, reference, false);
    }

    @Test(dataProvider="testingData")
    public void testFileToFileWithMD5(String fileIn, String extOut, String reference) throws Exception {
        doFileToFile(fileIn, extOut, reference, true);
    }

    @DataProvider(name="testingData")
    public Object[][] testingData() {
        return new String[][]{
                {"print_reads.sam", ".sam", null},
                {"print_reads.sam", ".bam", null},
                {"print_reads.sam", ".cram", "print_reads.fasta"},
                {"print_reads.bam", ".sam", null},
                {"print_reads.bam", ".bam", null},
                {"print_reads.bam", ".cram", "print_reads.fasta"},
                {"print_reads.cram", ".sam", "print_reads.fasta"},
                {"print_reads.cram", ".bam", "print_reads.fasta"},
                {"print_reads.cram", ".cram", "print_reads.fasta"},

                {"print_reads.sorted.sam", ".sam", null},
                {"print_reads.sorted.sam", ".bam", null},
                {"print_reads.sorted.sam", ".cram", "print_reads.fasta"},
                {"print_reads.sorted.bam", ".sam", null},
                {"print_reads.sorted.bam", ".bam", null},
                {"print_reads.sorted.bam", ".cram", "print_reads.fasta"},
                {"print_reads.sorted.cram", ".sam", "print_reads.fasta"},
                {"print_reads.sorted.cram", ".bam", "print_reads.fasta"},
                {"print_reads.sorted.cram", ".cram", "print_reads.fasta"},

                {"print_reads.sorted.queryname.sam", ".sam", null},
                {"print_reads.sorted.queryname.sam", ".bam", null},
                {"print_reads.sorted.queryname.sam", ".cram", "print_reads.fasta"},
                {"print_reads.sorted.queryname.bam", ".sam", null},
                {"print_reads.sorted.queryname.bam", ".bam", null},
                {"print_reads.sorted.queryname.bam", ".cram", "print_reads.fasta"},
                {"print_reads.sorted.queryname.cram", ".sam", "print_reads.fasta"},
                {"print_reads.sorted.queryname.cram", ".bam", "print_reads.fasta"},
                {"print_reads.sorted.queryname.cram", ".cram", "print_reads.fasta"},

                //test queryname-sorted crams with multiref containers in GATK:
                //print_reads.sorted.queryname_htsjdk_2.1.0.cram was generated from print_reads.sam
                //using gatk4 PrintReads/htsjdk.2.1.0, which includes changes to support
                //multireference containers
                {"print_reads.sorted.queryname.htsjdk-2.1.0.cram", ".cram", "print_reads.fasta"},
                {"print_reads.sorted.queryname.htsjdk-2.1.0.cram", ".sam", "print_reads.fasta"}
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