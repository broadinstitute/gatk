package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public final class ReplaceSamHeaderIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = getTestDataDir();
    private static final File BAM_FILE1 = new File(TEST_DATA_DIR, "picard/sam/ReplaceSamHeader/first10.bam");
    private static final File SAM_FILE1 = new File(TEST_DATA_DIR, "picard/sam/ReplaceSamHeader/first10.sam");
    private static final File BAM_FILE2 = new File(TEST_DATA_DIR, "picard/sam/ReplaceSamHeader/file2.bam");
    private static final File SAM_FILE2 = new File(TEST_DATA_DIR, "picard/sam/ReplaceSamHeader/file2.sam");
    private static final File SAM_FILE2_QUERYSORT = new File(TEST_DATA_DIR, "picard/sam/ReplaceSamHeader/file2NameSorted.sam");

    public String getTestedClassName() {
        return ReplaceSamHeader.class.getSimpleName();
    }

    @Test
    public void testSamSamSam() throws Exception {
        final File outputFile = BaseTest.createTempFile("ReplaceSamHeader.testBasic", ".sam");
        runIt(SAM_FILE1, SAM_FILE2, outputFile);
        final File expectedOut = new File(TEST_DATA_DIR, "picard/sam/ReplaceSamHeader/samSam.sam");

        //TODO fails due to picard bug https://github.com/broadinstitute/picard/issues/284
//        SamAssertionUtils.assertSamsEqual(expectedOut, outputFile, ValidationStringency.DEFAULT_STRINGENCY);

        //using text comparison
        IntegrationTestSpec.assertEqualTextFiles(expectedOut, outputFile);
    }

    @Test(expectedExceptions = UserException.class)
    public void testSamSamSam_clashingSort() throws Exception {
        final File outputFile = BaseTest.createTempFile("ReplaceSamHeader.testBasic", ".sam");
        runIt(SAM_FILE1, SAM_FILE2_QUERYSORT, outputFile);
    }

    @Test(enabled = false) //disabled due to picard bug https://github.com/broadinstitute/picard/issues/284
    public void testSamSamBam() throws Exception {
        final File outputFile = BaseTest.createTempFile("ReplaceSamHeader.testBasic", ".bam");
        runIt(SAM_FILE1, SAM_FILE2, outputFile);
        final File expectedOut = new File(TEST_DATA_DIR, "picard/sam/ReplaceSamHeader/samSam.bam");

        SamAssertionUtils.assertSamsEqual(expectedOut, outputFile, ValidationStringency.DEFAULT_STRINGENCY);
    }

    @Test
    public void testSamBamSam() throws Exception {
        final File outputFile = BaseTest.createTempFile("ReplaceSamHeader.testBasic", ".sam");
        runIt(SAM_FILE1, BAM_FILE2, outputFile);
        final File expectedOut = new File(TEST_DATA_DIR, "picard/sam/ReplaceSamHeader/samBam.sam");

        //TODO fails due to picard bug https://github.com/broadinstitute/picard/issues/284
//        SamAssertionUtils.assertSamsEqual(expectedOut, outputFile, ValidationStringency.SILENT);

        //using text comparison
        IntegrationTestSpec.assertEqualTextFiles(expectedOut, outputFile);
    }

    @Test(enabled = false) //disabled due to picard bug https://github.com/broadinstitute/picard/issues/284
    public void testSamBamBam() throws Exception {
        final File outputFile = BaseTest.createTempFile("ReplaceSamHeader.testBasic", ".bam");
        runIt(SAM_FILE1, BAM_FILE2, outputFile);
        final File expectedOut = new File(TEST_DATA_DIR, "picard/sam/ReplaceSamHeader/samBam.bam");

        SamAssertionUtils.assertSamsEqual(expectedOut, outputFile, ValidationStringency.DEFAULT_STRINGENCY);
    }

    @Test(enabled = false) //disabled due to picard bug https://github.com/broadinstitute/picard/issues/284
    public void testBamBamBam() throws Exception {
        final File outputFile = BaseTest.createTempFile("ReplaceSamHeader.testBasic", ".bam");
        runIt(BAM_FILE1, BAM_FILE2, outputFile);
        final File expectedOut = new File(TEST_DATA_DIR, "picard/sam/ReplaceSamHeader/bamBam.bam");

        SamAssertionUtils.assertSamsEqual(expectedOut, outputFile, ValidationStringency.DEFAULT_STRINGENCY);
    }

    private void runIt(final File inputFile, final File inputHeaderFile, final File outputFile) {
        final List<String> args = new ArrayList<>(Arrays.asList(
                "--input", inputFile.getAbsolutePath(),
                "--HEADER", inputHeaderFile.getAbsolutePath(),
                "--output", outputFile.getAbsolutePath()
                )
        );
        runCommandLine(args);
    }

}
