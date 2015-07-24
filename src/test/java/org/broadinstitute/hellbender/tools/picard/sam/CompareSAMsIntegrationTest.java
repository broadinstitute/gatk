package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMFormatException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;

/**
 * Test for CompareSams tool.
 */
public final class CompareSAMsIntegrationTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/sam/CompareSAMs");

    @Override
    public String getTestedClassName() {
        return CompareSAMs.class.getSimpleName();
    }

    @DataProvider(name="testDataValidFormats")
    public Object[][] testDataValid() {
        return new Object[][]{
                {"NA12878.chr17_69k_70k.dictFix.bam", "NA12878.chr17_69k_70k.dictFix.bam", ValidationStringency.SILENT, true},
                {"unmapped_first.sam", "unmapped_second.sam", ValidationStringency.SILENT, false},
        };
    }

    @DataProvider(name="testDataInvalidFormats")
    public Object[][] testingDataInvalid() {
        return new Object[][]{
                {"NA12878.chr17_69k_70k.dictFix.bam", "NA12878.chr17_69k_70k.dictFix.bam", ValidationStringency.STRICT, false},
        };
    }

    private void compareSAMHelper(
            String fileName1, String fileName2, ValidationStringency stringency, boolean expectedResult) throws Exception {
        final File testFile1 = new File(TEST_DATA_DIR, fileName1);
        final File testFile2 = new File(TEST_DATA_DIR, fileName2);
        final String[] args = new String[]{
                testFile1.getPath(), testFile2.getPath(), "--VALIDATION_STRINGENCY", stringency.name(),
        };
        Assert.assertEquals(runCommandLine(args), expectedResult);
    }

    @Test(dataProvider="testDataValidFormats")
    public void testCompareSAM(
            String fileName1, String fileName2, ValidationStringency stringency, boolean expectedResult) throws Exception {
        compareSAMHelper(fileName1,  fileName2,  stringency, expectedResult);
    }

     @Test(dataProvider="testDataInvalidFormats", expectedExceptions=SAMFormatException.class)
        public void testCompareInvalidSAMFormats(
                String fileName1, String fileName2, ValidationStringency stringency, boolean expectedResult) throws Exception {
         compareSAMHelper(fileName1,  fileName2,  stringency, expectedResult);
    }
}
