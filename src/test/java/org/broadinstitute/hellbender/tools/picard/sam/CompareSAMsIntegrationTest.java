package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFormatException;
import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

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
                {"multigroup_valid.cram", "multigroup_valid.cram", "basic.fasta", ValidationStringency.SILENT, true},
                {"NA12878.chr17_69k_70k.dictFix.bam", "NA12878.chr17_69k_70k.dictFix.bam", null, ValidationStringency.SILENT, true},
                {"unmapped_first.sam", "unmapped_second.sam", null, ValidationStringency.SILENT, false},
        };
    }

    @DataProvider(name="testDataInvalidFormats")
    public Object[][] testingDataInvalid() {
        return new Object[][]{
                {"NA12878.chr17_69k_70k.dictFix.bam", "NA12878.chr17_69k_70k.dictFix.bam", ValidationStringency.STRICT, false},
        };
    }

    private void compareSAMHelper(
            final String fileName1,
            final String fileName2,
            final String referenceFileName,
            final ValidationStringency stringency,
            final boolean expectedResult) throws Exception {
        final File testFile1 = new File(TEST_DATA_DIR, fileName1);
        final File testFile2 = new File(TEST_DATA_DIR, fileName2);
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add(testFile1.getPath());
        args.add(testFile2.getPath());
        if (null != referenceFileName) {
            args.add("--R");
            args.add(new File(TEST_DATA_DIR, referenceFileName).getAbsoluteFile());
        }
        args.add("--VALIDATION_STRINGENCY"); args.add(stringency.name());
        Assert.assertEquals(runCommandLine(args.getArgsArray()), expectedResult);
    }

    @Test(dataProvider="testDataValidFormats")
    public void testCompareSAM(
            final String fileName1,
            final String fileName2,
            final String referenceFileName,
            final ValidationStringency stringency,
            final boolean expectedResult) throws Exception {
        compareSAMHelper(fileName1, fileName2, referenceFileName, stringency, expectedResult);
    }

     @Test(dataProvider="testDataInvalidFormats", expectedExceptions=SAMFormatException.class)
        public void testCompareInvalidSAMFormats(
                String fileName1, String fileName2, ValidationStringency stringency, boolean expectedResult) throws Exception {
         compareSAMHelper(fileName1,  fileName2,  null, stringency, expectedResult);
    }

    @Test(expectedExceptions = SAMException.class)
    public void testCompareSAMMissingFile() throws Exception {
        compareSAMHelper("foo", "bar", null, ValidationStringency.DEFAULT_STRINGENCY, false);
    }

}
