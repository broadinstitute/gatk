package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.read.SamAssertionUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public class PrintReadsIntegrationTest extends CommandLineProgramTest{

    private static final File TEST_DATA_DIR = getTestDataDir();

    @Override
    public String getTestedClassName() {
        return PrintReads.class.getSimpleName();
    }


    @Test(dataProvider="testingData")
    public void testFileToFile(String fileIn, String extOut) throws Exception {
        String samFile= fileIn;
        final File outFile = File.createTempFile(samFile + ".", extOut);
        outFile.deleteOnExit();
        File ORIG_BAM = new File(TEST_DATA_DIR, samFile);
        final String[] args = new String[]{
                "--input" , ORIG_BAM.getAbsolutePath(),
                "--output", outFile.getAbsolutePath()
        };
        Assert.assertEquals(runCommandLine(args), null);
        SamAssertionUtils.assertSamsEqual(ORIG_BAM, outFile);
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

}