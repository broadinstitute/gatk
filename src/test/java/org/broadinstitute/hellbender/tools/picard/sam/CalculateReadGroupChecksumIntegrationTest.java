package org.broadinstitute.hellbender.tools.picard.sam;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public final class CalculateReadGroupChecksumIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/sam/CalculateReadGroupChecksum");

    @DataProvider(name="calculateReadGroupData")
    public Object[][] calculateReadGroupData() {
        return new Object[][] {
                { "first5000a.bam", null, "first5000a.CalculateReadGroupChecksum.txt"},     //   file created using picard 1.130
                { "multigroup_valid.bam", null, "multigroup_valid_result.txt"},           // result generated with picard 1.140
                { "multigroup_valid.cram", "basic.fasta", "multigroup_valid_result.txt"} // result generated with picard 1.140
        };
    }

    @Test(dataProvider="calculateReadGroupData")
    public void test(final String inputFileName, final String referenceFileName, final String expectedFileName) throws IOException {
        final File inputFile = new File(TEST_DATA_DIR, inputFileName);
        final File outfile = BaseTest.createTempFile("testCalculateReadGroupChecksum", ".txt");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--input"); args.add(inputFile.getAbsolutePath());
        args.add("--output"); args.add(outfile.getAbsolutePath());
        if (null != referenceFileName) {
            args.add("--R");
            args.add(new File(TEST_DATA_DIR, referenceFileName).getAbsolutePath());
        }
        runCommandLine(args.getArgsArray());

        final File expectedFile = new File(TEST_DATA_DIR, expectedFileName);
        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile);
    }

    @Test
    public void testNoOutputName() throws IOException {
        final File input = new File(TEST_DATA_DIR, "first5000a.bam");
        final File expectedFile = new File(TEST_DATA_DIR, "first5000a.CalculateReadGroupChecksum.txt");   //file created using picard 1.130
        final File outfile =  new File(input.getParentFile(), CalculateReadGroupChecksum.getOutputFileName(input));
        outfile.deleteOnExit();
        final String[] args = new String[]{
                "--input", input.getAbsolutePath()
        };
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile);
    }

    @Test
    public void testName() {
        final String name = "first5000a.bam";
        final File f = new File(TEST_DATA_DIR, name);
        final String outputFileName = CalculateReadGroupChecksum.getOutputFileName(f);
        Assert.assertEquals(outputFileName, name + ".read_group_md5");
    }
}
