package org.broadinstitute.hellbender.tools.picard.sam;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public final class CalculateReadGroupChecksumIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/sam/CalculateReadGroupChecksum");

    @Test
    public void test() throws IOException {
        final File input = new File(TEST_DATA_DIR, "first5000a.bam");
        final File expectedFile = new File(TEST_DATA_DIR, "first5000a.CalculateReadGroupChecksum.txt");   //file created using picard 1.130
        final File outfile = BaseTest.createTempFile("testCalculateReadGroupChecksum", ".txt");
        final String[] args = new String[]{
                "--INPUT", input.getAbsolutePath(),
                "--OUTPUT", outfile.getAbsolutePath(),
        };
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile);
    }

    @Test
    public void testNoOutputName() throws IOException {
        final File input = new File(TEST_DATA_DIR, "first5000a.bam");
        final File expectedFile = new File(TEST_DATA_DIR, "first5000a.CalculateReadGroupChecksum.txt");   //file created using picard 1.130
        final File outfile =  new File(input.getParentFile(), CalculateReadGroupChecksum.getOutputFileName(input));
        outfile.deleteOnExit();
        final String[] args = new String[]{
                "--INPUT", input.getAbsolutePath()
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
