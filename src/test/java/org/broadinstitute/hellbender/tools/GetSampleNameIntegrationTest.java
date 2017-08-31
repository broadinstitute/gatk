package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.junit.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

public class GetSampleNameIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/copynumber/allelic";
    private static final File NORMAL_BAM_FILE = new File(TEST_SUB_DIR, "collect-allelic-counts-normal.bam");
    private static final String TEST_SUB_DIR2 = publicTestDir + "org/broadinstitute/hellbender/tools";
    private static final File MS_BAD_BAM_FILE = new File(TEST_SUB_DIR2, "multi_sample_bam_header.bam");

    @Test
    public void testBasicUsage() throws IOException {

        final File outputFile = createTempFile("get-sample-name", ".txt");
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, NORMAL_BAM_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath()
        };
        runCommandLine(arguments);
        Assert.assertTrue(outputFile.exists());
        Assert.assertTrue(outputFile.length() > 0);
        Assert.assertTrue(Files.readAllLines(outputFile.toPath()).stream().count() == 1);
        Assert.assertTrue(Files.readAllLines(outputFile.toPath()).stream().filter(n -> n.equals("20")).count() == 1);
    }

    @Test(expectedExceptions = UserException.class)
    public void testMultiSampleBam() {
        final File outputFile = createTempFile("get-sample-name-ms", ".txt");
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, MS_BAD_BAM_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath()
        };
        runCommandLine(arguments);
    }
}
