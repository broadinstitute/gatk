package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public final class PrintReadsHeaderIntegrationTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = getTestDataDir();

    @Test
    public void testPrintReadsHeader() throws IOException {
        final File inFile = new File(TEST_DATA_DIR, "print_reads.sorted.bam");
        final File outFile = createTempFile("PrintReadsHeaderIntegrationTest_testPrintBAMHeader", ".txt");
        final String[] args = new String[] {
                "--input" , inFile.getAbsolutePath(),
                "--output", outFile.getAbsolutePath()
        };
        runCommandLine(args);

        //Make sure contents are the same
        IntegrationTestSpec.assertEqualTextFiles(outFile, new File(TEST_DATA_DIR, "print_reads.sorted.bam.header.txt"));
    }
}