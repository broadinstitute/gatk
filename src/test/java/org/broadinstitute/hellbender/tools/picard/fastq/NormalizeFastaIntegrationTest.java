package org.broadinstitute.hellbender.tools.picard.fastq;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public final class NormalizeFastaIntegrationTest extends CommandLineProgramTest {

    private static final File TEST_DATA_PATH = new File(getTestDataDir(), "picard/fastq/NormalizeFasta");

    public String getTestedClassName() {
        return NormalizeFasta.class.getSimpleName();
    }

    @Test
    public void testNormalize() throws IOException {
        final File input = new File(TEST_DATA_PATH, "testfasta.fasta");
        final File expectedFile = new File(TEST_DATA_PATH, "testFASTA_WS_4.fasta");
        final File outfile = BaseTest.createTempFile("normalized", ".fasta");

        final String[] args = {
                "--input", input.getAbsolutePath(),
                "--output", outfile.getAbsolutePath(),
                "--TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE", "TRUE",
                "--LINE_LENGTH", "5",
        };

        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testFailNormalize() {
        final File input = new File(TEST_DATA_PATH, "testfasta.fasta");

        final String[] args = {
                "--input", input.getAbsolutePath(),
                "--output", input.getAbsolutePath(),       //will blow up
                "--TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE", "TRUE",
                "--LINE_LENGTH", "5",
        };

        runCommandLine(args);

    }
}

