package org.broadinstitute.hellbender.cmdline;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.PicardNonZeroExitException;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;


// Test that we can successfully run Picard tools.
public final class PicardCommandLineProgramExecutorIntegrationTest extends CommandLineProgramTest {

    // Use the Picard tool "NormalizeFasta" as a simple test case.
    private static final File TEST_DATA_PATH = new File(getTestDataDir(), "picard/fastq/NormalizeFasta");

    @Override
    public String getTestedClassName() {
        return picard.reference.NormalizeFasta.class.getSimpleName();
    }

    @Override
    public List<String> injectDefaultVerbosity(final List<String> args) {
        // override/suppress GATK-specific argument injection
        return args;
    }

    @Test
    public void testPicardNormalizeFasta() throws IOException {
        final File input = new File(TEST_DATA_PATH, "testfasta.fasta");
        final File expectedFile = new File(TEST_DATA_PATH, "testFASTA_WS_4.fasta");
        final File outfile = createTempFile("normalized", ".fasta");

        final String[] args = {
                "-I", input.getAbsolutePath(),
                "-O", outfile.getAbsolutePath(),
                "--TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE", "TRUE",
                "--LINE_LENGTH", "5",
        };

        Assert.assertEquals(runCommandLine(args), 0);

        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile);
    }

    @Test(expectedExceptions=PicardNonZeroExitException.class)
    public void testPicardNormalizeFastaWithBadArgs() throws IOException {
        final File input = new File(TEST_DATA_PATH, "testfasta.fasta");
        final File outfile = createTempFile("normalized", ".fasta");

        // Use GATK-style lower case argument names, which are rejected by Picard
        // because it uses upper cased argument names (--INPUT/--OUTPUT)
        final String[] args = {
                "--input", input.getAbsolutePath(),
                "--output", outfile.getAbsolutePath(),
                "--TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE", "TRUE",
                "--LINE_LENGTH", "5",
        };

        // The tool should fail due to the lower case argument name; Picard doesn't throw but
        // does return a non-zero return code
        Assert.assertNotEquals(runCommandLine(args), 0);
    }

}
