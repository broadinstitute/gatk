package org.broadinstitute.hellbender.tools.picard.fastq;

/**
 * Created by dkling on 8/24/15.
 */


import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;


public final class NormalizeFastaIntegrationTest extends CommandLineProgramTest {

    private static final File TEST_DATA_PATH = new File(getTestDataDir(), "picard/fasta/" );

    public String getTestedClassName() {
        return NormalizeFasta.class.getSimpleName();
    }

    @Test
    public void testNormalize () throws IOException {

        final File input = new File(TEST_DATA_PATH, "exampleFASTA.fasta");
        final File expectedFile = new File(TEST_DATA_PATH, "Normalizedexp.fasta");
        final File outfile = BaseTest.createTempFile("normalized", ".fasta");

        final String[] args = {
                "--INPUT", input.getAbsolutePath(),
                "--OUTPUT", outfile.getAbsolutePath(),
        };

        runCommandLine(args);
        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");

    }

}
