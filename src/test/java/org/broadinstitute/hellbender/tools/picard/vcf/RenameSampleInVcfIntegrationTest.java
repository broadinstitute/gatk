package org.broadinstitute.hellbender.tools.picard.vcf;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public final class RenameSampleInVcfIntegrationTest extends CommandLineProgramTest {

    private static final File TEST_DATA_PATH = new File(getTestDataDir(), "picard/vcf/RenameSampleInVcf");

    public String getTestedClassName() {
        return RenameSampleInVcf.class.getSimpleName();
    }

    @Test
    public void testRename () throws IOException {

        final File input = new File(TEST_DATA_PATH, "input.vcf");
        final File expectedFile = new File(TEST_DATA_PATH, "expected_output.vcf");
        final File outfile = BaseTest.createTempFile("renamed", ".vcf");

        final String[] args = {
                "--input", input.getAbsolutePath(),
                "--output", outfile.getAbsolutePath(),
                "--OLD_SAMPLE_NAME", "NA12878",
                "--NEW_SAMPLE_NAME", "FRED"
        };

        runCommandLine(args);
        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");

    }

}
