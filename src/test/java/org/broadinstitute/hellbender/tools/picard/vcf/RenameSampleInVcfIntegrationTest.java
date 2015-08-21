package org.broadinstitute.hellbender.tools.picard.vcf;

import com.jcraft.jsch.UserInfo;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

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
                "--INPUT", input.getAbsolutePath(),
                "--OUTPUT", outfile.getAbsolutePath(),
                "--OLD_SAMPLE_NAME", "NA12878",
                "--NEW_SAMPLE_NAME", "FRED"
        };

        runCommandLine(args);
        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");

    }

}
