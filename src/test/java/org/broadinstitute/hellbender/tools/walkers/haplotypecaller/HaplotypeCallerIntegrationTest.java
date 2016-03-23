package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;

public class HaplotypeCallerIntegrationTest extends CommandLineProgramTest {

    public static final String TEST_FILES_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/haplotypecaller/";

    @Test
    public void testVCFMode() throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testVCFMode", ".vcf");
        final File expected = new File(TEST_FILES_DIR + "expected.testVCFMode.vcf");

        final String[] args = new String[] {
            "-I", NA12878_20_21_WGS_bam,
            "-R", b37_reference_20_21,
            "-O", output.getAbsolutePath(),
            "--" + StandardArgumentDefinitions.VERBOSITY_NAME, Log.LogLevel.INFO.toString()
        };

        runCommandLine(args);
        IntegrationTestSpec.assertEqualTextFiles(output, expected);
    }

    @Test
    public void testGVCFMode() throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testGVCFMode", ".g.vcf");
        final File expected = new File(TEST_FILES_DIR + "expected.testGVCFMode.g.vcf");

        final String[] args = new String[] {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-O", output.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.VERBOSITY_NAME, Log.LogLevel.INFO.toString(),
                "-ERC", "GVCF"
        };

        runCommandLine(args);
        IntegrationTestSpec.assertEqualTextFiles(output, expected);
    }
}
