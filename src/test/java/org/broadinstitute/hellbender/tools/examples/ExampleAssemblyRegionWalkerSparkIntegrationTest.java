package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;

public class ExampleAssemblyRegionWalkerSparkIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_OUTPUT_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/tools/examples/";

    @Test
    public void testExampleAssemblyRegionWalker() throws Exception {
        final File out = File.createTempFile("out", ".txt");
        out.delete();
        out.deleteOnExit();
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--input");
        args.add(NA12878_20_21_WGS_bam);
        args.add("--output");
        args.add(out.getAbsolutePath());
        args.add("--reference");
        args.add(b37_reference_20_21);
        args.add("-knownVariants " + dbsnp_138_b37_20_21_vcf);
        args.add("-L 20:10000000-10050000");
        this.runCommandLine(args.getArgsArray());
        File expected = new File(TEST_OUTPUT_DIRECTORY, "expected_ExampleAssemblyRegionWalkerIntegrationTest_output.txt");
        IntegrationTestSpec.assertEqualTextFiles(new File(out, "part-00000"), expected);
    }
}
