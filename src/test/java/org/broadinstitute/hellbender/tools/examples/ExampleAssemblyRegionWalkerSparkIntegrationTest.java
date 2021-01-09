package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class ExampleAssemblyRegionWalkerSparkIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_OUTPUT_DIRECTORY = exampleTestDir;

    @Test()
    public void testExampleAssemblyRegionWalkerStrict() throws IOException {
        final File out = IOUtils.createTempFile("out", ".txt");
        out.delete();
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("--input");
        args.addRaw(NA12878_20_21_WGS_bam);
        args.addRaw("--output");
        args.addRaw(out.getAbsolutePath());
        args.addRaw("--reference");
        args.addRaw(b37_reference_20_21);
        args.addRaw("-knownVariants " + dbsnp_138_b37_20_21_vcf);
        args.addRaw("-L 20:10000000-10050000");
        args.addRaw("--strict"); // needed to match walker version
        this.runCommandLine(args.getArgsArray());
        File expected = new File(TEST_OUTPUT_DIRECTORY, "expected_ExampleAssemblyRegionWalkerIntegrationTest_output.txt");
        IntegrationTestSpec.assertEqualTextFiles(new File(out, "part-00000"), expected);
    }

    @Test()
    public void testExampleAssemblyRegionWalkerNonStrict() throws IOException {
        final File out = IOUtils.createTempFile("out", ".txt");
        out.delete();
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("--input");
        args.addRaw(NA12878_20_21_WGS_bam);
        args.addRaw("--output");
        args.addRaw(out.getAbsolutePath());
        args.addRaw("--reference");
        args.addRaw(b37_reference_20_21);
        args.addRaw("-knownVariants " + dbsnp_138_b37_20_21_vcf);
        args.addRaw("-L 20:10000000-10050000");
        this.runCommandLine(args.getArgsArray());
        File expected = new File(TEST_OUTPUT_DIRECTORY, "expected_ExampleAssemblyRegionWalkerSparkIntegrationTest_non_strict_output.txt");
        IntegrationTestSpec.assertEqualTextFiles(new File(out, "part-00000"), expected);
    }
}
