package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Collections;

public final class ExampleVariantWalkerSparkIntegrationTest extends CommandLineProgramTest {

    private static final String TEST_DATA_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final String TEST_OUTPUT_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/tools/examples/";

    @Test
    public void testExampleVariantWalker() throws IOException {
        final File out = File.createTempFile("out", ".txt");
        out.delete();
        out.deleteOnExit();
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("-L 1:100-200");
        // note that joining with reads is not currently supported
        args.add("-V");
        args.add(TEST_DATA_DIRECTORY + "example_variants_withSequenceDict.vcf");
        args.add("-auxiliaryVariants");
        args.add(TEST_DATA_DIRECTORY + "feature_data_source_test.vcf");
        args.add("--output");
        args.add(out.getAbsolutePath());
        args.add("--reference");
        args.add(hg19MiniReference);
        this.runCommandLine(args.getArgsArray());
        File expected = new File(TEST_OUTPUT_DIRECTORY, "expected_ExampleVariantWalkerSparkIntegrationTest_output.txt");
        IntegrationTestSpec.assertEqualTextFiles(new File(out, "part-00000"), expected);
    }
    
    @Test
    public void testLongVariants() throws Exception {
        final File testOutput = createTempFile("output", ".txt");
        testOutput.delete();
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                        " -V " + TEST_OUTPUT_DIRECTORY + "input_long_variants.vcf" +
                        " -O " + testOutput, Collections.emptyList()
        );
        testSpec.executeTest("testExampleAssemblyRegionWalker", this);
        Assert.assertTrue(testOutput.exists());
    }
}
