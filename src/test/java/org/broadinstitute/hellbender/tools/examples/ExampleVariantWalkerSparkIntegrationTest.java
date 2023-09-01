package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public final class ExampleVariantWalkerSparkIntegrationTest extends CommandLineProgramTest {

    private static final String TEST_DATA_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final String TEST_OUTPUT_DIRECTORY = exampleTestDir;

    @Test
    public void testExampleVariantWalker() throws IOException {
        final File out = IOUtils.createTempFile("out", ".txt");
        Assert.assertTrue(out.delete(), "failed to perform necessary deletion during test setup");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("-L 1:100-200");
        // note that joining with reads is not currently supported
        args.addRaw("-V");
        args.addRaw(TEST_DATA_DIRECTORY + "example_variants_withSequenceDict.vcf");
        args.addRaw("-auxiliaryVariants");
        args.addRaw(TEST_DATA_DIRECTORY + "feature_data_source_test.vcf");
        args.addRaw("--output");
        args.addRaw(out.getAbsolutePath());
        args.addRaw("--reference");
        args.addRaw(hg19MiniReference);
        this.runCommandLine(args.getArgsArray());
        File expected = new File(TEST_OUTPUT_DIRECTORY, "expected_ExampleVariantWalkerSparkIntegrationTest_output.txt");
        IntegrationTestSpec.assertEqualTextFiles(new File(out, "part-00000"), expected);
    }

}
