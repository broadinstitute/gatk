package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public final class ExampleReadWalkerWithVariantsSparkIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_DATA_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final String TEST_OUTPUT_DIRECTORY = exampleTestDir;

    @Test
    public void testExampleReadWalkerWithVariantsSpark() throws IOException {
        final File out = IOUtils.createTempFile("out", ".txt");
        out.delete();
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("--input");
        args.addRaw(TEST_DATA_DIRECTORY + "reads_data_source_test1.bam");
        args.addRaw("--output");
        args.addRaw(out.getAbsolutePath());
        args.addRaw("--reference");
        args.addRaw(hg19MiniReference);
        args.addRaw("--variant ");
        args.addRaw(TEST_DATA_DIRECTORY + "example_variants_withSequenceDict.vcf");
        this.runCommandLine(args.getArgsArray());
        File expected = new File(TEST_OUTPUT_DIRECTORY, "expected_ExampleReadWalkerWithVariantsIntegrationTest_output.txt");
        IntegrationTestSpec.assertEqualTextFiles(new File(out, "part-00000"), expected);
    }
}
