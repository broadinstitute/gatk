package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.TestResources;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public final class ExampleReadWalkerWithReferenceSparkIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_DATA_DIRECTORY = TestResources.publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final String TEST_OUTPUT_DIRECTORY = TestResources.publicTestDir + "org/broadinstitute/hellbender/tools/examples/";

    @Test
    public void testExampleIntervalWalker() throws IOException {
        final File out = File.createTempFile("out", ".txt");
        out.delete();
        out.deleteOnExit();
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--input");
        args.add(TEST_DATA_DIRECTORY + "reads_data_source_test1.bam");
        args.add("--output");
        args.add(out.getAbsolutePath());
        args.add("--reference");
        args.add(TestResources.hg19MiniReference);
        this.runCommandLine(args.getArgsArray());
        File expected = new File(TEST_OUTPUT_DIRECTORY, "expected_ExampleReadWalkerWithReferenceIntegrationTest_output.txt");
        IntegrationTestSpec.assertEqualTextFiles(new File(out, "part-00000"), expected);
    }
}
