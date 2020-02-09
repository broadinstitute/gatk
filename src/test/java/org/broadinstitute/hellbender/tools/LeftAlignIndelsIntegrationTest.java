package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class LeftAlignIndelsIntegrationTest extends CommandLineProgramTest {

    // test that a bam whose indels were manually unaligned is correctly restored to left alignment
    @Test
    public void testLeftAlignIndels() throws IOException {
        final File input = new File(getToolTestDataDir(), "non-aligned.sam");
        final File expected = new File(getToolTestDataDir(), "left-aligned.sam");
        final File output = createTempFile("output", ".sam");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addReference(b37Reference)
                .addInput(input)
                .addOutput(output);

        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expected, "@");
    }
}