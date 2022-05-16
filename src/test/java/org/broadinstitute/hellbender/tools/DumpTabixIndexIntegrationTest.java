package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.sv.SiteDepthtoBAF;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Collections;

public class DumpTabixIndexIntegrationTest extends CommandLineProgramTest {
    @Test
    public void trivialDumpTabixTest() throws IOException {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.INPUT_SHORT_NAME, toolsTestDir + "test.tbi");
        argsBuilder.add(StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "%s");
        final IntegrationTestSpec testSpec =
                new IntegrationTestSpec(argsBuilder.getString(),
                        Collections.singletonList(toolsTestDir + "test.tbi.expected.txt"));
        testSpec.executeTest("dump tabix index", this);
    }
}
