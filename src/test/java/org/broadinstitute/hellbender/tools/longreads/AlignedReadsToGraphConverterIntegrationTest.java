package org.broadinstitute.hellbender.tools.longreads;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.annotations.Test;

import java.io.File;

public class AlignedReadsToGraphConverterIntegrationTest extends CommandLineProgramTest {

    private static final String aligned_bam_file = toolsTestDir + "longreads" + File.separator + "NA19240.chr20_19294351-19304184.corrected.bam";

    @Test
    public void testSerializeToDotFile() {

        final File tmpDirName = createTempDir("alignedReadsGraphTest");

        final ArgumentsBuilder arguments = new ArgumentsBuilder();
        arguments.addInput(new File(aligned_bam_file));
        arguments.addArgument("output-file-base-name", tmpDirName.getAbsolutePath());

        // Run the tool with our args:
        runCommandLine(arguments);
    }

}
