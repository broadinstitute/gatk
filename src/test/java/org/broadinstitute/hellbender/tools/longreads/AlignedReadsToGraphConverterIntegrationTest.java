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
        arguments.addBooleanArgument("create-dot-files", true);

        // Run the tool with our args:
        runCommandLine(arguments);
    }

    @Test
    public void testSerializeToGexfFile() {

        final File tmpDirName = createTempDir("alignedReadsGraphTest");

        final ArgumentsBuilder arguments = new ArgumentsBuilder();
        arguments.addInput(new File(aligned_bam_file));
        arguments.addArgument("output-file-base-name", tmpDirName.getAbsolutePath());
        arguments.addBooleanArgument("create-gexf-files", true);

        // Run the tool with our args:
        runCommandLine(arguments);
    }

    @Test
    public void testSerializeToGfaFile() {

        final File tmpDirName = createTempDir("alignedReadsGraphTest");

        final ArgumentsBuilder arguments = new ArgumentsBuilder();
        arguments.addInput(new File(aligned_bam_file));
        arguments.addArgument("output-file-base-name", tmpDirName.getAbsolutePath());

        // Run the tool with our args:
        runCommandLine(arguments);
    }

    @Test
    public void testSerialization() {

        final File tmpDirName = createTempDir("alignedReadsGraphTest");

        final String graphOutPath = tmpDirName.getAbsolutePath() + "/" + "graph_output.sdo";

        final ArgumentsBuilder graphCreationRunArgs = new ArgumentsBuilder();
        graphCreationRunArgs.addInput(new File(aligned_bam_file));
        graphCreationRunArgs.addArgument("output-file-base-name", tmpDirName.getAbsolutePath());
        graphCreationRunArgs.addArgument("graph-out", graphOutPath);
        graphCreationRunArgs.addBooleanArgument("skip-zip", true);

        // Run the tool with our args to create the kryo file:
        runCommandLine(graphCreationRunArgs);

        //==============================================================================

        final ArgumentsBuilder graphIngestRunArgs = new ArgumentsBuilder();
        graphIngestRunArgs.addInput(new File(aligned_bam_file));
        graphIngestRunArgs.addArgument("output-file-base-name", tmpDirName.getAbsolutePath());
        graphIngestRunArgs.addArgument("graph-in", graphOutPath);

        // Run the tool with our args to create the kryo file:
        runCommandLine(graphIngestRunArgs);
    }

}
