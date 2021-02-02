package org.broadinstitute.hellbender.tools.dragstr;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class CalibrateDragstrModelIntegrationTest extends CommandLineProgramTest {

    private static final File HUMAN_B37_20_21_TBL_GZ =
            new File("src/test/resources/large/org/broadinstitute/hellbender/tools/dragstr/compose-str-table-file-human_g1k_v37_20.21.dragen.zip");


    // Asserting that ReadsDataSourcePool passes the reference through to the underlying reader
    @Test()
    public void testReferenceGetsProvidedToCramInputWhenThreaded() throws IOException {
        final File output = File.createTempFile("test-output", ".zip");
        output.deleteOnExit();

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add(StandardArgumentDefinitions.REFERENCE_SHORT_NAME, b37_reference_20_21);
        args.add(CalibrateDragstrModel.STR_TABLE_PATH_FULL_NAME, HUMAN_B37_20_21_TBL_GZ);
        args.addFlag(CalibrateDragstrModel.PARALLEL_FULL_NAME);
        args.add(CalibrateDragstrModel.THREADS_FULL_NAME, 10);
        args.add(CalibrateDragstrModel.SHARD_SIZE_FULL_NAME,  CalibrateDragstrModel.DEFAULT_SHARD_SIZE);
        args.add(CalibrateDragstrModel.DOWN_SAMPLE_SIZE_FULL_NAME, CalibrateDragstrModel.DEFAULT_DOWN_SAMPLE_SIZE);
        args.addInput(NA12878_20_21_WGS_cram);
        args.addOutput(output);

        runCommandLine(args);

    }
}