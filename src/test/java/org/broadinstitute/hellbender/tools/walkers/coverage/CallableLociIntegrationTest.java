package org.broadinstitute.hellbender.tools.walkers.coverage;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class CallableLociIntegrationTest extends CommandLineProgramTest {

    public static final String testDataDir = largeFileTestDir + "/org/broadinstitute/hellbender/tools/walkers/coverage/callableloci/";

    @Test
    public void testBasicOperation() throws Exception {
        final File outputFile = createTempFile("callableLoci", ".bed");
        final File summaryFile = createTempFile("callableLoci", ".summary.txt");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(new File(largeFileTestDir + "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam"))
            .addReference(new File(b37_reference_20_21))
            .addOutput(outputFile)
            .add("summary", summaryFile)
            .add("format", "BED");

        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(
            outputFile,
            new File(testDataDir + "callable_loci.testBasicOperation.expected.bed"),
            "#");

        IntegrationTestSpec.assertEqualTextFiles(
            summaryFile,
            new File(testDataDir + "callable_loci.testBasicOperation.expected.summary.txt"),
            "#");
    }
}