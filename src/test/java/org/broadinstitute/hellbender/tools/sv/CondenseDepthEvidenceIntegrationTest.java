package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Collections;

public class CondenseDepthEvidenceIntegrationTest  extends CommandLineProgramTest {
    public static final String printEvidenceTestDir = toolsTestDir + "walkers/sv/printevidence/";
    public static final String sequenceDict = largeFileTestDir + "human_g1k_v37.20.21.dict";

    @Test
    public void positiveTest() throws IOException {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VERBOSITY_NAME, Log.LogLevel.ERROR.name());
        argsBuilder.add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, sequenceDict);
        argsBuilder.add(CondenseDepthEvidence.MIN_INTERVAL_SIZE_ARGNAME, 101);
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_SHORT_NAME, printEvidenceTestDir + "output.test_hg38.rd.txt");
        argsBuilder.add(StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "%s");
        final IntegrationTestSpec testSpec =
                new IntegrationTestSpec(argsBuilder.getString(),
                        Collections.singletonList(printEvidenceTestDir + "merged.rd.txt"));
        testSpec.setOutputFileExtension("rd.txt");
        testSpec.executeTest("positive test", this);
    }

    @Test
    public void negativeTest() throws IOException {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VERBOSITY_NAME, Log.LogLevel.ERROR.name());
        argsBuilder.add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, sequenceDict);
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_SHORT_NAME, printEvidenceTestDir + "no_condense.rd.txt");
        argsBuilder.add(StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "%s");

        // output equal to input, because this file contains rows that can't be merged for various reasons
        final IntegrationTestSpec testSpec =
                new IntegrationTestSpec(argsBuilder.getString(),
                        Collections.singletonList(printEvidenceTestDir + "no_condense.rd.txt"));
        testSpec.setOutputFileExtension("rd.txt");
        testSpec.executeTest("negative test", this);
    }
}
