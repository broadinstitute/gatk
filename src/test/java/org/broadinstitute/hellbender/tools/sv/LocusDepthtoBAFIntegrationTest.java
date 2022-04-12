package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Collections;

public class LocusDepthtoBAFIntegrationTest extends CommandLineProgramTest {
    public static final String ld2bafTestDir = toolsTestDir + "walkers/sv/LocusDepthtoBAF/";

    @Override public String getTestedToolName() { return "LocusDepthtoBAF"; }

    @Test
    public void testStdDevTooBig() throws IOException {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VERBOSITY_NAME, Log.LogLevel.ERROR.name());
        argsBuilder.add(LocusDepthtoBAF.MIN_HET_LIKELIHOOD, "0.");
        argsBuilder.add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, hg19_chr1_1M_dict);
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_SHORT_NAME, ld2bafTestDir + "test1.ld.txt");
        argsBuilder.add(StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "%s");
        final IntegrationTestSpec testSpec =
                new IntegrationTestSpec(argsBuilder.getString(),
                        Collections.singletonList(ld2bafTestDir + "test1.baf.txt"));
        testSpec.setOutputFileExtension("baf.txt");
        testSpec.executeTest("test1: standard deviation too big", this);
    }

    @Test
    public void testAOK() throws IOException {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VERBOSITY_NAME, Log.LogLevel.ERROR.name());
        argsBuilder.add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, hg19_chr1_1M_dict);
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_SHORT_NAME, ld2bafTestDir + "test2.ld.txt");
        argsBuilder.add(StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "%s");
        final IntegrationTestSpec testSpec =
                new IntegrationTestSpec(argsBuilder.getString(),
                        Collections.singletonList(ld2bafTestDir + "test2.baf.txt"));
        testSpec.setOutputFileExtension("baf.txt");
        testSpec.executeTest("test2: a-ok", this);
    }

    @Test
    public void testAdjustMedian() throws IOException {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VERBOSITY_NAME, Log.LogLevel.ERROR.name());
        argsBuilder.add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, hg19_chr1_1M_dict);
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_SHORT_NAME, ld2bafTestDir + "test3.ld.txt");
        argsBuilder.add(StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "%s");
        final IntegrationTestSpec testSpec =
                new IntegrationTestSpec(argsBuilder.getString(),
                        Collections.singletonList(ld2bafTestDir + "test3.baf.txt"));
        testSpec.setOutputFileExtension("baf.txt");
        testSpec.executeTest("test3: adjust median", this);
    }
}
