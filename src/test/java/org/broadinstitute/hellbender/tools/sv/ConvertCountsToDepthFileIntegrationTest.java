package org.broadinstitute.hellbender.tools.sv;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.codecs.DepthEvidenceCodec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Collections;

public class ConvertCountsToDepthFileIntegrationTest extends CommandLineProgramTest {
    public static final String converterTestDir = toolsTestDir + "walkers/sv/convertCountsToDepthFile";

    @Test
    public void testSuccessfullConvert() throws IOException {
        final String inputFilename = converterTestDir + "/NA12878.counts.tsv";
        final String expectedOutput = converterTestDir + "/expected_output" + DepthEvidenceCodec.FORMAT_SUFFIX;

        final String args =
                " --counts-filename " + inputFilename +
                " --" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME + " " + FULL_HG38_DICT +
                " --sample-name test" +
                " -" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME + " %s";

        final IntegrationTestSpec testSpec = new IntegrationTestSpec(args, Collections.singletonList(expectedOutput));
        testSpec.setOutputFileExtension(DepthEvidenceCodec.FORMAT_SUFFIX);

        testSpec.executeTest("test successful convert counts to depth file", this);
    }
}
