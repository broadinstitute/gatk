package org.broadinstitute.hellbender.tools.sv;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.codecs.DepthEvidenceCodec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Collections;

public class ConvertCountsToDepthFileIntegrationTest extends CommandLineProgramTest {
    public static final String converterTestDir = toolsTestDir + "walkers/sv/convertCountsToDepthFile";

    public record TestCase(
            String testLabel,
            String inputFilename,
            String expectedOutputFilename,
            String sampleLabel)
    { }

    @DataProvider
    public TestCase[] cases()
    {
        return new TestCase[]
            {
                new TestCase(
                    "Test successful convert counts to depth file",
                    converterTestDir + "/NA12878.counts.tsv",
                    converterTestDir + "/expected_output" + DepthEvidenceCodec.FORMAT_SUFFIX + ".gz",
                    "NA12878"),
                new TestCase(
                    "Test successful read sample name from file header",
                    converterTestDir + "/NA12878.counts.tsv",
                    converterTestDir + "/expected_output" + DepthEvidenceCodec.FORMAT_SUFFIX + ".gz",
                    null)
            };
    }
    @Test(dataProvider="cases")
    public void testSuccessfulConvert(final TestCase testCase) throws IOException {
        final String args =
            " --" + ConvertCountsToDepthFile.COUNTS_FILE_PATH_ARG_FULL_NAME + " " + testCase.inputFilename +
            " --" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME + " " + FULL_HG38_DICT +
            ((testCase.inputFilename == null) ?
                    "" :
                    (" --" + ConvertCountsToDepthFile.SAMPLE_NAME_ARG_FULL_NAME + " " + testCase.sampleLabel)) +
            " -" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME + " %s";

        final IntegrationTestSpec testSpec = new IntegrationTestSpec(args, Collections.singletonList(testCase.expectedOutputFilename));
        testSpec.setOutputFileExtension(DepthEvidenceCodec.FORMAT_SUFFIX + ".gz");
        testSpec.executeTest(testCase.testLabel, this);
    }
}
