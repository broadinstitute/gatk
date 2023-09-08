package org.broadinstitute.hellbender.tools.walkers.groundtruth;

import org.apache.commons.math3.util.Precision;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.walkers.variantrecalling.FlowTestConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

public class AddFlowBaseQualityIntegrationTest extends CommandLineProgramTest {

    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;
    public static final String OUTPUT_FILENAME = "add_flow_base_quality_output.sam";

    private static String testDir = publicTestDir + FlowTestConstants.GROUND_TRUTH_DATA_DIR;

    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled");
    }

    @Test
    public void testBasic() throws IOException {

        final File outputDir = createTempDir("testGroundTruthTest");
        final File expectedFile = new File(testDir + "/" + OUTPUT_FILENAME);
        final File outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : new File(outputDir + "/" + OUTPUT_FILENAME);

        final String[] args = buildCommonArgs(outputFile);

        runCommandLine(args);  // no assert, just make sure we don't throw

        // make sure we've generated the otuput file
        Assert.assertTrue(outputFile.exists());

        // walk the output and expected files, compare non-comment lines
        if ( !UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedFile, "@");
        }
    }

    private String[] buildCommonArgs(final File outputFile) {

        return new String[] {
                "-R", largeFileTestDir + "/Homo_sapiens_assembly38.fasta.gz",
                "-I", testDir + "/gt_scorer_input.bam",
                "--output", outputFile.getAbsolutePath(),
                "--intervals", "chr9:71000-74000",
        };
    }

    @DataProvider(name = "generateHmerBaseErrorProbabilities")
    public Object[][] getGenerateHmerBaseErrorProbabilities() {
        return new Object[][] {
                {
                    new int[] { 1, 0, 0, 1, 0, 1, 1 }, // key
                        new double[][] {
                                { 0.1, 0.1, 0.01, 0.1, 0.001, 0.05, 0.1},       // 1 less
                                { 0.9, 0.9, 0.99, 0.9, 0.999, 0.95, 0.9},       // key
                                { 0.1, 0.1, 0.01, 0.1, 0.001, 0.05, 0.1}        // 1 more
                        }, // error prob bands
                        3, // flow
                        4, // flowOrderLength
                        new double[] {0.00112, 0}, // result
                        5 // result precision
                },
                {
                        new int[] { 1, 0, 0, 2, 0, 1, 1 }, // key
                        new double[][] {
                                { 0.1, 0.1, 0.01, 0.1, 0.001, 0.05, 0.1},       // 1 less
                                { 0.9, 0.9, 0.99, 0.9, 0.999, 0.95, 0.9},       // key
                                { 0.1, 0.1, 0.01, 0.1, 0.001, 0.05, 0.1}        // 1 more
                        }, // error prob bands
                        3, // flow
                        4, // flowOrderLength
                        new double[] {0.02516, 0.00592}, // result (errorProbs)
                        5 // result precision
                }
        };
    }

    @Test(dataProvider = "generateHmerBaseErrorProbabilities")
    public void testGenerateHmerBaseErrorProbabilities(final int[] key, final double[][] errorProbBands, final int flow, final int flowOrderLength, final double[] result, final int resultPrecision) {

        final double[]        errorProbs = AddFlowBaseQuality.generateHmerBaseErrorProbabilities(key, errorProbBands, flow, flowOrderLength);

        Assert.assertEquals(Arrays.stream(errorProbs).map(v -> Precision.round(v, resultPrecision)).toArray(), result);
    }
}
