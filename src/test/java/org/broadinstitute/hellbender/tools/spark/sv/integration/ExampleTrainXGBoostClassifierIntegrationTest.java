package org.broadinstitute.hellbender.tools.spark.sv.integration;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.spark.sv.utils.ExampleTrainXGBoostClassifier;
import org.broadinstitute.hellbender.tools.spark.sv.utils.MachineLearningUtils;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

public class ExampleTrainXGBoostClassifierIntegrationTest extends CommandLineProgramTest {
    // TOOL_NAME = ExampleTrainXGBoostClassifier", but this allows IntelliJ to find usages
    private static final String TOOL_NAME = ExampleTrainXGBoostClassifier.class.getSimpleName();
    private static final String SV_UTILS_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/spark/sv/utils/";
    private static final String INPUT_FILE_PATH = SV_UTILS_DIR + "agaricus-integers.csv.gz";
    private static final MachineLearningUtils.TruthSet TRUTH_SET = MachineLearningUtils.loadCsvFile(INPUT_FILE_PATH);
    private static final double HYPERPARAMETER_TUNING_PROPORTION = 0.05;

    @Test(groups = "sv")
    protected void testTrainAndSaveClassifier() throws IOException {
        // require perfect accuracy of final classifier (this is easy, it's trained on whole data set)
        final File outputFile = runTool(INPUT_FILE_PATH);
        final MachineLearningUtils.GATKClassifier classifier = MachineLearningUtils.GATKClassifier.load(outputFile.getAbsolutePath());
        final int[] predictedLabels = classifier.predictClassLabels(TRUTH_SET.features);
        assertArrayEquals(predictedLabels, TRUTH_SET.classLabels, "predicted class labels should exactly match actual labels.");
    }

    private File runTool(final String inputPath) {
        final File outputFile = BaseTest.createTempFile("sv_xgboost_example_model", ".bin");
        outputFile.deleteOnExit();

        final ArgumentsBuilder argumentsBuilder = new ArgumentsBuilder();
        argumentsBuilder.addInput(new File(inputPath));
        argumentsBuilder.addOutput(outputFile);
        argumentsBuilder.addArgument(
                "hyperparameter-tuning-proportion", String.valueOf(HYPERPARAMETER_TUNING_PROPORTION)
        );
        argumentsBuilder.addBooleanArgument("auto-select-nthread", true);

        runCommandLine(Arrays.asList(argumentsBuilder.getArgsArray()), TOOL_NAME);
        return outputFile;
    }

    private static void assertArrayEquals(final int[] actual, final int[] expected, final String message) {
        Assert.assertEquals(expected.length, actual.length, message + ": array lengths not equal");
        for(int i = 0; i < actual.length; ++i) {
            Assert.assertEquals(actual[i], expected[i], message);
        }
    }
}
