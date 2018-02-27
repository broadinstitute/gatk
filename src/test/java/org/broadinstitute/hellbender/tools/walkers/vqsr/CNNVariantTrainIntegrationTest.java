package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

/**
 * Created by sam on 2/7/18.
 */
public class CNNVariantTrainIntegrationTest extends CommandLineProgramTest {


    @Test(groups = {"python"})
    public void testTrainingReferenceModel() throws IOException{
        CNNVariantWriteTensorsIntegrationTest wt = new CNNVariantWriteTensorsIntegrationTest();
        wt.generateReferenceTensors();
        final String dataDir = largeFileTestDir + "VQSR/reference_tensors/";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument("input-data-dir", dataDir)
                .addArgument("tensor-name", TensorMapEnum.reference.name())
                .addArgument("epochs", "1")
                .addArgument("training-steps", "50")
                .addArgument("model-name", "test_reference_model")
                .addArgument("output-dir", dataDir);

        runCommandLine(argsBuilder);
    }

    @Test(groups = {"python"})
    public void testTrainingReadModel() throws IOException{
        CNNVariantWriteTensorsIntegrationTest wt = new CNNVariantWriteTensorsIntegrationTest();
        wt.generateReadTensors();
        final String dataDir = largeFileTestDir + "VQSR/read_tensors/";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument("input-data-dir", dataDir)
                .addArgument("tensor-name", TensorMapEnum.read_tensor.name())
                .addArgument("epochs", "1")
                .addArgument("training-steps", "1")
                .addArgument("validation-steps", "1")
                .addArgument("model-name", "test_read_tensor_model")
                .addArgument("output-dir", dataDir)
                .addArgument("channels-last", "false");

        runCommandLine(argsBuilder);
    }

}
