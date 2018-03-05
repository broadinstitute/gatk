package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.apache.commons.io.FileUtils;
import org.testng.annotations.Test;
import java.io.IOException;
import java.io.File;

/**
 * Created by sam on 2/7/18.
 */
public class CNNVariantTrainIntegrationTest extends CommandLineProgramTest {

    @Test(groups = {"python"}, dependsOnGroups = {"writeTensors"})
    public void testTrainingReferenceModel() throws IOException{
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument("input-data-dir", CNNVariantWriteTensorsIntegrationTest.referenceTensorDir)
                .addArgument("tensor-name", TensorMapEnum.reference.name())
                .addArgument("epochs", "1")
                .addArgument("training-steps", "30")
                .addArgument("model-name", "test_reference_model")
                .addArgument("output-dir", CNNVariantWriteTensorsIntegrationTest.referenceTensorDir);

        runCommandLine(argsBuilder);
        FileUtils.deleteDirectory(new File(CNNVariantWriteTensorsIntegrationTest.referenceTensorDir));
    }

    @Test(groups = {"python"}, dependsOnGroups = {"writeTensors"})
    public void testTrainingReadModel() throws IOException{
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument("input-data-dir", CNNVariantWriteTensorsIntegrationTest.readTensorDir)
                .addArgument("tensor-name", TensorMapEnum.read_tensor.name())
                .addArgument("epochs", "1")
                .addArgument("training-steps", "5")
                .addArgument("validation-steps", "2")
                .addArgument("model-name", "test_read_tensor_model")
                .addArgument("output-dir", CNNVariantWriteTensorsIntegrationTest.readTensorDir)
                .addArgument("channels-last", "true");

        runCommandLine(argsBuilder);
        FileUtils.deleteDirectory(new File(CNNVariantWriteTensorsIntegrationTest.readTensorDir));
    }

}
