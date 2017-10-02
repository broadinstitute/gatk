package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;

/**
 * Created by sam on 9/25/17.
 */
public class ApplyNeuralNetUnitTest extends BaseTest {

    @Test(groups = {"PYTHON_VQSR"})
    public void testExecuteAsScriptWithScriptArguments() {
        final PythonScriptExecutor pythonExecutor = new PythonScriptExecutor(true);
        final Resource pythonScriptResource = new Resource("ApplyNeuralNetModel.py", NeuralNetExecutor.class);
        final boolean pythonReturnCode = pythonExecutor.executeScript(
                pythonScriptResource,
                null,
                Arrays.asList(
                        "--architecture",
                        "/dsde/working/sam/palantir_cnn/Analysis/vqsr_cnn/weights/m__base_quality_mode_phot__channels_last_False__id_g94982_no_qual_train2__window_size_128__read_limit_128__random_seed_12878__tensor_map_2d_mapping_quality__mode_ref_read_anno.hd5",
                        "--tensors",
                        "/Users/sam/vqsr_data/tensors/not_indel_subset/"
                )
        );

        Assert.assertTrue(pythonReturnCode, "Python exec failed");
    }
}
