package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.hellbender.utils.io.IOUtils;
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
        final Resource cnnSerializedResource = new Resource("cnn_1d_annotations.hd5", NeuralNetExecutor.class);
        final String inputFile = largeFileTestDir + "VQSR/phase1.projectConsensus.chr20.1M-10M.raw.snps.vcf";
        final File tempResourceFile = IOUtils.writeTempResource(cnnSerializedResource);

        final boolean pythonReturnCode = pythonExecutor.executeScript(
                pythonScriptResource,
                null,
                Arrays.asList(
                        "--architecture",
                        tempResourceFile.getAbsolutePath(),
                        "--input_vcf",
                        inputFile
                )
        );

        Assert.assertTrue(pythonReturnCode, "Python exec failed");
    }
}
