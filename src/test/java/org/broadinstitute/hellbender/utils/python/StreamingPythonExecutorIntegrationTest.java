package org.broadinstitute.hellbender.utils.python;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.Test;

import java.util.Collections;

public class StreamingPythonExecutorIntegrationTest extends GATKBaseTest {

    // This test validates that the StreamingPythonScriptExecutor will throw if the Python environment
    // is not activated, and will not pass if the environment is activated. It has to be an integration
    // test because unit tests are run on the Docker image, which has the Python environment activated.

    @Test()
    public void testRequirePythonEnvironment() {
        // This test is deliberately left out of the "python" test group in order to ensure that
        // it only executes when the Python environment has *NOT* been properly established. Also,
        // skip this test if we're running on the Docker because the Python environment is always
        // activated there.
        if (isGATKDockerContainer()) {
            throw new SkipException("Python environment validation test must be skipped when running on the Docker");
        }

        // validate that we throw if the GATK Python environment is not active
        final RuntimeException rte = Assert.expectThrows(RuntimeException.class, ()-> {
            final StreamingPythonScriptExecutor<String> streamingPythonExecutor =
                    new StreamingPythonScriptExecutor<>(PythonScriptExecutor.PythonExecutableName.PYTHON3, true);
            streamingPythonExecutor.start(Collections.emptyList());
        });

        // make sure that the underlying cause is actually a PythonScriptExecutorException
        Assert.assertEquals(rte.getCause().getClass(), PythonScriptExecutorException.class);
    }
}
