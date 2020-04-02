package org.broadinstitute.hellbender.utils.python;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Collections;

/**
 * Tests to nominally validate that the GATK conda environment is activated and that its dependencies are accessible.
 */
public class PythonEnvironmentIntegrationTest {
    private final static String NL = System.lineSeparator();

    @DataProvider(name="supportedPythonPackages")
    public Object[][] getSupportedPythonPackages() {
        return new Object[][] {
                // names of base packages that we should be able to import from within the GATK conda environment
                // NOTE: these must be kept in sync with the versions in gatkcondaenv.yml.template
                { "mkl",            "2.3.0" },
                { "numpy",          "1.17.5" },
                { "theano",         "1.0.4" },
                { "tensorflow",     "1.15.0" },
                { "scipy",          "1.0.0" },
                { "pymc3",          "3.1" },
                { "keras",          "2.2.4" },
                { "sklearn",        "0.22.2.post1" },
                { "matplotlib",     "3.2.1" },
                { "pandas",         "1.0.3" },
                { "argparse",       null },
                { "gcnvkernel",     null }
        };
    }

    @Test(groups = {"python"}, dataProvider="supportedPythonPackages")
    public void testGATKPythonEnvironmentPackagePresent(final String packageName, final String expectedVersion) {
        // Sanity check to ensure that we can import these packages, and that conda is resolving them to
        // the specific version we requested in the conda environment definition.
        final StreamingPythonScriptExecutor<String> streamingPythonExecutor =
                new StreamingPythonScriptExecutor<>(true);
        Assert.assertNotNull(streamingPythonExecutor);
        Assert.assertTrue(streamingPythonExecutor.start(Collections.emptyList(), true, null));

        final String failedDependencyMessage = String.format(
                "The installed version of %s does not match the %s version that was requested. " +
                        "Check the build log to see the actual version that was resolved by conda.",
                packageName,
                expectedVersion);
        try {
            streamingPythonExecutor.sendSynchronousCommand(String.format("import %s" + NL, packageName));
        } catch (final PythonScriptExecutorException e) {
            throw new RuntimeException(failedDependencyMessage, e);
        }

        // for some dependencies, we also validate that the conda environment got the version that we asked for
        if (expectedVersion != null) {
            try {
                streamingPythonExecutor.sendSynchronousCommand(
                    String.format("assert(%s.__version__ == '%s')" + NL,
                            packageName,
                            expectedVersion));
            } catch (final PythonScriptExecutorException e) {
                throw new RuntimeException(failedDependencyMessage, e);
            }
        }
    }
}
