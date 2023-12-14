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

    @DataProvider(name="dataPackagePresent")
    public Object[][] getDataPackagePresent() {
        return new Object[][] {
                // names of base packages that we should be able to import from within the GATK conda environment
                // NOTE: these must be kept in sync with the versions in gatkcondaenv.yml.template
                { "mkl",            "2.4.0" },
                { "numpy",          "1.26.2" },
                { "pytensor",       "2.18.1" },
                { "torch",          "2.1.0.post100" },
                { "scipy",          "1.11.4" },
                { "pymc",           "5.10.0" },
                { "h5py",           "3.10.0" },
                { "sklearn",        "1.3.2" },
                { "matplotlib",     "3.8.2" },
                { "pandas",         "2.1.3" },
                { "argparse",       null },
                { "gcnvkernel",     null },

                // R packages
                // Commented out since we can't check versions of R packages in this test
                // { "r-backports",    "1.1.10" },
        };
    }

    /**
     * See documentation in scripts/gatkcondaenv.yml.template.
     */
    @DataProvider(name="dataPackageMKLEnabled")
    public Object[][] getDataPackageMKLEnabled() {
        return new Object[][] {
                // { "numpy",          "numpy.__config__.get_info('blas_mkl_info') != {} and numpy.__config__.get_info('lapack_mkl_info') != {}" },  TODO this check might need to be done differently for conda-forge numpy 1.26.2, we remove it for now
                { "pytensor",       "'-lmkl_rt' in pytensor.config.blas__ldflags" },
                { "torch",          "'BLAS_INFO=mkl' in torch.__config__.show() and 'USE_MKL=ON' in torch.__config__.show()" }
        };
    }

    @Test(groups = {"python"}, dataProvider="dataPackagePresent")
    public void testGATKPythonEnvironmentPackagePresent(final String packageName,
                                                        final String expectedVersion) {
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
            Assert.fail(failedDependencyMessage);
        }

        // for some dependencies, we also validate that the conda environment got the version that we asked for
        if (expectedVersion != null) {
            try {
                streamingPythonExecutor.sendSynchronousCommand(
                        String.format("assert %s.__version__ == '%s'" + NL,
                                packageName,
                                expectedVersion));
            } catch (final PythonScriptExecutorException e) {
                Assert.fail(failedDependencyMessage);
            }
        }
    }

    /**
     * See documentation in scripts/gatkcondaenv.yml.template.
     */
    @Test(groups = {"python"}, dataProvider="dataPackageMKLEnabled")
    public void testGATKPythonEnvironmentPackageMKLEnabled(final String packageName,
                                                           final String assertMKLEnabled) {
        final StreamingPythonScriptExecutor<String> streamingPythonExecutor =
                new StreamingPythonScriptExecutor<>(true);
        Assert.assertNotNull(streamingPythonExecutor);
        Assert.assertTrue(streamingPythonExecutor.start(Collections.emptyList(), true, null));

        final String failedAssertMKLEnabledMessage = String.format("The installed version of %s does not have MKL enabled.", packageName);
        try {
            streamingPythonExecutor.sendSynchronousCommand(String.format("import %s; assert(%s)" + NL, packageName, assertMKLEnabled));
        } catch (final PythonScriptExecutorException e) {
            Assert.fail(failedAssertMKLEnabledMessage);
        }
    }
}
