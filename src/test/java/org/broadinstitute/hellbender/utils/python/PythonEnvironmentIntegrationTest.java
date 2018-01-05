package org.broadinstitute.hellbender.utils.python;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

// Tests to nominally validate that the GATK conda environment is activated, and that it's
// dependencies are accessible.
//
// We could validate the version numbers of the dependencies, but that would require
// changes every time we change the gatkcondaenv.yml file.
//
public class PythonEnvironmentIntegrationTest {
    private final static String NL = System.lineSeparator();

    @DataProvider(name="supportedPythonPackages")
    public Object[][] getSupportedPythonPackages() {
        return new Object[][] {
                // names of base packages that we should be able to import from within the GATK conda environment
                // this list isn't exhaustive
                { "numpy"},
                { "scipy"},
                { "tensorflow"},
                { "theano" },
                { "keras" },
                { "pymc3" },
                { "argparse" },
                { "gcnvkernel" }
        };
    }

    @Test(groups = {"python"}, dataProvider="supportedPythonPackages")
    public void testGATKPythonEnvironmentPackagePresent(final String packageName) {
        // Do a basic sanity check of the GATK Python conda environment. This test should only be run on
        // the GATK docker image, or if the conda environment has been activated manually.

        // We use the default python executable name ("python"), which in the activated gatk conda env should be Python 3.6.1
        final PythonScriptExecutor pythonExecutor = new PythonScriptExecutor(true);
        Assert.assertTrue(pythonExecutor.executeCommand(String.format("import %s", packageName) + NL,null,null));
    }

}
