package org.broadinstitute.hellbender.utils.python;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class PythonScriptExecutorUnitTest extends GATKBaseTest {
    final String HELLO_WORLD_SCRIPT = "print (\"hello, world\")";
    final String EXIT_WITH_CODE_SCRIPT = "import sys;sys.exit(239)";

    @Test(groups = "python")
    public void testPythonExists() {
        Assert.assertTrue(
                new PythonScriptExecutor(true).externalExecutableExists(),
                "python not found in environment ${PATH}"
        );
    }

    @Test(groups = "python", dependsOnMethods = "testPythonExists")
    public void testExecuteAsCommand() {
        final PythonScriptExecutor pythonExecutor = new PythonScriptExecutor(true);
        final boolean ret = pythonExecutor.executeCommand(HELLO_WORLD_SCRIPT, null, null);

        Assert.assertTrue(ret, "Python exec failed");
    }

    @Test(groups = "python", dependsOnMethods = "testPythonExists")
    public void testExecuteAsModule() {
        // use builtin module "random"
        final PythonScriptExecutor executor = new PythonScriptExecutor(true);
        final boolean ret = executor.executeModule("random", null, null);

        Assert.assertTrue(ret, "Python exec failed");
    }

    @Test(groups = "python", dependsOnMethods = "testPythonExists")
    public void testExecuteAsScript() {
        final File scriptFile = writeTemporaryScriptFile(HELLO_WORLD_SCRIPT, PythonScriptExecutor.PYTHON_EXTENSION);

        final PythonScriptExecutor executor = new PythonScriptExecutor(true);
        final boolean ret = executor.executeScript(scriptFile.getAbsolutePath(), null, null);

        Assert.assertTrue(ret, "Python exec failed");
    }

    @Test(groups = "python", dependsOnMethods = "testPythonExists")
    public void testExecuteAsScriptAndGetProcessOutput() {
        final File scriptFile = writeTemporaryScriptFile(HELLO_WORLD_SCRIPT, PythonScriptExecutor.PYTHON_EXTENSION);

        final PythonScriptExecutor executor = new PythonScriptExecutor(true);
        final ProcessOutput po = executor.executeScriptAndGetOutput(scriptFile.getAbsolutePath(), null, null);

        Assert.assertEquals(po.getExitValue(), 0);
    }

    @Test(groups = "python", dependsOnMethods = "testPythonExists")
    public void testExecuteAsScriptAndGetProcessOutputWithFailureCode() {
        final File scriptFile = writeTemporaryScriptFile(EXIT_WITH_CODE_SCRIPT, PythonScriptExecutor.PYTHON_EXTENSION);

        final PythonScriptExecutor executor = new PythonScriptExecutor(true);
        final ProcessOutput po = executor.executeScriptAndGetOutput(scriptFile.getAbsolutePath(), null, null);

        Assert.assertEquals(po.getExitValue(), 239);
    }

    @Test(groups = "python", dependsOnMethods = "testPythonExists")
    public void testExecuteAsRawArgs() {
        final PythonScriptExecutor pythonExecutor = new PythonScriptExecutor(true);
        final boolean ret = pythonExecutor.executeArgs(new ArrayList<String>(Arrays.asList("-c", HELLO_WORLD_SCRIPT)));

        Assert.assertTrue(ret, "Python exec failed");
    }

    @Test(groups = "python", dependsOnMethods = "testPythonExists")
    public void testExecuteAsRawArgsAndGetProcessOutput() {
        final PythonScriptExecutor pythonExecutor = new PythonScriptExecutor(true);
        final ProcessOutput po = pythonExecutor.executeArgsAndGetOutput(new ArrayList<String>(Arrays.asList("-c", HELLO_WORLD_SCRIPT)));

        Assert.assertEquals(po.getExitValue(), 0);
    }

    @Test(groups = "python", dependsOnMethods = "testPythonExists")
    public void testExecuteAsRawArgsAndGetProcessOutputWithFailureCode() {
        final PythonScriptExecutor pythonExecutor = new PythonScriptExecutor(true);
        final ProcessOutput po = pythonExecutor.executeArgsAndGetOutput(new ArrayList<String>(Arrays.asList("-c", EXIT_WITH_CODE_SCRIPT)));

        Assert.assertEquals(po.getExitValue(), 239);
    }

    @Test(groups = "python", dependsOnMethods = "testPythonExists")
    public void testExecuteAsRawArgsSerial() {
        final PythonScriptExecutor pythonExecutor = new PythonScriptExecutor(true);
        final List<String> args = new ArrayList<>(Arrays.asList("-c", HELLO_WORLD_SCRIPT));

        Assert.assertTrue(pythonExecutor.executeArgs(args), "First Python exec failed");
        Assert.assertTrue(pythonExecutor.executeArgs(args), "Second Python exec failed");
    }

    @Test(groups = "python", dependsOnMethods = "testPythonExists")
    public void testExecuteAsRawArgsWithPythonArguments() {
        final PythonScriptExecutor pythonExecutor = new PythonScriptExecutor(true);
        final boolean ret = pythonExecutor.executeArgs(new ArrayList<>(Arrays.asList("-v", "-c", HELLO_WORLD_SCRIPT )));

        Assert.assertTrue(ret, "Python exec failed");
    }

    @Test(groups = "python", dependsOnMethods = "testPythonExists")
    public void testExecuteAsScriptWithScriptArguments() {
        final String SCRIPT_WITH_ARGUMENTS =
                "import fileinput\n" +
                "for line in fileinput.input():\n" +
                "    print(line)\n";
        final File scriptFile = writeTemporaryScriptFile(SCRIPT_WITH_ARGUMENTS, PythonScriptExecutor.PYTHON_EXTENSION);
        final PythonScriptExecutor executor = new PythonScriptExecutor(true);
        final boolean ret = executor.executeScript(
                scriptFile.getAbsolutePath(),
                null,
                Collections.singletonList(scriptFile.getAbsolutePath()));

       Assert.assertTrue(ret, "Python exec failed");
    }

    @Test(groups = "python", dependsOnMethods = "testPythonExists", expectedExceptions = PythonScriptExecutorException.class)
    public void testNonExistentScriptException() throws IOException {
        final PythonScriptExecutor executor = new PythonScriptExecutor(true);
        executor.executeScript(
                BaseTest.getSafeNonExistentFile("Nonexistent" + PythonScriptExecutor.PYTHON_EXTENSION).getAbsolutePath(),
                null,
                null);
    }

    @Test(groups = "python", dependsOnMethods = "testPythonExists")
    public void testNonExistentScriptNoException() {
        final PythonScriptExecutor executor = new PythonScriptExecutor(true);
        executor.setIgnoreExceptions(true);
        Assert.assertFalse(
                executor.executeScript(
                        BaseTest.getSafeNonExistentFile("Nonexistent" + PythonScriptExecutor.PYTHON_EXTENSION).getAbsolutePath(),
                        null,
                        null
                ),
                "Exec should have returned false when the job failed");
    }

    private File writeTemporaryScriptFile(final String content, final String suffix) {
        final File tempScriptFile = IOUtils.writeTempFile(content, "myTestScript", suffix);
        tempScriptFile.deleteOnExit();
        return tempScriptFile;
    }

}
