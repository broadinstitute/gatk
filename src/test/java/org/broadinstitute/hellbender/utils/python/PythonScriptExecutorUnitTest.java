package org.broadinstitute.hellbender.utils.python;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class PythonScriptExecutorUnitTest extends BaseTest {
    final String HELLO_WORLD_SCRIPT = "print \"hello, world\"";

    @Test(groups = {"Python"})
    public void testPythonExists() {
        Assert.assertTrue(
                new PythonScriptExecutor(true).getExternalExecutorExists(),
                "python not found in environment ${PATH}"
        );
    }

    @Test(groups = {"Python"}, dependsOnMethods = "testPythonExists")
    public void testExecuteAsCommand() {
        final PythonScriptExecutor pythonExecutor = new PythonScriptExecutor(true);
        final boolean ret = pythonExecutor.executeCommand(HELLO_WORLD_SCRIPT, null, null);

        Assert.assertTrue(ret, "Python exec failed");
    }

    @Test(groups = {"Python"}, dependsOnMethods = "testPythonExists")
    public void testExecuteAsModule() {
        // use builtin module "random"
        final PythonScriptExecutor executor = new PythonScriptExecutor(true);
        final boolean ret = executor.executeModule("random", null, null);

        Assert.assertTrue(ret, "Python exec failed");
    }

    @Test(groups = {"Python"}, dependsOnMethods = "testPythonExists")
    public void testExecuteAsScript() {
        final File scriptFile = writeTemporaryScriptFile(HELLO_WORLD_SCRIPT, PythonScriptExecutor.pyExtension);

        final PythonScriptExecutor executor = new PythonScriptExecutor(true);
        final boolean ret = executor.executeScript(scriptFile.getAbsolutePath(), null, null);

        Assert.assertTrue(ret, "Python exec failed");
    }

    @Test(groups = {"Python"}, dependsOnMethods = "testPythonExists")
    public void testExecuteAsRawArgs() {
        final PythonScriptExecutor pythonExecutor = new PythonScriptExecutor(true);
        final boolean ret = pythonExecutor.executeArgs(new ArrayList<String>(Arrays.asList("-c", HELLO_WORLD_SCRIPT)));

        Assert.assertTrue(ret, "Python exec failed");
    }

    @Test(groups = {"Python"}, dependsOnMethods = "testPythonExists")
    public void testExecuteAsRawArgsSerial() {
        final PythonScriptExecutor pythonExecutor = new PythonScriptExecutor(true);
        final List<String> args = new ArrayList<>(Arrays.asList("-c", HELLO_WORLD_SCRIPT));

        Assert.assertTrue(pythonExecutor.executeArgs(args), "First Python exec failed");
        Assert.assertTrue(pythonExecutor.executeArgs(args), "Second Python exec failed");
    }

    @Test(groups = {"Python"}, dependsOnMethods = "testPythonExists")
    public void testExecuteAsRawArgsWithPythonArguments() {
        final PythonScriptExecutor pythonExecutor = new PythonScriptExecutor(true);
        final boolean ret = pythonExecutor.executeArgs(new ArrayList<>(Arrays.asList("-v", "-c", HELLO_WORLD_SCRIPT )));

        Assert.assertTrue(ret, "Python exec failed");
    }

    @Test(groups = {"Python"}, dependsOnMethods = "testPythonExists")
    public void testExecuteAsScriptWithScriptArguments() {
        final String SCRIPT_WITH_ARGUMENTS =
                "import fileinput\n" +
                "for line in fileinput.input():\n" +
                "    print(line)\n";
        final File scriptFile = writeTemporaryScriptFile(SCRIPT_WITH_ARGUMENTS, PythonScriptExecutor.pyExtension);
        final PythonScriptExecutor executor = new PythonScriptExecutor(true);
        final boolean ret = executor.executeScript(
                scriptFile.getAbsolutePath(),
                null,
                Collections.singletonList(scriptFile.getAbsolutePath()));

       Assert.assertTrue(ret, "Python exec failed");
    }

    @Test(groups = {"Python"}, dependsOnMethods = "testPythonExists", expectedExceptions = PythonScriptExecutorException.class)
    public void testNonExistentScriptException() throws IOException {
        final PythonScriptExecutor executor = new PythonScriptExecutor(true);
        executor.executeScript(
                BaseTest.getSafeNonExistentFile("Nonexistent" + PythonScriptExecutor.pyExtension).getAbsolutePath(),
                null,
                null);
    }

    @Test(groups = {"Python"}, dependsOnMethods = "testPythonExists")
    public void testNonExistentScriptNoException() {
        final PythonScriptExecutor executor = new PythonScriptExecutor(true);
        executor.setIgnoreExceptions(true);
        Assert.assertFalse(
                executor.executeScript(
                        BaseTest.getSafeNonExistentFile("Nonexistent" + PythonScriptExecutor.pyExtension).getAbsolutePath(),
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
