package org.broadinstitute.hellbender.utils.runtime;

import htsjdk.samtools.util.BufferedLineReader;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.concurrent.TimeoutException;

// NOTE: some Python implementations have bugs where sometimes the prompt is printed to stdout,
// and sometimes to stderr:
//
//  https://bugs.python.org/issue17620
//  https://bugs.python.org/issue1927
//
// NOTE: TestNG has a bug where it throws ArrayIndexOutOfBoundsException instead of TimeoutException
// exception when the test time exceeds the timeOut threshold. This is fixed but not yet released:
//
//  https://github.com/cbeust/testng/issues/1493
//
public class StreamingProcessControllerUnitTest extends BaseTest {
    private final static String NL = System.lineSeparator();
    private final static String PYTHON_PROMPT = ">>> ";

    @Test(groups = "python", timeOut = 10000)
    public void testSerialCommands() throws TimeoutException, IOException {
        // start an interactive Python session with unbuffered IO
        final ProcessSettings processSettings = new ProcessSettings(new String[] {"python", "-i", "-u"});

        final StreamingProcessController controller = new StreamingProcessController(processSettings, PYTHON_PROMPT);
        Assert.assertTrue(controller.start());

        // consume the Python startup banner debris, but don't try validate it
        controller.getProcessOutputByPrompt();

        controller.writeProcessInput("x = 37" + NL);
        ProcessOutput response = controller.getProcessOutputByPrompt();
        StreamingPythonTestUtils.assertPythonPrompt(response, PYTHON_PROMPT);

        controller.writeProcessInput("x = x + 2" + NL);
        response = controller.getProcessOutputByPrompt();
        StreamingPythonTestUtils.assertPythonPrompt(response, PYTHON_PROMPT);

        final File tempFile = createTempFile("testPythonStdoutSerial", "txt");
        final String writeOutput = String.format("fd=open('%s', 'w')\nfd.write(str(x) + '\\n')\nfd.close()", tempFile.getAbsolutePath());

        controller.writeProcessInput(writeOutput + NL);
        response = controller.getProcessOutputByPrompt();
        StreamingPythonTestUtils.assertPythonPrompt(response, PYTHON_PROMPT);

        try (final FileInputStream fis = new FileInputStream(tempFile);
             final BufferedLineReader br = new BufferedLineReader(fis)) {
            Assert.assertEquals(br.readLine(), "39");
        }

        controller.terminate();
        Assert.assertFalse(controller.getProcess().isAlive());
    }

    @Test(timeOut = 10000)
    public void testMultipleParallelStreamingControllers() throws TimeoutException {
        final ProcessSettings catProcessSettings = new ProcessSettings(new String[] {"cat"});
        final StreamingProcessController catController = new StreamingProcessController(catProcessSettings);
        Assert.assertTrue(catController.start());

        final ProcessSettings teeProcessSettings = new ProcessSettings(new String[] {"tee"});
        final StreamingProcessController teeController = new StreamingProcessController(teeProcessSettings);
        Assert.assertTrue(teeController.start());

        final String catString = "send to cat\n";
        final String teeString = "send to tee\n";

        catController.writeProcessInput(catString);
        teeController.writeProcessInput(teeString);

        final ProcessOutput catResponseLine = catController.getProcessOutputByLine();
        final ProcessOutput teeResponseLine = teeController.getProcessOutputByLine();

        StreamingPythonTestUtils.assertResponseOutput(catResponseLine, catString, true);
        StreamingPythonTestUtils.assertResponseOutput(teeResponseLine, teeString, true);

        catController.terminate();
        teeController.terminate();

        Assert.assertFalse(catController.getProcess().isAlive());
        Assert.assertFalse(teeController.getProcess().isAlive());
    }

    @Test(groups = "python", timeOut = 10000)
    public void testStartupCommandExecution() throws TimeoutException, IOException {
        final String writeOutTemplate = "fd=open('%s', 'w')\nfd.write('some output\\n')\nfd.close()\n";
        final File tempFile = createTempFile("streamingControllerStartupCommand", ".txt");
        final String writeOutScript = String.format(writeOutTemplate, tempFile.getAbsolutePath());

        final ProcessSettings processSettings = new ProcessSettings(new String[] {"python", "-i", "-u", "-c", writeOutScript});

        final StreamingProcessController controller = new StreamingProcessController(processSettings, PYTHON_PROMPT);

        Assert.assertTrue(controller.start());
        controller.terminate();
        Assert.assertFalse(controller.getProcess().isAlive());

        try (final FileInputStream fis = new FileInputStream(tempFile);
             final BufferedLineReader br = new BufferedLineReader(fis)) {
            Assert.assertEquals(br.readLine(), "some output");
        }
    }

    // using invocationTimeOut seems to be buggy and causes this to stop after one invocation
    //@Test(invocationCount = 10, invocationTimeOut = 10000)
    @Test(groups = "python", timeOut = 10000) // make sure the test timeout exceeds the controller timeout, since we want to trigger that
    public void testIsOutputAvailable() throws TimeoutException, InterruptedException {
        // start an interactive Python session with unbuffered IO
        final ProcessSettings processSettings = new ProcessSettings(new String[] {"python", "-i", "-u"});

        final StreamingProcessController controller = new StreamingProcessController(processSettings, PYTHON_PROMPT);
        Assert.assertTrue(controller.start());
        // consume the Python startup banner debris, but don't try validate it
        controller.getProcessOutputByPrompt();

        Assert.assertFalse(controller.isOutputAvailable());
        controller.writeProcessInput("a = 3" + NL);
        while (!controller.isOutputAvailable()) {
            // TestNG will throw if we exceed the timeout
        };
        Assert.assertTrue(controller.isOutputAvailable());

        controller.terminate();
        Assert.assertFalse(controller.getProcess().isAlive());
    }

    @Test(groups = "python", timeOut = 10000)
    public void testStderrOutput() throws TimeoutException {
        // test write to stderr from python
        // start an interactive Python session with unbuffered IO
        final ProcessSettings processSettings = new ProcessSettings(new String[] {"python", "-i", "-u"});

        final StreamingProcessController controller = new StreamingProcessController(processSettings, PYTHON_PROMPT);
        Assert.assertTrue(controller.start());
        // consume the Python startup banner debris, but don't try validate it
        controller.getProcessOutputByPrompt();

        controller.writeProcessInput("import sys" + NL + "sys.stderr.write('error output to stderr\\n')" + NL);
        ProcessOutput po = controller.getProcessOutputByPrompt();

        Assert.assertNotNull(po.getStderr());
        Assert.assertNotNull(po.getStderr().getBufferString());
        Assert.assertNotNull(po.getStderr().getBufferString().contains("error output to stderr"));

        controller.terminate();
        Assert.assertFalse(controller.getProcess().isAlive());
    }

    @Test(groups = "python", timeOut = 10000)
    public void testStderrRedirect() throws TimeoutException {
        // test write to stderr from python
        // start an interactive Python session with unbuffered IO
        final ProcessSettings processSettings = new ProcessSettings(new String[] {"python", "-i", "-u"});

        // redirect the process' stderr to stdout
        processSettings.setRedirectErrorStream(true);

        final StreamingProcessController controller = new StreamingProcessController(processSettings, PYTHON_PROMPT);
        Assert.assertNotNull(controller.start());
        // consume the Python startup banner debris, but don't try validate it
        controller.getProcessOutputByPrompt();

        // write to stderr, but we expect to get it from stdout due to redirection
        controller.writeProcessInput("import sys" + NL + "sys.stderr.write('error output to stderr\\n')" + NL);
        ProcessOutput po = controller.getProcessOutputByPrompt();

        Assert.assertNotNull(po.getStdout());
        Assert.assertNotNull(po.getStdout().getBufferString());
        Assert.assertNotNull(po.getStdout().getBufferString().contains("error output to stderr"));

        controller.terminate();
        Assert.assertFalse(controller.getProcess().isAlive());
    }

    @Test(timeOut = 10000)
    public void testFIFOLifetime() throws TimeoutException {
        // cat is a red herring here; we're just testing that a FIFO is created, and then deleted after termination
        final ProcessSettings catProcessSettings = new ProcessSettings(new String[] {"cat"});
        final StreamingProcessController catController = new StreamingProcessController(catProcessSettings);
        catController.start();

        final File fifo = catController.createFIFO();
        Assert.assertTrue(fifo.exists());

        catController.terminate();

        final File fifoParent = fifo.getParentFile();
        Assert.assertFalse(fifo.exists());
        Assert.assertFalse(fifoParent.exists());
        Assert.assertFalse(catController.getProcess().isAlive());
    }

    // this needs to have a testNG timeOut that is longer than the timeout built in to the StreamingProcessController,
    // since we want to trigger the latter first, to ensure that we hit the builtin one first
    @Test(groups = "python", timeOut = 10000, expectedExceptions=TimeoutException.class)
    public void testPromptTimeout() throws TimeoutException {
        // start an interactive Python session with unbuffered IO
        final ProcessSettings processSettings = new ProcessSettings(new String[] {"python", "-i", "-u"});

        final StreamingProcessController controller = new StreamingProcessController(processSettings, "bogus");
        Assert.assertTrue(controller.start());

        try {
            // this will hang waiting for a prompt to appear, until timeout, since the prompt is bogus
            controller.getProcessOutputByPrompt();
        }
        finally {
            controller.terminate();
            Assert.assertFalse(controller.getProcess().isAlive());
        }
    }

    @Test(groups = "python", expectedExceptions=IllegalStateException.class)
    public void testInvalidPromptSynchronization() throws TimeoutException {
        final ProcessSettings pythonProcessSettings = new ProcessSettings(new String[] {"python", "-i", "-u"});

        final StreamingProcessController pythonController = new StreamingProcessController(pythonProcessSettings);
        pythonController.start();

        try {
            pythonController.getProcessOutputByPrompt();
        }
        finally {
            pythonController.terminate();
            Assert.assertFalse(pythonController.getProcess().isAlive());
        }
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testRedundantStart() throws TimeoutException {
        final ProcessSettings catProcessSettings = new ProcessSettings(new String[] {"cat"});
        final StreamingProcessController catController = new StreamingProcessController(catProcessSettings);
        Assert.assertTrue(catController.start());

        try {
            catController.start();
        } finally {
            catController.terminate();
            Assert.assertFalse(catController.getProcess().isAlive());
        }
    }

}
