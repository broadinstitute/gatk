package org.broadinstitute.hellbender.utils.python;

import htsjdk.samtools.util.BufferedLineReader;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.runtime.AsynchronousStreamWriterService;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;
import org.broadinstitute.hellbender.utils.runtime.StreamingPythonTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

// NOTE: Python has bugs where it sometimes prints the prompt to stdout, and sometimes to stderr:
//
// https://bugs.python.org/issue17620
// https://bugs.python.org/issue1927
//
// Beware TestNG has a bug where it throws ArrayIndexOutOfBoundsException instead of TimeoutException
// exception when the test time exceeds the timeOut threshold. This is fixed but not yet released:
//
//  https://github.com/cbeust/testng/issues/1493
//
public class StreamingPythonScriptExecutorUnitTest extends GATKBaseTest {
    private static final String NL = System.lineSeparator();

    @DataProvider(name="supportedPythonVersions")
    public Object[][] getSupprtedPythonExecutableNames() {
        return new Object[][] {
                { PythonScriptExecutor.PythonExecutableName.PYTHON },
                { PythonScriptExecutor.PythonExecutableName.PYTHON3 },
        };
    }

    @Test(groups = "python", dataProvider="supportedPythonVersions")
    public void testPythonExists(final PythonScriptExecutor.PythonExecutableName executableName) {
        Assert.assertTrue(
                new StreamingPythonScriptExecutor(executableName,true).externalExecutableExists(),
                "Python not found in environment ${PATH}"
        );
    }

    @Test(groups = "python", dataProvider="supportedPythonVersions", dependsOnMethods = "testPythonExists", timeOut=10000)
    public void testExecuteCommand(final PythonScriptExecutor.PythonExecutableName executableName) throws IOException {
        // create a temporary output file
        final File tempFile = createTempFile("pythonExecuteCommandTest", "txt");
        final String WRITE_FILE_SCRIPT =
                String.format(
                        "with open('%s', 'w') as tempFile:\n    tempFile.write('Hello world')" + NL + NL,
                        tempFile.getAbsolutePath());
        final String CLOSE_FILE_SCRIPT = "tempFile.close()" + NL;

        final StreamingPythonScriptExecutor streamingPythonExecutor =
                new StreamingPythonScriptExecutor(executableName,true);
        Assert.assertNotNull(streamingPythonExecutor);

        Assert.assertTrue(streamingPythonExecutor.start(Collections.emptyList()));

        try {
            streamingPythonExecutor.sendSynchronousCommand(WRITE_FILE_SCRIPT);
            streamingPythonExecutor.sendSynchronousCommand(CLOSE_FILE_SCRIPT);
        }
        finally {
            streamingPythonExecutor.terminate();
            Assert.assertFalse(streamingPythonExecutor.getProcess().isAlive());
        }

        // read the temp file in and validate
        try (final FileInputStream fis= new FileInputStream(tempFile);
             final BufferedLineReader br = new BufferedLineReader(fis)) {
            Assert.assertEquals(br.readLine(), "Hello world");
        }

    }

    @Test(groups = "python", dataProvider="supportedPythonVersions", dependsOnMethods = "testPythonExists")
    public void testStderrOutput(final PythonScriptExecutor.PythonExecutableName executableName) {
        // test write to stderr from python
        final StreamingPythonScriptExecutor streamingPythonExecutor =
                new StreamingPythonScriptExecutor(executableName,true);
        Assert.assertNotNull(streamingPythonExecutor);
        Assert.assertTrue(streamingPythonExecutor.start(Collections.emptyList()));

        try {
            streamingPythonExecutor.sendSynchronousCommand("import sys" + NL);
            final ProcessOutput po = streamingPythonExecutor.sendSynchronousCommand(
                    "sys.stderr.write('error output to stderr\\n')" + NL);

            Assert.assertNotNull(po.getStderr());
            Assert.assertNotNull(po.getStderr().getBufferString());
            Assert.assertNotNull(po.getStderr().getBufferString().contains("error output to stderr"));
        }
        finally {
            streamingPythonExecutor.terminate();
            Assert.assertFalse(streamingPythonExecutor.getProcess().isAlive());
        }
    }

    @Test(groups = "python", dataProvider="supportedPythonVersions", dependsOnMethods = "testPythonExists")
    public void testTerminateWhilePythonBlocked(final PythonScriptExecutor.PythonExecutableName executableName) {
        // Test termination on a Python process that is blocked on I/O to ensure that we don't leave
        // zombie Python processes around. This unfinished code block will cause Python to hang with
        // a continuation prompt writing for more input, and must be sent to Python asynchronously.
        final String UNTERMINATED_SCRIPT = "with open('nonexistent.txt', 'w') as foo:" + NL;
        final StreamingPythonScriptExecutor streamingPythonExecutor =
                new StreamingPythonScriptExecutor(executableName, true);
        Assert.assertTrue(streamingPythonExecutor.start(Collections.emptyList()));

        try {
            // send the output asynchronously so *we* don't block
            streamingPythonExecutor.sendAsynchronousCommand(UNTERMINATED_SCRIPT);
        }
        finally {
            streamingPythonExecutor.terminate();
            Assert.assertFalse(streamingPythonExecutor.getProcess().isAlive());
        }
    }

    @Test(groups = "python", dataProvider="supportedPythonVersions", dependsOnMethods = "testPythonExists", timeOut = 10000)
    public void testIsOutputAvailable(final PythonScriptExecutor.PythonExecutableName executableName) {

        final String UNTERMINATED_SCRIPT = "with open('nonexistent.txt', 'w') as foo:" + NL;
        final StreamingPythonScriptExecutor streamingPythonExecutor =
                new StreamingPythonScriptExecutor(executableName, true);
        Assert.assertTrue(streamingPythonExecutor.start(Collections.emptyList()));

        try {
            // send the output asynchronously so *we* don't block
            streamingPythonExecutor.sendAsynchronousCommand(UNTERMINATED_SCRIPT);
            Assert.assertFalse(streamingPythonExecutor.isOutputAvailable());

            // now terminate the line; we should get a prompt
            streamingPythonExecutor.sendAsynchronousCommand("\n\n");

            // loop waiting for output; TestNG will timeout if this spins too long
            while (!streamingPythonExecutor.isOutputAvailable()) { }
            Assert.assertTrue(streamingPythonExecutor.isOutputAvailable());
        }
        finally {
            streamingPythonExecutor.terminate();
            Assert.assertFalse(streamingPythonExecutor.getProcess().isAlive());
        }
    }

    @Test(groups = "python", dataProvider="supportedPythonVersions", dependsOnMethods = "testPythonExists", expectedExceptions = RuntimeException.class)
    public void testTryTerminatedController(final PythonScriptExecutor.PythonExecutableName executableName) {
        final StreamingPythonScriptExecutor streamingPythonExecutor =
                new StreamingPythonScriptExecutor(executableName, true);
        Assert.assertTrue(streamingPythonExecutor.start(Collections.emptyList()));

        streamingPythonExecutor.terminate();
        Assert.assertFalse(streamingPythonExecutor.getProcess().isAlive());

        streamingPythonExecutor.sendAsynchronousCommand("bogus command = controller is already terminated\n");
    }

    @Test(groups = "python", dataProvider="supportedPythonVersions", dependsOnMethods = "testPythonExists", expectedExceptions = IllegalArgumentException.class)
    public void testRequireNewlineTerminatedCommand(final PythonScriptExecutor.PythonExecutableName executableName) {
        final String NO_NEWLINE_SCRIPT = "print 'hello, world'";
        final StreamingPythonScriptExecutor streamingPythonExecutor =
                new StreamingPythonScriptExecutor(executableName, true);
        Assert.assertTrue(streamingPythonExecutor.start(Collections.emptyList()));
        try {
            streamingPythonExecutor.sendAsynchronousCommand(NO_NEWLINE_SCRIPT);
        }
        finally {
            streamingPythonExecutor.terminate();
            Assert.assertFalse(streamingPythonExecutor.getProcess().isAlive());
        }
    }

    @Test(groups = "python", dataProvider="supportedPythonVersions", dependsOnMethods = "testPythonExists", timeOut=10000, invocationCount = 20)
    public void testMultipleFIFORoundTrips(final PythonScriptExecutor.PythonExecutableName executableName) throws IOException {
        // Python script statements to open, read, and write to and from the FIFO and temporary files
        final String PYTHON_INITIALIZE_COUNT    = "count = 0" + NL;
        final String PYTHON_OPEN_FIFO           = "fifoFile = open('%s', 'r')" + NL;
        final String PYTHON_OPEN_TEMP_FILE      = "tempFile = open('%s', 'w')" + NL;
        final String PYTHON_TRANSFER_FIFO_TO_TEMP_FILE = "tempFile.write(fifoFile.readline())" + NL;
        final String PYTHON_UPDATE_COUNT        = "count = count + 1" + NL;
        final String PYTHON_WRITE_COUNT_TO_TEMP_FILE = "tempFile.write(str(count))" + NL;
        final String PYTHON_CLOSE_TEMP_FILE     = "tempFile.close()" + NL;
        final String PYTHON_CLOSE_FIFO          = "fifoFile.close()" + NL;

        final String ROUND_TRIP_DATA = "round trip data (%s)" + NL;

        final List<String> linesWrittenToFIFO = new ArrayList<>();

        final StreamingPythonScriptExecutor streamingPythonExecutor =
                new StreamingPythonScriptExecutor(executableName, true);
        Assert.assertNotNull(streamingPythonExecutor);
        // ignore/drop the startup banner
        Assert.assertTrue(streamingPythonExecutor.start(Collections.emptyList()));

        // create a temporary file
        final File tempFile = createTempFile("pythonRoundTripTest", "txt");

        try {
            // initialize a counter in Python to keep track of # of round trips
            ProcessOutput response = streamingPythonExecutor.sendSynchronousCommand(PYTHON_INITIALIZE_COUNT);
            StreamingPythonTestUtils.assertPythonPrompt(response, StreamingPythonScriptExecutor.PYTHON_PROMPT);

            // ask python to open the temp file
            response = streamingPythonExecutor.sendSynchronousCommand(String.format(PYTHON_OPEN_TEMP_FILE, tempFile.getAbsolutePath()));
            StreamingPythonTestUtils.assertPythonPrompt(response, StreamingPythonScriptExecutor.PYTHON_PROMPT);

            // obtain a FIFO and ask Python to open it; the Python process will block waiting for a writer
            // to open the FIFO, so send the command asynchronously, and don't request output until after we
            // open FIFO for writing
            final File fifoFile = streamingPythonExecutor.getFIFOForWrite();
            streamingPythonExecutor.sendAsynchronousCommand(String.format(PYTHON_OPEN_FIFO, fifoFile.getAbsolutePath()));

            // now open the FIFO for write; this will also block until Python executes the line above and opens the FIFO
            try (final FileOutputStream fileWriter = new FileOutputStream(fifoFile)) {

                // synchronize on the prompt resulting from the python fifo open statement
                response = streamingPythonExecutor.getAccumulatedOutput();
                StreamingPythonTestUtils.assertPythonPrompt(response, StreamingPythonScriptExecutor.PYTHON_PROMPT);

                // now write data to the fifo and round trip to python ROUND_TRIP_COUNT times
                final int ROUND_TRIP_COUNT = 100;
                for (int i = 0; i < ROUND_TRIP_COUNT; i++) {

                    // write some data to the fifo
                    final String roundTripLine = String.format(ROUND_TRIP_DATA, i);
                    fileWriter.write(roundTripLine.getBytes());
                    fileWriter.flush();

                    // tell Python to read the line we just wrote to the FIFO, and write it out to the temp file;
                    // synchronize on a prompt only since python 3.x writes the # of characters written to the file
                    // to stdout, but Python 2.7 does not
                    streamingPythonExecutor.sendAsynchronousCommand(PYTHON_TRANSFER_FIFO_TO_TEMP_FILE);
                    response = streamingPythonExecutor.getAccumulatedOutput();
                    StreamingPythonTestUtils.assertPythonPrompt(response, StreamingPythonScriptExecutor.PYTHON_PROMPT);

                    // tell Python to update it's count
                    streamingPythonExecutor.sendAsynchronousCommand(PYTHON_UPDATE_COUNT);
                    response = streamingPythonExecutor.getAccumulatedOutput();
                    StreamingPythonTestUtils.assertPythonPrompt(response, StreamingPythonScriptExecutor.PYTHON_PROMPT);

                    // keep track of everything we write to the FIFO so we can validate at the end of the test
                    linesWrittenToFIFO.add(roundTripLine);
                }
            }

            // tell Python to write the final count to the temp file
            streamingPythonExecutor.sendAsynchronousCommand(PYTHON_WRITE_COUNT_TO_TEMP_FILE);
            response = streamingPythonExecutor.getAccumulatedOutput();
            StreamingPythonTestUtils.assertPythonPrompt(response, StreamingPythonScriptExecutor.PYTHON_PROMPT);

            // tell Python to close the temp file
            streamingPythonExecutor.sendAsynchronousCommand(PYTHON_CLOSE_TEMP_FILE);
            response = streamingPythonExecutor.getAccumulatedOutput();
            StreamingPythonTestUtils.assertPythonPrompt(response, StreamingPythonScriptExecutor.PYTHON_PROMPT);

            // tell Python to close the FIFO file
            streamingPythonExecutor.sendAsynchronousCommand(PYTHON_CLOSE_FIFO);
            response = streamingPythonExecutor.getAccumulatedOutput();
            StreamingPythonTestUtils.assertPythonPrompt(response, StreamingPythonScriptExecutor.PYTHON_PROMPT);
        }
        finally {
            streamingPythonExecutor.terminate();
            Assert.assertFalse(streamingPythonExecutor.getProcess().isAlive());
        }

        // read the temp file in and make sure everything written to the FIFO was round-tripped by Python,
        // including the round-trip count
        try (final FileInputStream fis= new FileInputStream(tempFile);
             final BufferedLineReader br = new BufferedLineReader(fis)) {
            linesWrittenToFIFO.forEach(expectedLine -> Assert.assertEquals(br.readLine() + NL, expectedLine));
        }
    }

    @Test(groups = "python", dataProvider="supportedPythonVersions", dependsOnMethods = "testPythonExists", timeOut=10000, invocationCount = 10)
    public void testAsyncWriteService(final PythonScriptExecutor.PythonExecutableName executableName) throws IOException, InterruptedException {
        // Python script statements to open, read, and write to and from the FIFO and temporary files
        final String PYTHON_OPEN_FIFO           = "fifoFile = open('%s', 'r')" + NL;
        final String PYTHON_OPEN_TEMP_FILE      = "tempFile = open('%s', 'w')" + NL;
        final String PYTHON_TRANSFER_FIFO_TO_TEMP_FILE = "for i in range(%s):\n    tempFile.write(fifoFile.readline())" + NL + NL;
        final String PYTHON_CLOSE_TEMP_FILE     = "tempFile.close()" + NL;
        final String PYTHON_CLOSE_FIFO          = "fifoFile.close()" + NL;

        final String ROUND_TRIP_DATA = "round trip data (%s)" + NL;

        final List<String> linesWrittenToFIFO = new ArrayList<>();

        final StreamingPythonScriptExecutor streamingPythonExecutor =
                new StreamingPythonScriptExecutor(executableName, false);
        Assert.assertNotNull(streamingPythonExecutor);
        // ignore/drop the startup banner
        Assert.assertTrue(streamingPythonExecutor.start(Collections.emptyList()));

        // create a temporary file
        final File tempFile = createTempFile("pythonRoundTripTest", "txt");

        AsynchronousStreamWriterService<String> FIFOWriteService = null;
        try {
            // ask python to open the temp file
            ProcessOutput response = streamingPythonExecutor.sendSynchronousCommand(String.format(PYTHON_OPEN_TEMP_FILE, tempFile.getAbsolutePath()));
            StreamingPythonTestUtils.assertPythonPrompt(response, StreamingPythonScriptExecutor.PYTHON_PROMPT);

            // obtain a FIFO and ask Python to open it; the Python process will block waiting for a writer
            // to open the FIFO, so send the command asynchronously, and don't request output until after we
            // open FIFO for writing
            final File fifoFile = streamingPythonExecutor.getFIFOForWrite();
            streamingPythonExecutor.sendAsynchronousCommand(String.format(PYTHON_OPEN_FIFO, fifoFile.getAbsolutePath()));

            final int ROUND_TRIP_COUNT = 100*1024;
            final int syncFrequency = 1024;
            // now open the FIFO for write; this will also block until Python executes the line above and opens the FIFO
            try (final FileOutputStream fileWriter = new FileOutputStream(fifoFile)) {

                // synchronize on the prompt resulting from the python fifo open statement
                response = streamingPythonExecutor.getAccumulatedOutput();
                StreamingPythonTestUtils.assertPythonPrompt(response, StreamingPythonScriptExecutor.PYTHON_PROMPT);

                // get an async writer; write data to the fifo and round trip to python ROUND_TRIP_COUNT times
                FIFOWriteService = streamingPythonExecutor.getAsynchronousStreamWriterService(fileWriter, AsynchronousStreamWriterService.stringSerializer);
                List<String> fifoData = new ArrayList<>(syncFrequency);
                int count = 0;

                for (int i = 0; i < ROUND_TRIP_COUNT; i++) {
                    if (i != 0 && (i % syncFrequency) == 0) {
                        // wait for the last batch to complete before we start a new one
                        Future<Integer> batchResult = FIFOWriteService.waitForPreviousBatchCompletion(1, TimeUnit.SECONDS);
                        streamingPythonExecutor.sendAsynchronousCommand(String.format(PYTHON_TRANSFER_FIFO_TO_TEMP_FILE, count));
                        FIFOWriteService.startAsynchronousBatchWrite(fifoData);
                        count = 0;
                        fifoData = new ArrayList<>(syncFrequency);
                    }

                    // save the line to be written to the fifo
                    final String roundTripLine = String.format(ROUND_TRIP_DATA, i);
                    fifoData.add(roundTripLine);
                    // keep track of everything we write to the FIFO so we can validate at the end of the test
                    linesWrittenToFIFO.add(roundTripLine);
                    count++;
                }

                // write the final batch
                FIFOWriteService.waitForPreviousBatchCompletion(1, TimeUnit.SECONDS);
                if (fifoData.size() != 0) {
                    streamingPythonExecutor.sendAsynchronousCommand(String.format(PYTHON_TRANSFER_FIFO_TO_TEMP_FILE, count));
                    FIFOWriteService.startAsynchronousBatchWrite(fifoData);
                    FIFOWriteService.waitForPreviousBatchCompletion(1, TimeUnit.SECONDS);
                }
                Assert.assertTrue(FIFOWriteService.terminate());
            }

            Assert.assertEquals(linesWrittenToFIFO.size(), ROUND_TRIP_COUNT);

            // tell Python to close the temp file
            streamingPythonExecutor.sendAsynchronousCommand(PYTHON_CLOSE_TEMP_FILE);
            response = streamingPythonExecutor.getAccumulatedOutput();
            StreamingPythonTestUtils.assertPythonPrompt(response, StreamingPythonScriptExecutor.PYTHON_PROMPT);

            // tell Python to close the FIFO file
            streamingPythonExecutor.sendAsynchronousCommand(PYTHON_CLOSE_FIFO);
            response = streamingPythonExecutor.getAccumulatedOutput();
            StreamingPythonTestUtils.assertPythonPrompt(response, StreamingPythonScriptExecutor.PYTHON_PROMPT);
        }
        finally {
            if (FIFOWriteService != null) {
                Assert.assertTrue(FIFOWriteService.terminate());
            }
            streamingPythonExecutor.terminate();
            Assert.assertFalse(streamingPythonExecutor.getProcess().isAlive());
        }

        // read the temp file in and make sure everything written to the FIFO was round-tripped by Python,
        // including the round-trip count
        try (final FileInputStream fis= new FileInputStream(tempFile);
             final BufferedLineReader br = new BufferedLineReader(fis)) {
            linesWrittenToFIFO.forEach(expectedLine -> Assert.assertEquals(br.readLine() + NL, expectedLine));
        }
    }

    @Test(groups = "python", dataProvider="supportedPythonVersions", dependsOnMethods = "testPythonExists", expectedExceptions = PythonScriptExecutorException.class)
    public void testPythonAssert(final PythonScriptExecutor.PythonExecutableName executableName) {
        final StreamingPythonScriptExecutor streamingPythonExecutor =
                new StreamingPythonScriptExecutor(executableName,true);
        Assert.assertNotNull(streamingPythonExecutor);
        Assert.assertTrue(streamingPythonExecutor.start(Collections.emptyList()));
        streamingPythonExecutor.getAccumulatedOutput();

        try {
            ProcessOutput po = streamingPythonExecutor.sendSynchronousCommand("assert 0" + NL);

            Assert.assertNotNull(po.getStderr());
            Assert.assertNotNull(po.getStderr().getBufferString());
        }
        finally {
            streamingPythonExecutor.terminate();
            Assert.assertFalse(streamingPythonExecutor.getProcess().isAlive());
        }
    }

}
