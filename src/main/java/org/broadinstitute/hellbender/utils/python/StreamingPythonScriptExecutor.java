package org.broadinstitute.hellbender.utils.python;

import com.google.common.annotations.VisibleForTesting;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.runtime.*;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Future;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Python executor used to interact with a cooperative, keep-alive Python process. The executor issues commands
 * to call Python functions in the {@code tool} module in {@code gatktool} Python package. These include functions
 * for managing an acknowledgement FIFO that is used to signal completion of Python commands, and a data FIFO that
 * can be used to stream data to Python.
 *
 *  - construct the executor
 *  - start the remote process ({@link #start}.
 *  - optionally call {@code #getStreamWriter} to initialize and create a data transfer fifo.
 *  - send one or more synchronous or asynchronous commands to be executed in Python
 *  - optionally send data one or more times of type {@ocde T} through the async writer
 *  - execute python code to close the data fifo
 *  - terminate the executor {@link #terminate}
 *
 * Guidelines for writing GATK tools that use Python interactively:
 *
 *   - Program correctness should not rely on consumption of anything written by Python to stdout/stderr. All
 *     data should be transferred through the stream writer or a file.
 *   - Python code should write errors to stderr.
 *   - Prefer single line commands that run a script, vs. multi-line Python code embedded in Java
 *   - Terminate commands with a newline.
 *   - Try not to be chatty (maximize use of the fifo buffer by writing to it in batches before reading from Python)
 *
 * @param <T> type of data that will be streamed to the Python process
 */
public class StreamingPythonScriptExecutor<T> extends PythonExecutorBase {
    private static final Logger logger = LogManager.getLogger(StreamingPythonScriptExecutor.class);
    private static final String NL = System.lineSeparator();

    private final List<String> curatedCommandLineArgs = new ArrayList<>();

    private StreamingProcessController spController;
    private ProcessSettings processSettings;

    private File dataTransferFIFOFile;
    private FileOutputStream dataTransferFIFOWriter;
    private AsynchronousStreamWriter<T> asyncWriter;

    private File profileResults;

    // Python commands that are executed in the companion python process. The functions called
    // here live in the {@code tool} module in {@code gatktool} Python package.
    private final static String PYTHON_IMPORT_GATK = "from gatktool import tool" + NL;
    private final static String PYTHON_INITIALIZE_GATK = "tool.initializeGATK('%s')" + NL;
    private final static String PYTHON_START_PROFILING = "tool.startProfiling()" + NL;
    private final static String PYTHON_TERMINATE_GATK = "tool.terminateGATK()" + NL;
    private final static String PYTHON_INITIALIZE_DATA_FIFO = "tool.initializeDataFIFO('%s')" + NL;
    private final static String PYTHON_CLOSE_DATA_FIFO = "tool.closeDataFIFO()" + NL;
    private final static String PYTHON_SEND_ACK_REQUEST = "tool.sendAck()" + NL;
    private final static String PYTHON_END_PROFILING = "tool.endProfiling('%s')" + NL;

    // keep track of when an ack request has been made and reject attempts to send another ack
    // request until the previous one has been handled
    private boolean isAckRequestOutstanding = false;

    /**
     * The start method must be called to actually start the remote executable.
     *
     * @param ensureExecutableExists throw if the python executable cannot be located
     */
    public StreamingPythonScriptExecutor(boolean ensureExecutableExists) {
        this(PythonExecutableName.PYTHON, ensureExecutableExists);
    }

    /**
     * The start method must be called to actually start the remote executable.
     *
     * @param pythonExecutableName name of the python executable to start
     * @param ensureExecutableExists throw if the python executable cannot be found
     */
    public StreamingPythonScriptExecutor(final PythonExecutableName pythonExecutableName, final boolean ensureExecutableExists) {
        super(pythonExecutableName, ensureExecutableExists);
    }

    /**
     * Start the Python process.
     *
     * @param pythonProcessArgs args to be passed to the python process
     * @return true if the process is successfully started
     */
    public boolean start(final List<String> pythonProcessArgs) {
        return start(pythonProcessArgs, false, null);
    }

    /**
     * Start the Python process.
     *
     * @param pythonProcessArgs args to be passed to the python process
     * @param enableJournaling true to enable Journaling, which records all interprocess IO to a file. This is
     *                         expensive and should only be used for debugging purposes.
     * @return true if the process is successfully started
     */
    public boolean start(final List<String> pythonProcessArgs, final boolean enableJournaling, final File profileResults) {
        this.profileResults = profileResults;
        final List<String> args = new ArrayList<>();
        args.add(externalScriptExecutableName);
        args.add("-u");
        args.add("-i");
        if (pythonProcessArgs != null) {
            args.addAll(pythonProcessArgs);
        }

        curatedCommandLineArgs.addAll(args);

        final InputStreamSettings isSettings = new InputStreamSettings();
        final OutputStreamSettings stdOutSettings = new OutputStreamSettings();
        stdOutSettings.setBufferSize(-1);
        final OutputStreamSettings stdErrSettings = new OutputStreamSettings();
        stdErrSettings.setBufferSize(-1);

        processSettings = new ProcessSettings(
                args.toArray(new String[args.size()]),
                false,  // redirect error
                null,           // directory
                null,         // env
                isSettings,
                stdOutSettings,
                stdErrSettings
        );

        // start the process, initialize the python code, and do the ack fifo handshake
        spController = new StreamingProcessController(processSettings, enableJournaling);
        final File ackFIFOFile = spController.start();
        if (ackFIFOFile == null) {
            return false;
        }
        initializeTool(ackFIFOFile);
        return true;
    }

    /**
     * Send a command to Python, and wait for an ack, returning all accumulated output
     * since the last call to either <link #sendSynchronousCommand/> or <line #getAccumulatedOutput/>
     * This is a blocking call - if no acknowledgment is received from the remote process, it will
     * block indefinitely. If an exception is raised in the Python code, or a negative acknowledgment
     * is received, an PythonScriptExecutorException will be thrown.
     *
     * The caller is required to terminate commands with the correct number of newline(s) as appropriate for
     * the command being issued. Since white space is significant in Python, failure to do so properly can
     * leave the Python parser blocked waiting for more newlines to terminate indented code blocks.
     *
     * @param line data to be sent to the remote process
     * @return ProcessOutput
     * @throws UserException if a timeout occurs
     */
    public ProcessOutput sendSynchronousCommand(final String line) {
        if (!line.endsWith(NL)) {
            throw new IllegalArgumentException(
                    "Python commands must be newline-terminated in order to be executed. " +
                            "Indented Python code blocks must be terminated with additional newlines");
        }
        spController.writeProcessInput(line);
        sendAckRequest();
        return waitForAck();
    }

    /**
     * Send a command to the remote process without waiting for a response. This method should only
     * be used for responses that will block the remote process.
     *
     * NOTE: Before executing further synchronous statements after calling this method, getAccumulatedOutput
     * should be called to enforce a synchronization point.
     *
     * The caller is required to terminate commands with the correct number of newline(s) as appropriate for
     * the command being issued. Since white space is significant in Python, failure to do so properly can
     * leave the Python parser blocked waiting for more newlines to terminate indented code blocks.
     *
     * @param line data to send to the remote process
     */
    public void sendAsynchronousCommand(final String line) {
        if (!line.endsWith(NL)) {
            throw new IllegalArgumentException("Python commands must be newline-terminated");
        }
        spController.writeProcessInput(line);
        sendAckRequest(); // but don't wait for it..the caller should subsequently call waitForAck
    }

    /**
     * Wait for an acknowledgement (which must have been previously requested).
     * @return true if a positive acknowledgement (ack) is received, false if negative (nck)
     */
    public ProcessOutput waitForAck() {
        if (!isAckRequestOutstanding) {
            throw new GATKException("No ack request is outstanding. An ack request must be issued first");
        }
        final boolean isAck = spController.waitForAck();
        isAckRequestOutstanding = false;
        // At every ack receipt, we want to retrieve the stdout/stderr output in case we're journaling
        final ProcessOutput po = getAccumulatedOutput();
        // if the ack was negative, throw, since the ack queue is no longer reliably in sync
        if (!isAck) {
            throw new PythonScriptExecutorException(
                    String.format(
                            "A nack was received from the Python process (most likely caused by a raised exception caused by): %s",
                            po.toString()));
        }
        return po;
    }

    /**
    /**
     * Return a (not necessarily executable) string representing the current command line for this executor
     * for error reporting purposes.
     * @return A string representing the command line used for this executor.
     */
    public String getApproximateCommandLine() {
        return curatedCommandLineArgs.stream().collect(Collectors.joining(" "));
    }

    /**
     * Obtain a stream writer that serializes and writes batches of items of type {@code T} on a background thread.
     * @param itemSerializer {@code Function} that  accepts items of type {@code T} and converts them to a
     *                                       {@code ByteArrayOutputStream} that is subsequently written to the stream
     * @return An {@link AsynchronousStreamWriter}
     */
    public void initStreamWriter(final Function<T, ByteArrayOutputStream> itemSerializer) {
        Utils.nonNull(itemSerializer, "An item serializer must be provided for the async writer service");

        dataTransferFIFOFile = spController.createDataFIFO();

        // Open the FIFO for writing. Opening a FIFO for read or write will block until there is a reader/writer
        // on the other end, so before we open it, send a non blocking, ASYNCHRONOUS command to the Python process
        // to open the FIFO for reading. The Python process will then block until we open the FIFO below.
        sendAsynchronousCommand(String.format(PYTHON_INITIALIZE_DATA_FIFO, dataTransferFIFOFile.getAbsolutePath()));
        try {
            dataTransferFIFOWriter = new FileOutputStream(dataTransferFIFOFile);
            asyncWriter = spController.getAsynchronousStreamWriter(dataTransferFIFOWriter, itemSerializer);
            // synchronize on an ack for the async command sent above before returning
            waitForAck();
        } catch ( IOException e ) {
            throw new GATKException("Failure opening FIFO for writing", e);
        }
    }

    /**
     * Request that a batch of items be written to the stream on a background thread. Any previously requested batch
     * must have already been completed and retrieved via {@link #waitForPreviousBatchCompletion}.
     *
     * @param pythonCommand command that will be executed asynchronously to cconsume the data written to the stream
     * @param batchList a list of items to be written
     */
    public void startBatchWrite(final String pythonCommand, final List<T> batchList) {
        Utils.nonNull(pythonCommand);
        Utils.nonNull(batchList);
        Utils.nonEmpty(batchList);
        sendAsynchronousCommand(pythonCommand);
        asyncWriter.startBatchWrite(batchList);
    }

    /**
     * Waits for a batch that was previously initiated via {@link #startBatchWrite(String, List)}}
     * to complete, flushes the target stream and returns the corresponding completed Future. The Future representing
     * a given batch can only be obtained via this method once. If no work is outstanding, and/or the previous batch
     * has already been retrieved, null is returned.
     * @return returns null if no previous work to complete, otherwise a completed Future
     */
    public Future<Integer> waitForPreviousBatchCompletion() {
        // wait for the batch queue to be completely written
        final Future<Integer> numberOfItemsWritten = asyncWriter.waitForPreviousBatchCompletion();
        if (numberOfItemsWritten != null) {
            // wait for the written items to be completely consumed
            waitForAck();
        }
        return numberOfItemsWritten;
    }

    /**
     * Get the Process object associated with this executor. For testing only.
     *
     * @return
     */
    @VisibleForTesting
    protected Process getProcess() {
        return spController.getProcess();
    }

    /**
     * Terminate the remote process, closing the fifo if any.
     */
    public void terminate() {
        if (profileResults != null) {
            spController.writeProcessInput(String.format(PYTHON_END_PROFILING, profileResults.getAbsolutePath()));
            sendAckRequest();
            waitForAck();
        }
        if (dataTransferFIFOWriter != null) {
            if (asyncWriter != null) {
                if(!asyncWriter.terminate()){
                    throw new GATKException("failed to close asyncWriter");
                }
            }
            spController.writeProcessInput(PYTHON_CLOSE_DATA_FIFO);
            sendAckRequest();
            waitForAck();
            try {
                dataTransferFIFOWriter.close();
                dataTransferFIFOWriter = null;
                dataTransferFIFOFile = null;
            } catch (IOException e) {
                throw new GATKException("IOException closing fifo", e);
            }
        }

        // we can't get an ack for this, since it closes down the ack fifo
        spController.writeProcessInput(PYTHON_TERMINATE_GATK);
        spController.terminate();
    }

    /**
     * Return all data accumulated since the last call to {@link #getAccumulatedOutput} (either directly, or
     * indirectly through {@link #sendSynchronousCommand}.
     *
     * Note that the output returned is somewhat non-deterministic, in that there is no guaranty that all of
     * the output from the previous command has been flushed at the time this call is made.
     *
     * @return ProcessOutput containing all accumulated output from stdout/stderr
     * @throws UserException if a timeout occurs waiting for output
     * @throws PythonScriptExecutorException if a traceback is detected in the output
     */
    public ProcessOutput getAccumulatedOutput() {
        return spController.getProcessOutput();
    }

    private void initializeTool(final File ackFIFOFile) {
        //  first we need to import the module; no ack expected yet as we haven't initialized the ack fifo
        spController.writeProcessInput(PYTHON_IMPORT_GATK); // no ack generated
        spController.writeProcessInput(String.format(PYTHON_INITIALIZE_GATK, ackFIFOFile.getAbsolutePath()));
        sendAckRequest(); // queue up an ack request
        // open the FIFO to unblock the remote caller (which should be blocked on open for read), and then
        // wait for the ack to be sent
        spController.openAckFIFOForRead();
        waitForAck();
        if (profileResults != null) {
            spController.writeProcessInput(PYTHON_START_PROFILING);
            sendAckRequest();
            waitForAck();
        }
    }

    private void sendAckRequest() {
        if (isAckRequestOutstanding) {
            throw new GATKException("An ack request is already outstanding. The previous ack request must be retrieved" +
                    " before a new ack request can be issued");
        }
        spController.writeProcessInput(PYTHON_SEND_ACK_REQUEST);
        isAckRequestOutstanding = true;
    }

}
