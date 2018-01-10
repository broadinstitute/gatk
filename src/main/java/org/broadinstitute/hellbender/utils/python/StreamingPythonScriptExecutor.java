package org.broadinstitute.hellbender.utils.python;

import com.google.common.annotations.VisibleForTesting;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.runtime.*;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeoutException;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Python executor used to interact with a keep-alive Python process. The lifecycle of an executor is typically:
 *
 *  - construct the executor
 *  - start the remote process ({@link #start}) and synchronize on the prompt ({@link #getAccumulatedOutput()}
 *  - create a fifo {@link #getFIFOForWrite}
 *  - execute asynchronous Python code {@link #sendAsynchronousCommand} to open the fifo for reading
 *  - execute local java code to open the fifo for writing
 *  - synchronize on the prompt resulting from the python code opening the fifo {@link #getAccumulatedOutput}
 *  - send/receive input one or more times (write/flush to the fifo)/synchronize on the prompt ({@link #getAccumulatedOutput) output
 *  - close the fifo locally
 *  - terminate the executor {@link #terminate}
 *
 * Guidelines for writing GATK tools that use Python interactively:
 *
 *   - Program correctness should not rely on consumption of anything written by Python to stdout/stderr other than
 *     the use the prompt for synchronization via {@link #getAccumulatedOutput}. All data should be transferred through
 *     a FIFO or file.
 *   - Always synchronize after starting the Python process (through {@link #start}, followed by
 *     {@link #getAccumulatedOutput}).
 *   - Python code should write errors to stderr.
 *   - The FIFO should always be flushed before executing Python code that reads from it. Failure to do so can result
 *     in the Python process being blocked.
 *   - Prefer single line commands that run a script, vs. multi-line Python code embedded in Java
 *   - Always terminated with newlines (otherwise Python will block)
 *   - Terminate commands with a newline.
 *   - Try not to be chatty (maximize use of the fifo buffer by writing to it in batches before reading from Python)
 *
 * NOTE: Python implementations are unreliable about honoring standard I/O stream redirection. Its not safe to
 * try to synchronize based on anything written to standard I/O streams, since Python sometimes prints the
 * prompt to stdout, and sometimes to stderr:
 *
 *  https://bugs.python.org/issue17620
 *  https://bugs.python.org/issue1927
 */
public class StreamingPythonScriptExecutor extends PythonExecutorBase {
    private static final Logger logger = LogManager.getLogger(StreamingPythonScriptExecutor.class);
    private static final String NL = System.lineSeparator();

    private final List<String> curatedCommandLineArgs = new ArrayList<>();

    private StreamingProcessController spController;
    private ProcessSettings processSettings;

    final public static String PYTHON_PROMPT = ">>> ";

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
        return start(pythonProcessArgs, false);
    }

    /**
     * Start the Python process.
     *
     * @param pythonProcessArgs args to be passed to the python process
     * @param enableJournaling true to enable Journaling, which records all interprocess IO to a file. This is
     *                         expensive and should only be used for debugging purposes.
     * @return true if the process is successfully started
     */
    public boolean start(final List<String> pythonProcessArgs, final boolean enableJournaling) {
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

        spController = new StreamingProcessController(processSettings, PYTHON_PROMPT, enableJournaling);
        return spController.start();
    }

    /**
     * Send a command to Python, and wait for a response prompt, returning all accumulated output
     * since the last call to either <link #sendSynchronousCommand/> or <line #getAccumulatedOutput/>
     * This is a blocking call, and should be used for commands that execute quickly and synchronously.
     * If no output is received from the remote process during the timeout period, an exception will be thrown.
     *
     * The caller is required to terminate commands with a newline. The executor doesn't do this
     * automatically since doing so would alter the number of prompts issued by the remote process.
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
        return getAccumulatedOutput();
    }

    /**
     * Send a command to the remote process without waiting for a response. This method should only
     * be used for responses that will block the remote process.
     *
     * NOTE: Before executing further synchronous statements after calling this method, getAccumulatedOutput
     * should be called to enforce a synchronization point.
     *
     * The caller is required to terminate commands with a newline. The executor doesn't do this
     * automatically since it can alter the number of prompts, and thus synchronization points, issued
     * by the remote process.
     *
     * @param line data to send to the remote process
     */
    public void sendAsynchronousCommand(final String line) {
        if (!line.endsWith(NL)) {
            throw new IllegalArgumentException("Python commands must be newline-terminated");
        }
        spController.writeProcessInput(line);
    }

    /**
     * See if any output is currently available. This is non-blocking, and can be used to determine if a blocking
     * call can be made; it is always safe to call getAccumulatedOutput if isOutputAvailable is true.
     * @return true if data is available from the remote process.
     */
    public boolean isOutputAvailable() {
        return spController.isOutputAvailable();
    }

    /**
     * Return all data accumulated since the last call to {@link #getAccumulatedOutput} (either directly, or
     * indirectly through {@link #sendSynchronousCommand}, collected until an output prompt is detected.
     *
     * Note that the output returned is somewhat non-deterministic, in that the only guaranty is that a prompt
     * was detected on either stdout or stderr. It is possible for the remote process to produce the prompt on
     * one stream (stderr or stdout), and additional output on the other; this method may detect the prompt before
     * detecting the additional output on the other stream. Such output will be retained, and returned as part of
     * the payload the next time output is retrieved.
     *
     * For this reason, program correctness should not rely on consuming data written by Python to standard streams.
     *
     * This should only be used for short, synchronous commands that produce output quickly. If no data has been
     * sent from the process, this call blocks waiting for data, or the timeout (default 5 seconds) to be reached.
     *
     * Longer-running, blocking commands can be executed using {@link #sendAsynchronousCommand}, in combination
     * with {@link #isOutputAvailable}.
     *
     * @return ProcessOutput containing all accumulated output from stdout/stderr
     * @throws UserException if a timeout occurs waiting for output
     * @throws PythonScriptExecutorException if a traceback is detected in the output
     */
    public ProcessOutput getAccumulatedOutput() {
        try {
            final ProcessOutput po = spController.getProcessOutputByPrompt();
            final StreamOutput stdErr = po.getStderr();
            if (stdErr != null) {
                final String stdErrText = stdErr.getBufferString();
                if (stdErrText != null && stdErrText.contains("Traceback")) {
                    throw new PythonScriptExecutorException("Traceback detected: " + stdErrText);
                }
            }
            return po;
        } catch (TimeoutException e) {
            throw new UserException("A timeout ocurred waiting for output from the remote Python command.", e);
        }
    }

    /**
     * Terminate the remote process, closing the fifo if any.
     */
    public void terminate() {
        spController.terminate();
    }

    /**
     * Obtain a temporary FIFO to be used to transfer data to Python. The FIFO is only valid for the
     * lifetime of the executor; it is destroyed when the executor is terminated.
     *
     * NOTE: Since opening a FIFO for write blocks until it is opened for read, the caller is responsible
     * for ensuring that a Python command to open the FIFO has been executed (asynchronously) before executing
     * code to open it for write. For this reason, the opening of the FIFO is left to the caller.
     *
     * @return
     */
    public File getFIFOForWrite() {
        return spController.createFIFO();
    }

    /**
     * Return a {@link AsynchronousStreamWriterService} to be used to write to an output stream, typically on a FIFO,
     * on a background thread.
     * @param streamWriter stream to which items should be written.
     * @param itemSerializer function that converts an item of type {@code T} to a {@code ByteArrayOutputStream} for serialization
     * @param <T> Type of items to be written to the stream.
     * @return {@link AsynchronousStreamWriterService}
     */
    public <T> AsynchronousStreamWriterService<T> getAsynchronousStreamWriterService(
            final OutputStream streamWriter,
            final Function<T, ByteArrayOutputStream> itemSerializer)
    {
        Utils.nonNull(streamWriter);
        Utils.nonNull(itemSerializer);

        return spController.getAsynchronousStreamWriterService(streamWriter, itemSerializer);
    }

    /**
     * Return a (not necessarily executable) string representing the current command line for this executor
     * for error reporting purposes.
     * @return A string representing the command line used for this executor.
     */
    public String getApproximateCommandLine() {
        return curatedCommandLineArgs.stream().collect(Collectors.joining(" "));
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
}
