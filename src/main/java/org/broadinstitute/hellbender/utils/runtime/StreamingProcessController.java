package org.broadinstitute.hellbender.utils.runtime;

import com.google.common.io.Files;
import org.apache.commons.io.IOUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.*;
import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;

/**
 * Facade to Runtime.exec() and java.lang.Process. Handles running of an interactive, keep-alive process that
 * uses a FIFO to indicate remote command success or failure. The controller is agnostic about what process is
 * being started; the only requirement is that it is capable of executing, at the callers demand, a command that
 * writes an acknowledgement message to the ack FIFO file that is managed by this controller, which can be
 * detected by {@link #waitForAck}. Creates separate threads for reading stdout and stderr and provides access
 * to stdout and stderr as strings.
 */
public final class StreamingProcessController extends ProcessControllerBase<CapturedStreamOutputSnapshot> {
    private static final Logger logger = LogManager.getLogger(StreamingProcessController.class);

    private final ProcessSettings settings;

    // Names for the ack and data fifo files used by this controller. Each controller instance creates a
    // unique temporary directory for these files, so the names don't have to be unique across instances.
    private static String ACK_FIFO_FILE_NAME = "gatkStreamingControllerAck.fifo";
    private static String DATA_FIFO_FILE_NAME = "gatkStreamingControllerData.fifo";

    private File fifoTempDir = null;
    private File ackFIFOFile = null;
    private File dataFIFOFile = null;

    // These strings must be kept in sync with the ones used in the Python package
    public static String ACK_MESSAGE = "ack";
    public static String NCK_MESSAGE = "nck";
    private static int ACK_MESSAGE_SIZE = 3; // "ack" or "nck"
    private static String ACK_LOG_MESSAGE = "Ack received\n\n";
    private static String NCK_LOG_MESSAGE = "Nck received\n\n";
    private InputStream ackFIFOInputStream;
    private Future<Boolean> ackFuture;

    // keep an optional journal of all IPC; disabled/no-op by default
    private ProcessJournal processJournal = new ProcessJournal();
    private OutputStream processStdinStream; // stream to which we send remote commands

    //Timeout used when terminating the remote process
    private static final int REMOTE_PROCESS_TERMINATION_TIMEOUT_SECONDS = 30;

    /**
     * @param settings Settings to be used.
     */
    public StreamingProcessController(final ProcessSettings settings) {
        this(settings, false);
    }

    /**
     * Create a controller using the specified settings.
     *
     * @param settings                 Settings to be run.
     * @param enableJournaling         Turn on process I/O journaling. This records all inter-process communication to a file.
     *                                 Journaling incurs a performance penalty, and should be used for debugging purposes only.
     */
    public StreamingProcessController(
            final ProcessSettings settings,
            final boolean enableJournaling) {
        Utils.nonNull(settings, "Process settings are required");
        this.settings = settings;

        if (enableJournaling) {
            processJournal.enable(settings.getCommandString());
        }
    }

    /**
     * Starts the remote process running based on the setting specified in the constructor.
     *
     * @return the FIFO File to be used for ack messages. This caller should pass the name
     * to the companion process to be used for writing ack/nck messages.
     */
    public File start() {
        if (process != null) {
            throw new IllegalStateException("This controller is already running a process");
        }
        process = launchProcess(settings);
        startListeners();
        processStdinStream = getProcess().getOutputStream();

        // create the fifo temp directory, and one FIFO to use for IPC signalling
        fifoTempDir = Files.createTempDir();
        ackFIFOFile = createFIFOFile(ACK_FIFO_FILE_NAME);
        return ackFIFOFile;
    }

    /**
     * Write some input to the remote process, without waiting for output.
     *
     * @param line data to be written to the remote process
     */
    public void writeProcessInput(final String line) {
        try {
            // Its possible that we already have completed futures from output from a previous command that
            // weren't consumed before the last synchronization point. If so, drop them and write any existing
            // output to the journal so it appears in the journal before the command we're about to send.
            if (stdErrFuture != null && stdErrFuture.isDone()) {
                processJournal.writeLogMessage("Dropping stale stderr output before send: \n" + stdErrFuture.get().getBufferString() + "\n");
                stdErrFuture = null;
            }
            if (stdOutFuture != null && stdOutFuture.isDone()) {
                processJournal.writeLogMessage("Dropping stale stdout output before send: \n" + stdOutFuture.get().getBufferString() + "\n");
                stdOutFuture = null;
            }
        } catch (InterruptedException e) {
            throw new GATKException(String.format("Interrupted retrieving stale future: " + line, e));
        } catch (ExecutionException e) {
            throw new GATKException(String.format("Execution exception retrieving stale future: " + line, e));
        }

        startListeners();

        try {
            // write and flush the output stream that is the process' input
            processStdinStream.write(line.getBytes());
            processStdinStream.flush();
            processJournal.writeOutbound(line);
        } catch (IOException e) {
            throw new GATKException(String.format("Error writing (%s) to stdin on command", line), e);
        }
    }

    /**
     * Attempt to retrieve any output from a (possibly) completed future that is immediately available without
     * blocking or waiting for the future to complete.
     * @param streamOutputFuture Future from which to retrieve output
     * @return {@link StreamOutput} if the future is complete, or null otherwise
     */
    private StreamOutput drainOutputStream(final Future<CapturedStreamOutputSnapshot> streamOutputFuture) {
        StreamOutput ret = null;

        if (streamOutputFuture != null) {
            try {
                if (streamOutputFuture.isDone()) {
                    ret = streamOutputFuture.get();
                }
            } catch (InterruptedException e) {
                throw new GATKException("InterruptedException attempting to retrieve output from remote process", e);
            } catch (ExecutionException e) {
                throw new GATKException("ExecutionException attempting to retrieve output from remote process", e);
            }
        }
        return ret;
    }

    /**
     * Wait for *any* output from the process (to stdout or stderr), by draining the stream(s) until no more
     * output is available. Uses a timeout to prevent the calling thread from hanging. Use isOutputAvailable
     * for non-blocking check.
     */
    public ProcessOutput getProcessOutput() {
        StreamOutput stdout;
        StreamOutput stderr;

        // Get whatever output is immediately available without blocking.
        stderr = drainOutputStream(stdErrFuture);
        if (stderr != null) {
            stdErrFuture = null;
        }
        stdout = drainOutputStream(stdOutFuture);
        if (stdout != null) {
            stdOutFuture = null;
        }

        // log the output, and restart the listeners for next time
        processJournal.writeInbound(stdout, stderr);
        startListeners();

        return new ProcessOutput(0, stdout, stderr);
    }

    /**
     * Open a stream on the ack FIFO for reading.
     *
     * NOTE: Before this is called, the caller must have requested execution of code in the remote process to
     * open the ackFIFOFile for write. Failure to do so would cause this call to block indefinitely.
     */
    public void openAckFIFOForRead() {
        try {
            ackFIFOInputStream = new FileInputStream(ackFIFOFile);
        } catch (FileNotFoundException e) {
            throw new GATKException("Can't open ack FIFO for read");
        }
    }

    /**
     * Wait for a previously requested acknowledgement to be received. The remote process can deliver a positive
     * ack to indicate successful command completion, or a negative ack to indicate command execution failure.
     * @return true if an positive acknowledgement (ACK) was received, false if a negative acknowledgement (NCK)
     * was received
     */
    public boolean waitForAck() {
        if (ackFuture != null) {
            throw new GATKException("An ack is already outstanding");
        }
        ackFuture = executorService.submit(
                () -> {
                    try {
                        byte[] ack = new byte[ACK_MESSAGE_SIZE];
                        int nBytes = ackFIFOInputStream.read(ack, 0, ACK_MESSAGE_SIZE);
                        if (nBytes != ACK_MESSAGE_SIZE) {
                            throw new GATKException(String.format("Failure reading ack message from ack fifo, ret: (%d)", nBytes));
                        } else if (Arrays.equals(ack, NCK_MESSAGE.getBytes())) {
                            return false;
                        } else if (Arrays.equals(ack, ACK_MESSAGE.getBytes())) {
                            return true;
                        } else {
                            logger.error("Unrecognized string written to ack fifo");
                            return false;
                        }
                    } catch (IOException e) {
                        throw new GATKException("IOException reading from ack fifo", e);
                    }
                }
        );

        try {
            // blocking call to wait for the ack
            boolean isAck = ackFuture.get();
            processJournal.writeLogMessage(isAck ? ACK_LOG_MESSAGE : NCK_LOG_MESSAGE);
            ackFuture = null;
            return isAck;
        } catch (InterruptedException | ExecutionException e) {
            throw new GATKException("Exception waiting for ack from Python: " + e.getMessage(), e);
        }
    }

    /**
     * Create a temporary FIFO suitable for sending output to the remote process. The FIFO is only valid for the
     * lifetime of the controller; the FIFO is destroyed when the controller is destroyed.
     *
     * @return a FIFO File
     */
    public File createDataFIFO() {
        if (dataFIFOFile != null) {
            throw new IllegalArgumentException("Only one data FIFO per controller is supported");
        }
        dataFIFOFile = createFIFOFile(DATA_FIFO_FILE_NAME);
        return dataFIFOFile;
    }

    private File createFIFOFile(final String fifoName) {
        final String fifoTempFileName = String.format("%s/%s", fifoTempDir.getAbsolutePath(), fifoName);
        return org.broadinstitute.hellbender.utils.io.IOUtils.createFifoFile(org.broadinstitute.hellbender.utils.io.IOUtils.getPath(fifoTempFileName), true);
    }

    /**
     * Return a {@link AsynchronousStreamWriter} to be used to write to a stream on a background thread.
     * @param outputStream stream to which items should be written.
     * @param itemSerializer
     * @param <T> Type of items to be written to the stream.
     * @return {@link AsynchronousStreamWriter}
     */
    public <T> AsynchronousStreamWriter<T> getAsynchronousStreamWriter(
            final OutputStream outputStream,
            final Function<T, ByteArrayOutputStream> itemSerializer) {
        Utils.nonNull(outputStream);
        Utils.nonNull(itemSerializer);
        return new AsynchronousStreamWriter<>(executorService, outputStream, itemSerializer);
    }

    /**
     * Close the FIFO; called on controller termination
     */
    private void closeFIFOs() {
        if (dataFIFOFile != null) {
            dataFIFOFile.delete();
        }
        if (ackFIFOFile != null) {
            ackFIFOFile.delete();
        }
        fifoTempDir.delete();
    }

    /**
     * Submit a task to start listening for output
     */
    private void startListeners() {
        // Submit runnable task to start capturing.
        if (stdOutFuture == null) {
            stdOutFuture = executorService.submit(
                    new OutputCapture(
                            new CapturedStreamOutputSnapshot(settings.getStdoutSettings(), process.getInputStream(), System.out),
                            ProcessStream.STDOUT,
                            this.getClass().getSimpleName(),
                            controllerId));
        }
        if (!settings.isRedirectErrorStream() && stdErrFuture == null) {
            // don't waste a callable on stderr if its just going to be redirected to stdout
            stdErrFuture = executorService.submit(
                    new OutputCapture(
                            new CapturedStreamOutputSnapshot(settings.getStderrSettings(), process.getErrorStream(), System.err),
                            ProcessStream.STDERR,
                            this.getClass().getSimpleName(),
                            controllerId));
        }
    }

    /**
     * Stops the process from running and tries to ensure the process is cleaned up properly.
     * NOTE: sub-processes started by process may be zombied with their parents set to pid 1.
     * NOTE: capture threads may block on read.
     */
    @Override
    protected void tryCleanShutdown() {
        if (stdErrFuture != null && !stdErrFuture.isDone()) {
            boolean isCancelled = stdErrFuture.cancel(true);
            if (!isCancelled) {
                logger.error("Failure cancelling stderr task");
            }
        }
        if (stdOutFuture != null && !stdOutFuture.isDone()) {
            boolean isCancelled = stdOutFuture.cancel(true);
            if (!isCancelled) {
                logger.error("Failure cancelling stdout task");
            }
        }
        if (process != null) {
            // terminate the app by closing the process' INPUT stream
            IOUtils.closeQuietly(process.getOutputStream());
        }
    }

    /**
     * Close the input stream, close the FIFO, and wait for the remote process to terminate
     * destroying it if necessary.
     */
    public void terminate() {
        closeFIFOs();
        tryCleanShutdown();
        boolean exited = false;
        try {
            exited = process.waitFor(REMOTE_PROCESS_TERMINATION_TIMEOUT_SECONDS, TimeUnit.SECONDS);
            processJournal.close();
            if (!exited) {
                // we timed out waiting for the process to exit; it may be in a blocking call, so just force
                // it to shutdown
                process.destroy();
                process.waitFor(REMOTE_PROCESS_TERMINATION_TIMEOUT_SECONDS, TimeUnit.SECONDS);
            }
        } catch (InterruptedException e) {
            logger.error(String.format("Interrupt exception waiting for process (%s) to terminate", settings.getCommandString()));
        }
        if (process.isAlive()) {
            throw new GATKException("Failure terminating remote process");
        }
    }

    /**
     * Keep a journal of all inter-process communication. Can incur significant runtime overhead, and should
     * be used for debugging only.
     */
    private class ProcessJournal {
        private File journalingFile = null;
        private FileWriter journalingFileWriter;

        public void enable(final String commandString) {
            final String journalingFileName = String.format("gatkStreamingProcessJournal-%d.txt", new Random().nextInt());
            journalingFile = new File(journalingFileName);
            try {
                journalingFileWriter = new FileWriter(journalingFile);
                journalingFileWriter.write("Initial process command line: ");
                journalingFileWriter.write(settings.getCommandString() + "\n\n");
            } catch (IOException e) {
                throw new GATKException(String.format("Error creating streaming process journaling file %s for command \"%s\"",
                        commandString,
                        journalingFile.getAbsolutePath()), e);
            }
            logger.info(String.format("Enabling streaming process journaling file %s", journalingFileName));
        }

        /**
         * Record outbound data being written to a remote process.
         * @param line line being written
         */
        public void writeOutbound(final String line) {
            try {
                if (journalingFileWriter != null) {
                    journalingFileWriter.write("Sending: \n[");
                    journalingFileWriter.write(line);
                    journalingFileWriter.write("]\n\n");
                    journalingFileWriter.flush();
                }
            } catch (IOException e) {
                throw new GATKException("Error writing to output to process journal", e);
            }
        }

        /**
         * Record inbound data being received from a remote process.
         * @param stdout Data received from stdout.
         * @param stderr Data received from stderr.
         */
        public void writeInbound(final StreamOutput stdout, final StreamOutput stderr) {

            if (journalingFileWriter != null) {
                try {
                    if (stdout != null) {
                        journalingFileWriter.write("Received from stdout: [");
                        journalingFileWriter.write(stdout.getBufferString());
                        journalingFileWriter.write("]\n");
                    }
                    if (stderr != null) {
                        journalingFileWriter.write("Received from stderr: [");
                        journalingFileWriter.write(stderr.getBufferString());
                        journalingFileWriter.write("]\n");
                    }
                    journalingFileWriter.write("\n");
                    journalingFileWriter.flush();
                } catch (IOException e) {
                    throw new GATKException(String.format("Error writing to journaling file %s", journalingFile.getAbsolutePath()), e);
                }
            }
        }

        /**
         * Write an advisory log message to the journal.
         * @param message Message to be written to the journal.
         */
        public void writeLogMessage(final String message) {
            if (journalingFileWriter != null) {
                try {
                    journalingFileWriter.write(message);
                    journalingFileWriter.flush();
                } catch (IOException e) {
                    throw new GATKException(String.format("Error writing to journaling file %s", journalingFile.getAbsolutePath()), e);
                }
            }
        }

        /**
         * Close the journal file writer.
         */
        public void close() {
            try {
                if (journalingFileWriter != null) {
                    writeLogMessage("Shutting down journal normally");
                    journalingFileWriter.flush();
                    journalingFileWriter.close();
                }
            } catch (IOException e) {
                throw new GATKException(String.format("Error closing streaming process journaling file %s", journalingFile.getAbsolutePath()), e);
            }
        }
    }

}
