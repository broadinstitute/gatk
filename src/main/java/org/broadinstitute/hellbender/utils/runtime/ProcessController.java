package org.broadinstitute.hellbender.utils.runtime;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.*;
import java.util.concurrent.ExecutionException;

/**
 * Facade to Runtime.exec() and java.lang.Process.  Handles
 * running a process to completion and returns stdout and stderr
 * as strings.  Creates separate threads for reading stdout and stderr,
 * then reuses those threads for each process most efficient use is
 * to create one of these and use it repeatedly.  Instances are not
 * thread-safe, however.
 *
 * TODO: java.io sometimes zombies the backround threads locking up on read().
 * Supposedly NIO has better ways of interrupting a blocked stream but will
 * require a little bit of refactoring.
 */
public final class ProcessController extends ProcessControllerBase<CapturedStreamOutput> {
    private static final Logger logger = LogManager.getLogger(ProcessController.class);

    /**
     * Returns a thread local ProcessController.
     * Should NOT be closed when finished so it can be reused by the thread.
     *
     * @return a thread local ProcessController.
     */
    public static ProcessController getThreadLocal() {
        // If the local controller was destroyed get a fresh instance.
        if (threadProcessController.get().destroyed)
            threadProcessController.remove();
        return threadProcessController.get();
    }

    /**
     * Thread local process controller container.
     */
    private static final ThreadLocal<ProcessController> threadProcessController =
            new ThreadLocal<ProcessController>() {
                @Override
                protected ProcessController initialValue() {
                    return new ProcessController();
                }
            };

    /**
     * Similar to Runtime.exec() but drains the output and error streams.
     *
     * @param command Command to run.
     * @return The result code.
     */
    public static int exec(String[] command) {
        ProcessController controller = ProcessController.getThreadLocal();
        return controller.exec(new ProcessSettings(command)).getExitValue();
    }

    /**
     * Executes a command line program with the settings and waits for it to return,
     * processing the output on a background thread.
     *
     * @param settings Settings to be run.
     * @return The output of the command.
     */
    public ProcessOutput exec(ProcessSettings settings) {
        StreamOutput stdout;
        StreamOutput stderr;
        int exitCode;

        launchProcess(settings);

        try {
            stdOutFuture = executorService.submit(
                    new OutputCapture(
                            new CapturedStreamOutput(settings.getStdoutSettings(), process.getInputStream(), System.out),
                            ProcessStream.STDOUT,
                            this.getClass().getSimpleName(),
                            controllerId));
            stdErrFuture = executorService.submit(
                    new OutputCapture(
                            new CapturedStreamOutput(settings.getStderrSettings(), process.getErrorStream(), System.err),
                            ProcessStream.STDERR,
                            this.getClass().getSimpleName(),
                            controllerId));

            // Write stdin content
            InputStreamSettings stdinSettings = settings.getStdinSettings();
            Set<StreamLocation> streamLocations = stdinSettings.getStreamLocations();
            if (!streamLocations.isEmpty()) {
                try {
                    OutputStream stdinStream = process.getOutputStream();
                    for (StreamLocation location : streamLocations) {
                        InputStream inputStream;
                        switch (location) {
                            case Buffer:
                                inputStream = new ByteArrayInputStream(stdinSettings.getInputBuffer());
                                break;
                            case File:
                                try {
                                    inputStream = FileUtils.openInputStream(stdinSettings.getInputFile());
                                } catch (IOException e) {
                                    throw new UserException.BadInput(e.getMessage());
                                }
                                break;
                            case Standard:
                                inputStream = System.in;
                                break;
                            default:
                                throw new GATKException("Unexpected stream location: " + location);
                        }
                        try {
                            IOUtils.copy(inputStream, stdinStream);
                        } finally {
                            if (location != StreamLocation.Standard)
                                IOUtils.closeQuietly(inputStream);
                        }
                    }
                    stdinStream.flush();
                } catch (IOException e) {
                    throw new GATKException("Error writing to stdin on command: " + settings.getCommandString(), e);
                }
            }

            // Wait for the process to complete.
            try {
                process.getOutputStream().close();
                process.waitFor();
                stdout = stdOutFuture.get();
                stdOutFuture = null;
                stderr = stdErrFuture.get();
                stdErrFuture = null;
            } catch (ExecutionException e) {
                throw new GATKException("Execution exception during process output retrieval", e);
            } catch (IOException e) {
                throw new GATKException("Unable to close stdin on command: " + settings.getCommandString(), e);
            } catch (InterruptedException e) {
                throw new GATKException("Process interrupted", e);
            }
        } finally {
            exitCode = process.exitValue();
            process = null;
            running.remove(this);
        }
        return new ProcessOutput(exitCode, stdout, stderr);

    }

    /**
     * Stops the process from running and tries to ensure process is cleaned up properly.
     * NOTE: sub-processes started by process may be zombied with their parents set to pid 1.
     * NOTE: capture threads may block on read.
     * TODO: Try to use NIO to interrupt streams.
     */
    @Override
    protected void tryCleanShutdown() {
        destroyed = true;

        if (stdErrFuture != null) {
            boolean isCancelled = stdErrFuture.cancel(true);
            if (!isCancelled) {
                logger.error("Failure cancelling stderr task");
            }
        }
        if (stdOutFuture != null) {
            boolean isCancelled = stdOutFuture.cancel(true);
            if (!isCancelled) {
                logger.error("Failure cancelling stdout task");
            }
        }
        if (process != null) {
            process.destroy();
            IOUtils.closeQuietly(process.getInputStream());
            IOUtils.closeQuietly(process.getErrorStream());
        }
    }

}
