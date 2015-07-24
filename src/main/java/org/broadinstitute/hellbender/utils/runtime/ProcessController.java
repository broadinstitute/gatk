package org.broadinstitute.hellbender.utils.runtime;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.*;

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
public final class ProcessController {
    private static final Logger logger = LogManager.getLogger(ProcessController.class);

    private static enum ProcessStream {Stdout, Stderr}

    // Tracks running processes.
    private static final Set<ProcessController> running = Collections.synchronizedSet(new HashSet<>());

    // Tracks this running process.
    private Process process;

    // Threads that capture stdout and stderr
    private final OutputCapture stdoutCapture;
    private final OutputCapture stderrCapture;

    // When a caller destroys a controller a new thread local version will be created
    private boolean destroyed = false;

    // Communication channels with output capture threads

    // Holds the stdout and stderr sent to the background capture threads
    private final Map<ProcessStream, CapturedStreamOutput> toCapture =
            new EnumMap<>(ProcessStream.class);

    // Holds the results of the capture from the background capture threads.
    // May be the content via toCapture or an StreamOutput.EMPTY if the capture was interrupted.
    private final Map<ProcessStream, StreamOutput> fromCapture =
            new EnumMap<>(ProcessStream.class);

    // Useful for debugging if background threads have shut down correctly
    private static int nextControllerId = 0;
    private final int controllerId;

    public ProcessController() {
        // Start the background threads for this controller.
        synchronized (running) {
            controllerId = nextControllerId++;
        }
        stdoutCapture = new OutputCapture(ProcessStream.Stdout, controllerId);
        stderrCapture = new OutputCapture(ProcessStream.Stderr, controllerId);
        stdoutCapture.start();
        stderrCapture.start();
    }

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
        if (destroyed)
            throw new IllegalStateException("This controller was destroyed");

        ProcessBuilder builder = new ProcessBuilder(settings.getCommand());
        builder.directory(settings.getDirectory());

        Map<String, String> settingsEnvironment = settings.getEnvironment();
        if (settingsEnvironment != null) {
            Map<String, String> builderEnvironment = builder.environment();
            builderEnvironment.clear();
            builderEnvironment.putAll(settingsEnvironment);
        }

        builder.redirectErrorStream(settings.isRedirectErrorStream());

        StreamOutput stdout = null;
        StreamOutput stderr = null;

        // Start the process running.

        try {
            synchronized (toCapture) {
                process = builder.start();
            }
            running.add(this);
        } catch (IOException e) {
            String message = String.format("Unable to start command: %s\nReason: %s",
                    StringUtils.join(builder.command(), " "),
                    e.getMessage());
            throw new GATKException(message);
        }

        int exitCode;

        try {
            // Notify the background threads to start capturing.
            synchronized (toCapture) {
                toCapture.put(ProcessStream.Stdout,
                        new CapturedStreamOutput(settings.getStdoutSettings(), process.getInputStream(), System.out));
                toCapture.put(ProcessStream.Stderr,
                        new CapturedStreamOutput(settings.getStderrSettings(), process.getErrorStream(), System.err));
                toCapture.notifyAll();
            }

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
                    throw new GATKException("Error writing to stdin on command: " + StringUtils.join(builder.command(), " "), e);
                }
            }

            // Wait for the process to complete.
            try {
                process.getOutputStream().close();
                process.waitFor();
            } catch (IOException e) {
                throw new GATKException("Unable to close stdin on command: " + StringUtils.join(builder.command(), " "), e);
            } catch (InterruptedException e) {
                throw new GATKException("Process interrupted", e);
            } finally {
                while (!destroyed && stdout == null || stderr == null) {
                    synchronized (fromCapture) {
                        if (fromCapture.containsKey(ProcessStream.Stdout))
                            stdout = fromCapture.remove(ProcessStream.Stdout);
                        if (fromCapture.containsKey(ProcessStream.Stderr))
                            stderr = fromCapture.remove(ProcessStream.Stderr);
                        try {
                            if (stdout == null || stderr == null)
                                fromCapture.wait();
                        } catch (InterruptedException e) {
                            // Log the error, ignore the interrupt and wait patiently
                            // for the OutputCaptures to (via finally) return their
                            // stdout and stderr.
                            logger.error(e);
                        }
                    }
                }

                if (destroyed) {
                    if (stdout == null)
                        stdout = StreamOutput.EMPTY;
                    if (stderr == null)
                        stderr = StreamOutput.EMPTY;
                }
            }
        } finally {
            synchronized (toCapture) {
                exitCode = process.exitValue();
                process = null;
            }
            running.remove(this);
        }

        return new ProcessOutput(exitCode, stdout, stderr);
    }

    /**
     * @return The set of still running processes.
     */
    public static Set<ProcessController> getRunning() {
        synchronized (running) {
            return new HashSet<>(running);
        }
    }

    /**
     * Stops the process from running and tries to ensure process is cleaned up properly.
     * NOTE: sub-processes started by process may be zombied with their parents set to pid 1.
     * NOTE: capture threads may block on read.
     * TODO: Try to use NIO to interrupt streams.
     */
    public void tryDestroy() {
        destroyed = true;
        synchronized (toCapture) {
            if (process != null) {
                process.destroy();
                IOUtils.closeQuietly(process.getInputStream());
                IOUtils.closeQuietly(process.getErrorStream());
            }
            stdoutCapture.interrupt();
            stderrCapture.interrupt();
            toCapture.notifyAll();
        }
    }

    @Override
    protected void finalize() throws Throwable {
        try {
            tryDestroy();
        } catch (Exception e) {
            logger.error(e);
        }
        super.finalize();
    }

    private class OutputCapture extends Thread {
        private final int controllerId;
        private final ProcessStream key;

        /**
         * Reads in the output of a stream on a background thread to keep the output pipe from backing up and freezing the called process.
         *
         * @param key The stdout or stderr key for this output capture.
         * @param controllerId Unique id of the controller.
         */
        public OutputCapture(ProcessStream key, int controllerId) {
            super(String.format("OutputCapture-%d-%s-%s-%d", controllerId, key.name().toLowerCase(),
                    Thread.currentThread().getName(), Thread.currentThread().getId()));
            this.controllerId = controllerId;
            this.key = key;
            setDaemon(true);
        }

        /**
         * Runs the capture.
         */
        @Override
        public void run() {
            while (!destroyed) {
                StreamOutput processStream = StreamOutput.EMPTY;
                try {
                    // Wait for a new input stream to be passed from this process controller.
                    CapturedStreamOutput capturedProcessStream = null;
                    while (!destroyed && capturedProcessStream == null) {
                        synchronized (toCapture) {
                            if (toCapture.containsKey(key)) {
                                capturedProcessStream = toCapture.remove(key);
                            } else {
                                toCapture.wait();
                            }
                        }
                    }

                    if (!destroyed) {
                        // Read in the input stream
                        processStream = capturedProcessStream;
                        capturedProcessStream.readAndClose();
                    }
                } catch (InterruptedException e) {
                    logger.info("OutputCapture interrupted, exiting");
                    break;
                } catch (IOException e) {
                    logger.error("Error reading process output", e);
                } finally {
                    // Send the string back to the process controller.
                    synchronized (fromCapture) {
                        fromCapture.put(key, processStream);
                        fromCapture.notify();
                    }
                }
            }
        }
    }
}
