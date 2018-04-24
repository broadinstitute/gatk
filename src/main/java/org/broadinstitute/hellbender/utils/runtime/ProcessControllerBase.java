package org.broadinstitute.hellbender.utils.runtime;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.*;

public abstract class ProcessControllerBase<CAPTURE_POLICY extends CapturedStreamOutput>
{
    private static final Logger logger = LogManager.getLogger(ProcessControllerBase.class);

    protected enum ProcessStream {STDOUT, STDERR}

    // Tracks running controllers.
    protected static final Set<ProcessControllerBase<?>> running = Collections.synchronizedSet(new LinkedHashSet<>());

    // We can't use a fixed thread pool size, since we could run out if multiple controllers were to make blocking calls
    final static protected ExecutorService executorService = Executors.newCachedThreadPool(
            new ThreadFactory() {
                @Override
                public Thread newThread(Runnable r) {
                    final Thread newThread = Executors.defaultThreadFactory().newThread(r);
                    newThread.setDaemon(true);
                    newThread.setName("GATKProcessController");
                    return newThread;
                }
            }
    );

    // Tracks the running process associated with this controller.
    protected Process process;

    // The capture policy dictates the termination condition when reading output data from the remote process,
    // either read everything (for ProcessController) or read while available (StreamingProcessController).
    protected Future<CAPTURE_POLICY> stdOutFuture;
    protected Future<CAPTURE_POLICY> stdErrFuture;

    // When a caller destroys a controller a new thread local version will be created
    protected boolean destroyed = false;

    // Useful for debugging if background threads have shut down correctly
    private static int nextControllerId = 0;
    protected final int controllerId;

    public ProcessControllerBase() {
        synchronized (running) {
            controllerId = nextControllerId++;
        }
    }

    /**
     * @return The set of still running processes.
     */
    public static Set<ProcessControllerBase<?>> getRunning() {
        synchronized (running) {
            return new LinkedHashSet<>(running);
        }
    }

    /**
     * Executes a command line program with the settings and waits for it to return,
     * processing the output on a background thread.
     *
     * @param settings Settings to be run.
     * @return The output of the command.
     */
    protected Process launchProcess(ProcessSettings settings) {
        Utils.validate(!destroyed, "This controller was destroyed");

        final ProcessBuilder builder = new ProcessBuilder(settings.getCommand());
        builder.directory(settings.getDirectory());

        final Map<String, String> settingsEnvironment = settings.getEnvironment();
        if (settingsEnvironment != null) {
            final Map<String, String> builderEnvironment = builder.environment();
            builderEnvironment.clear();
            builderEnvironment.putAll(settingsEnvironment);
        }

        builder.redirectErrorStream(settings.isRedirectErrorStream());

        // Start the process running.
        try {
            process = builder.start();
            running.add(this);
        } catch (IOException e) {
            final String message = String.format("Unable to start command: %s\nReason: %s",
                    StringUtils.join(builder.command(), " "),
                    e.getMessage());
            throw new GATKException(message, e);
        }
        return process;
    }

    @VisibleForTesting
    public Process getProcess() {
        return process;
    }

    /**
     * Stops the process from running and tries to ensure process is cleaned up properly.
     * NOTE: sub-processes started by process may be zombied with their parents set to pid 1.
     * NOTE: capture threads may block on read.
     * TODO: Try to use NIO to interrupt streams.
     */
    abstract void tryCleanShutdown();

    @Override
    protected void finalize() throws Throwable {
        try {
            tryCleanShutdown();
        } catch (Exception e) {
            logger.error(e);
        }
        super.finalize();
    }

    protected class OutputCapture implements Callable<CAPTURE_POLICY> {
        private final CAPTURE_POLICY capturedProcessStream;
        private final ProcessStream key; // un-referenced, but retained for debugging
        private final String controllerContextName;  // un-referenced, but retained for debugging

        /**
         * Reads in the output of a stream on a background thread to keep the output pipe from backing up and freezing the called process.
         *
         * @param key The stdout or stderr key for this output capture.
         * @param controllerId Unique id of the controller.
         */
        public OutputCapture(
                final CAPTURE_POLICY capturedProcessStream,
                ProcessStream key,
                final String className,
                int controllerId) {
            this.key = key;
            this.capturedProcessStream = capturedProcessStream;
            this.controllerContextName = String.format(
                    "OutputCapture-%s-%d-%s", className, controllerId, key.name().toLowerCase()
            );
        }

        /**
         * Runs the capture.
         */
        @Override
        public CAPTURE_POLICY call() {
            try {
                // Delegate to the capture stream and let it's policy dictate how much to read
                capturedProcessStream.read();
            } catch (IOException e) {
                logger.error("Error reading process output", e);
            }
            return capturedProcessStream;
        }
    }
}

