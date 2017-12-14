package org.broadinstitute.hellbender.utils.runtime;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;

/**
 * Base class for executors that find and run scripts in an external script engine process (R, Python, etc).
 *
 * Subclasses must implement:
 *
 *   {@link #getApproximateCommandLine}
 *   {@link #getScriptException}
 */
public abstract class ScriptExecutor {
    private static final Logger logger = LogManager.getLogger(org.broadinstitute.hellbender.utils.runtime.ScriptExecutor.class);

    protected final String externalScriptExecutableName;    // external program to run; e.g. "RScript" or "python"
    protected boolean ignoreExceptions = false;

    private final File externalScriptExecutablePath;        // File for path to externalScriptExecutable

    /**
     * @param externalScriptExecutableName Name of the script engine to run (i.e. "RScript" or "python")
     */
    public ScriptExecutor(final String externalScriptExecutableName) {
        Utils.nonNull(externalScriptExecutableName);
        this.externalScriptExecutableName = externalScriptExecutableName;
        this.externalScriptExecutablePath = RuntimeUtils.which(externalScriptExecutableName);
    }

    /**
     * @return true if the executable exists and can be found on the path
     */
    public boolean externalExecutableExists() {
        return externalScriptExecutablePath != null;
    }

    /**
     * Set to true to have the ScriptExecutor catch and ignore GATK exceptions.
     *
     * @param ignoreExceptions
     */
    public void setIgnoreExceptions(final boolean ignoreExceptions) {
        this.ignoreExceptions = ignoreExceptions;
    }

    /**
     * Return a (not necessarily executable) string representing the command line for this executor for error
     * reporting purposes.
     * @return Command line string.
     */
    public abstract String getApproximateCommandLine();

    protected void executableMissing() {
        throw new UserException.CannotExecuteScript(
                externalScriptExecutableName,
                String.format("Please add the %s directory to your environment ${PATH}", externalScriptExecutableName));
    }

    /**
     * Called by the script executor when error is encountered. Subclasses should override and return an exception
     * that is a subclass of ScriptExecutorException.
     *
     * @param message String with the cause of the exception.
     * @return a {#ScriptExecutorException}-derived exception object
     */
    public abstract ScriptExecutorException getScriptException(final String message);

    /**
     * Execute the script represented by the arguments in {@code commandLineArguments}.
     *
     * @param commandLineArguments
     * @return true if the command executed successfully, otherwise false
     */
    protected boolean executeCuratedArgs(final String[] commandLineArguments) {
        if (!externalExecutableExists()) {
            if (!ignoreExceptions) {
                executableMissing();
            } else {
                logger.warn("Skipping: " + getApproximateCommandLine());
                return false;
            }
        }

        try {
            final ProcessSettings processSettings = new ProcessSettings(commandLineArguments);
            //if debug is enabled, output the stdout and stderr, otherwise capture it to a buffer
            if (logger.isDebugEnabled()) {
                processSettings.getStdoutSettings().printStandard(true);
                processSettings.getStderrSettings().printStandard(true);
            } else {
                processSettings.getStdoutSettings().setBufferSize(8192);
                processSettings.getStderrSettings().setBufferSize(8192);
            }

            final ProcessController controller = ProcessController.getThreadLocal();

            if (logger.isDebugEnabled()) {
                logger.debug("Executing:");
                for (final String arg: commandLineArguments) {
                    logger.debug("  " + arg);
                }
            }
            final ProcessOutput po = controller.exec(processSettings);
            final int exitValue = po.getExitValue();
            logger.debug("Result: " + exitValue);

            if (exitValue != 0){
                final StringBuilder message = new StringBuilder();
                message.append(
                        String.format("\n%s exited with %d\nCommand Line: %s",
                                externalScriptExecutableName,
                                exitValue,
                                String.join(" ", commandLineArguments)));
                //if debug was enabled the stdout/error were already output somewhere
                if (!logger.isDebugEnabled()){
                    message.append(String.format("\nStdout: %s\nStderr: %s",
                            po.getStdout().getBufferString(),
                            po.getStderr().getBufferString()));
                }
                throw getScriptException(message.toString());
            }

            return true;
        } catch (final GATKException e) {
            if (!ignoreExceptions) {
                throw e;
            } else {
                logger.warn(e.getMessage());
                return false;
            }
        }
    }

}
