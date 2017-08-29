package org.broadinstitute.hellbender.utils.python;

import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.runtime.ScriptExecutor;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Generic service for executing Python Scripts.
 *
 * <ul><li>
 * All tools that use PythonScriptExecutor must have a Java-based front-end, with standard GATK (Barclay-based) arguments.
 * <li>
 * Minimize the amount of code written in Python -- as much of each tool's work as possible should be done in Java. In
 * particular, reading/writing final inputs and outputs should happen in Java. This is important to ensure the goal of
 * universal GCS support, consistent Google authentication handling, etc.
 *<li>
 * The names of any files that are created by Python code should be passed in to the python code from Java.
 *<li>
 *<li>
 * All dependencies (Python and native) of Python libraries used should be clearly documented, and included in the default
 * GATK docker image.
 * </li></ul>
 *
 * NOTE: Serial executions on the same PythonScriptExecutor are each run in a new process. No python state is retained
 * between command/script/module execution. Using -i doesn't buy you anything (for this version of the executor, at
 * least) since the process is terminated after each command completes.
 */
public class PythonScriptExecutor extends ScriptExecutor<PythonScriptExecutorException> {
    private static final Logger logger = LogManager.getLogger(PythonScriptExecutor.class);

    private static String python = "python";    // default executor name
    private static String python3 = "python3";
    public static String pyExtension = ".py";

    private final List<String> curatedCommandLineArgs = new ArrayList<>();

    /**
     * @param ensureExists throw if the python executor cannot be located
     */
    public PythonScriptExecutor(boolean ensureExists) {
        this(python, ensureExists);
    }

    /**
     * @param pythonExecutorName name of the python executable to start
     * @param ensureScriptEngineExists throw if the python executor cannot be found
     */
    public PythonScriptExecutor(final String pythonExecutorName, final boolean ensureScriptEngineExists) {
        super(pythonExecutorName);
        if (!pythonExecutorName.equals(python) && !pythonExecutorName.equals(python3)) {
            throw new IllegalArgumentException("python executable name must be either 'python' or 'python3'");
        }
        if (ensureScriptEngineExists && !getExternalExecutorExists()) {
            executorMissing();
        }
    }

    /**
     * Execute a python command (-c). No intermediate shell is created.
     *
     * @param command python command to be executed
     * @param pythonProcessArgs args to be passed to the python process
     * @param scriptArgs args to be passed to the python code
     * @return true if the command succeeds, otherwise false
     */
    public boolean executeCommand(final String command, final List<String> pythonProcessArgs, final List<String> scriptArgs) {
        Utils.nonNull(command, "Command string cannot be null");

        final List<String> args = new ArrayList<>();
        if (pythonProcessArgs != null) {
            args.addAll(pythonProcessArgs);
        }
        args.add("-c");
        args.add(command);
        if (scriptArgs != null) {
            args.addAll(scriptArgs);
        }
        return executeArgs(args);
    }

    /**
     * Execute a python module (-m). Modules must be on sys.path
     *
     * @param moduleName name of the module to execute
     * @param pythonProcessArgs args to be passed to the python process
     * @param scriptArgs args to be passed to the python code
     * @return true if the command succeeds, otherwise false
     */
    public boolean executeModule(final String moduleName, final List<String> pythonProcessArgs, final List<String> scriptArgs) {
        Utils.nonNull(moduleName, "module name cannot be null");
        if (moduleName.endsWith(pyExtension)) {
            throw new IllegalArgumentException(String.format("\"%s\" suffix should not be included to run a Python module", pyExtension));
        }

        final List<String> args = new ArrayList<>();
        if (pythonProcessArgs != null) {
            args.addAll(pythonProcessArgs);
        }
        args.add("-m");
        args.add(moduleName);
        if (scriptArgs != null) {
            args.addAll(scriptArgs);
        }
        return executeArgs(args);
    }

    /**
     * Execute a python script from a Resource file.
     *
     * @param scriptResource {@link Resource} for the script to execute
     * @param pythonProcessArgs args to be passed to the python process
     * @param scriptArgs args to be passed to the python code
     * @return true if the command succeeds, otherwise false
     */
    public boolean executeScript(final Resource scriptResource, final List<String> pythonProcessArgs, final List<String> scriptArgs) {
        Utils.nonNull(scriptResource, "script resource cannot be null");
        final File tempResourceFile = IOUtils.writeTempResource(scriptResource);

        try {
            return executeScript(tempResourceFile.getAbsolutePath(), pythonProcessArgs, scriptArgs);
        } finally {
            FileUtils.deleteQuietly(tempResourceFile);
        }
    }

    /**
     * Execute a python script.
     *
     * @param scriptName full path name of the script to execute
     * @param pythonProcessArgs args to be passed to the python process
     * @param scriptArgs args to be passed to the python code
     * @return true if the command succeeds
     */
    public boolean executeScript(final String scriptName, final List<String> pythonProcessArgs, final List<String> scriptArgs) {
        Utils.nonNull(scriptName, "script name cannot be null");
        if (!scriptName.endsWith(pyExtension)) {
            throw new IllegalArgumentException(String.format("Python script name (%s) must end with \"%s\"",
                    scriptName,
                    pyExtension));
        }

        final List<String> args = new ArrayList<>();
        if (pythonProcessArgs != null) {
            args.addAll(pythonProcessArgs);
        }
        args.add(scriptName);
        if (scriptArgs != null) {
            args.addAll(scriptArgs);
        }
        return executeArgs(args);
    }

    /**
     * Executes the Python executor using the values in {@code rawArgs}
     *
     * @param rawArgs raw command line arguments to be passed to the Python process
     * @return true if the command succeeds, otherwise false
     */
    public boolean executeArgs(final List<String> rawArgs) {
        Utils.nonNull(rawArgs, "Raw args cannot be null");

        // executor name first, followed by rawArgs
        curatedCommandLineArgs.clear();
        curatedCommandLineArgs.add(externalScriptExecutorName);
        rawArgs.forEach(curatedCommandLineArgs::add);

        try {
            // actually run the script
            return executeCuratedArgs(curatedCommandLineArgs.toArray(new String[curatedCommandLineArgs.size()]));
        } catch (final GATKException e) {
            if (!ignoreExceptions) {
                throw e;
            } else {
                logger.warn(e.getMessage());
                return false;
            }
        }
    }

    /**
     * Return an exception specific to this executor type, to be thrown on error conditions.
     * @param message
     */
    @Override
    public PythonScriptExecutorException getScriptException(final String message) {
        return new PythonScriptExecutorException(message.toString());
    }

    /**
     * Return a (not necessarily executable) string representing the current command line for this executor
     * for error reporting purposes.
     * @return Command line string.
     */
    public String getApproximateCommandLine() {
        return curatedCommandLineArgs.stream().collect(Collectors.joining(" "));
    }

}
