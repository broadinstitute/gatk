package org.broadinstitute.hellbender.utils.python;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.runtime.ScriptExecutor;
import org.broadinstitute.hellbender.utils.runtime.ScriptExecutorException;

/**
 * Base class for services for executing Python Scripts.
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
 */
public abstract class PythonExecutorBase extends ScriptExecutor {
    private static final Logger logger = LogManager.getLogger(PythonExecutorBase.class);

    /**
     * Enum of possible executables that can be launched by this executor.
     */
    public enum PythonExecutableName {

        PYTHON("python"),
        PYTHON3("python3");

        private final String executableName;

        PythonExecutableName (final String executableName) {
            this.executableName = executableName;
        }

        public String getExecutableName() {
            return executableName;
        }
    }

    // File extension used for python scripts
    public static final String PYTHON_EXTENSION = ".py";

    /**
     * @param ensureExecutableExists throw if the python executable cannot be located
     */
    public PythonExecutorBase(boolean ensureExecutableExists) {
        this(PythonExecutableName.PYTHON, ensureExecutableExists);
    }

    /**
     * @param pythonExecutableName name of the python executable to start
     * @param ensureExecutableExists throw if the python executable cannot be found
     */
    public PythonExecutorBase(final PythonExecutableName pythonExecutableName, final boolean ensureExecutableExists) {
        super(pythonExecutableName.getExecutableName());
        if (ensureExecutableExists && !externalExecutableExists()) {
            executableMissing();
        }
    }

    /**
     * Return an exception specific to this executor type, to be thrown on error conditions.
     * @param message
     */
    @Override
    public ScriptExecutorException getScriptException(final String message) {
        return new PythonScriptExecutorException(message);
    }

    /**
     * Return a (not necessarily executable) string representing the current command line for this executor
     * for error reporting purposes.
     * @return Command line string.
     */
    public abstract String getApproximateCommandLine();

}
