package org.broadinstitute.hellbender.utils.python;

import org.broadinstitute.hellbender.utils.runtime.ScriptExecutorException;

/**
 * Python script execution exception.
 */
public class PythonScriptExecutorException extends ScriptExecutorException {

    private static final long serialVersionUID = 0L;

    public PythonScriptExecutorException(String msg) {
        super(msg);
    }
}
