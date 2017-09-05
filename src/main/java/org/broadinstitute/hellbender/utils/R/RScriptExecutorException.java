package org.broadinstitute.hellbender.utils.R;

import org.broadinstitute.hellbender.utils.runtime.ScriptExecutorException;

public final class RScriptExecutorException extends ScriptExecutorException {
    private static final long serialVersionUID = 0L;

    public RScriptExecutorException(String msg) {
        super(msg);
    }
}
