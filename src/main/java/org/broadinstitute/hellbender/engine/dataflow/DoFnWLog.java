package org.broadinstitute.hellbender.engine.dataflow;

import com.google.cloud.dataflow.sdk.transforms.DoFn;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import org.broadinstitute.hellbender.utils.logging.BunnyLog;

/**
 * A DoFn that automatically logs its start and end.
 * When subclassing, if you want to write your own start/finishBundle then
 * make sure to call super.startBundle() and super.finishBundle().
 */
public abstract class DoFnWLog<T,U> extends DoFn<T,U> {
    private static final long serialVersionUID = 1L;

    protected String className;
    protected String opName;
    protected transient BunnyLog bunny;
    protected transient Logger log = null;

    public DoFnWLog(String name) {
        this.opName = this.className = name;
    }

    public DoFnWLog(String className, String operationName) {
        this.className = className;
        this.opName = operationName;
    }

    @Override
    public void startBundle(DoFn<T, U>.Context c) throws Exception {
        log = LogManager.getLogger(className);
        bunny = new BunnyLog(log);
        bunny.start(opName);
    }
    
    @Override
    public void finishBundle(DoFn<T, U>.Context c) throws Exception {
        bunny.end();
    }
}
