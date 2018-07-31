package org.broadinstitute.hellbender.utils.logging;

import com.google.common.annotations.VisibleForTesting;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * A logger wrapper class which only outputs the first warning provided to it
 */
public class OneShotLogger {
    @VisibleForTesting
    Logger logger;
    private boolean hasWarned = false;

    public OneShotLogger(final Class<?> clazz) {
        logger = LogManager.getLogger(clazz);
    }

    public OneShotLogger(final Logger logger) {
        this.logger = logger;
    }

    /*
     * Will write a warning only once for an instance, otherwise it doesn't emit a message
     */
    public void warn(String message) {
        if (!hasWarned) {
            logger.warn(message);
            hasWarned=true;
        }
    }
}
