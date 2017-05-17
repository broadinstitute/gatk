package org.broadinstitute.hellbender.utils;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport;

import java.util.Collections;
import java.util.List;

public final class RetryUtils {
    private static final Logger logger = LogManager.getLogger(RetryUtils.class);

    private RetryUtils(){}

    public static void runWithRetries(final Runnable toRun, final int maxRetries, final int baseRetryDelay, final List<Class<? extends Exception>> retryableExceptions){
        Utils.containsNoNull(retryableExceptions, "Null exception type not allowed");
        Utils.nonEmpty(retryableExceptions);
        Utils.validate(maxRetries > 0, "MaxRetries can't be zero");
        int tries = 0;
        while( tries < maxRetries){
            try{
                toRun.run();
                break;
            } catch (final Exception e){
                if (retryableExceptions.stream()
                        .anyMatch( retryableException -> e.getClass().isAssignableFrom(retryableException))){
                    tries += 1;
                    final int delay = calculateDelay(baseRetryDelay, tries);
                    logger.warn("Encountered retryable failure. Try " + tries + "/" + maxRetries + ".");
                    logger.warn("Error was:\n" + ExceptionUtils.getStackTrace(e));
                    logger.warn("Waiting for " + delay + " before trying again.");
                    try {
                        Thread.sleep(delay);
                    } catch (final InterruptedException interruptException) {
                        throw new RuntimeException("Interrupted while waiting to retry.", interruptException);
                    }
                } else {
                    throw e;
                }
            }
        }
    }

    private static int calculateDelay(int baseRetryDelay, int iteration) {
        return baseRetryDelay * iteration ^ 2;
    }
}
