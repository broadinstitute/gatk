package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.util.Log;
import com.google.common.collect.EnumHashBiMap;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.core.LoggerContext;
import org.apache.logging.log4j.core.config.Configuration;
import org.apache.logging.log4j.core.config.LoggerConfig;
import org.broadinstitute.hellbender.exceptions.GATKException;

/**
 * Logging utilities.
 *
 * Hellbender tools use the Picard Log.LogLevel enum as the type for VERBOSITY command line arguments (ideally
 * we would use a log4j level enum, but each log4j level is represented by a static object so there is no
 * proper built-in enum that is compatible with the command line argument framework). Therefore we use the
 * Picard enum as the Hellbender currency and convert back and forth between that and the log4j namespace
 * as necessary. Note that the two namespaces have very similar, but not identical level names, and log4j
 * has a wider set of level values than those supported by the Picard enum.
 */
 public class LoggingUtils {

    // Map between the logging level used throughout Hellbender code (which is the Picard Log.LogLevel enum),
    // and the log4j log Level values.
    private static EnumHashBiMap<Log.LogLevel, Level> loggingLevelNamespaceMap = null;
    private static String badLevelValue = "Unrecognized verbosity level. Must be one of (DEBUG/ERROR/INFO/WARNING)";

    private static EnumHashBiMap<Log.LogLevel, Level> getLoggingLevelNamespaceMap() {
        if (null == loggingLevelNamespaceMap) {
            loggingLevelNamespaceMap = EnumHashBiMap.create(Log.LogLevel.class);

            loggingLevelNamespaceMap.put(Log.LogLevel.DEBUG, Level.DEBUG);
            loggingLevelNamespaceMap.put(Log.LogLevel.ERROR, Level.ERROR);
            loggingLevelNamespaceMap.put(Log.LogLevel.WARNING, Level.WARN);
            loggingLevelNamespaceMap.put(Log.LogLevel.INFO, Level.INFO);
        }
        return loggingLevelNamespaceMap;
    }

    // Package-private for unit test access
    static Log.LogLevel levelFromLog4jLevel(Level log4jLevel) {
        try {
            return getLoggingLevelNamespaceMap().inverse().get(log4jLevel);
        }
        catch (IllegalArgumentException e) {
            throw new GATKException.ShouldNeverReachHereException(badLevelValue, e);
        }
    }

    private static Level levelToLog4jLevel(Log.LogLevel picardLevel) {
        try {
            return getLoggingLevelNamespaceMap().get(picardLevel);
        }
        catch (IllegalArgumentException e) {
            throw new GATKException.ShouldNeverReachHereException(badLevelValue, e);
        }
    }

    /**
     * Propagate the verbosity level to both Picard and log4j.
     */
    public static void setLoggingLevel(final Log.LogLevel VERBOSITY) {

        // Call the Picard API to establish the logging level used by Picard
        Log.setGlobalLogLevel(VERBOSITY);

        // Now establish the logging level used by log4j by propagating the requested
        // logging level to all loggers associated with our logging configuration.
        final LoggerContext loggerContext = (LoggerContext) LogManager.getContext(false);
        final Configuration loggerContextConfig = loggerContext.getConfiguration();
        final String contextClassName = LoggingUtils.class.getName();
        final LoggerConfig loggerConfig = loggerContextConfig.getLoggerConfig(contextClassName);

        loggerConfig.setLevel(levelToLog4jLevel(VERBOSITY));
        loggerContext.updateLoggers();
    }
}
