package org.broadinstitute.hellbender.utils;

import com.google.common.collect.BiMap;
import com.google.common.collect.EnumHashBiMap;
import htsjdk.samtools.util.Log;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.core.LoggerContext;
import org.apache.logging.log4j.core.config.Configuration;
import org.apache.logging.log4j.core.config.LoggerConfig;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.logging.ConsoleHandler;
import java.util.logging.Handler;
import java.util.logging.Logger;

/**
 * Logging utilities.
 *
 * Hellbender tools use the Picard Log.LogLevel enum as the type for VERBOSITY command line arguments (ideally
 * we would use a log4j level enum, but each log4j level is represented by a static object so there is no
 * proper built-in enum that is compatible with the command line argument framework). Therefore we use the
 * Picard enum as the Hellbender currency and convert back and forth between that and the log4j namespace
 * as necessary. Note that the two namespaces have very similar, but not identical level names, and log4j
 * has a wider set of level values than those supported by the Picard enum.  Likewise the java built in logger has
 * a different set of level values which we chose a mapping to.
 */
public class LoggingUtils {

    // Map between the logging level used throughout Hellbender code (which is the Picard Log.LogLevel enum),
    // and the log4j log Level values.
    private static BiMap<Log.LogLevel, Level> loggingLevelNamespaceMap;
    static {
        loggingLevelNamespaceMap = EnumHashBiMap.create(Log.LogLevel.class);
        loggingLevelNamespaceMap.put(Log.LogLevel.ERROR, Level.ERROR);
        loggingLevelNamespaceMap.put(Log.LogLevel.WARNING, Level.WARN);
        loggingLevelNamespaceMap.put(Log.LogLevel.INFO, Level.INFO);
        loggingLevelNamespaceMap.put(Log.LogLevel.DEBUG, Level.DEBUG);
    }

    private static BiMap<Log.LogLevel, java.util.logging.Level> javaUtilLevelNamespaceMap;
    static {
        javaUtilLevelNamespaceMap = EnumHashBiMap.create(Log.LogLevel.class);
        javaUtilLevelNamespaceMap.put(Log.LogLevel.ERROR, java.util.logging.Level.SEVERE);
        javaUtilLevelNamespaceMap.put(Log.LogLevel.WARNING, java.util.logging.Level.WARNING);
        javaUtilLevelNamespaceMap.put(Log.LogLevel.INFO, java.util.logging.Level.INFO);
        javaUtilLevelNamespaceMap.put(Log.LogLevel.DEBUG, java.util.logging.Level.FINEST);
    }


    // Package-private for unit test access
    static Log.LogLevel levelFromLog4jLevel(Level log4jLevel) {
        return loggingLevelNamespaceMap.inverse().get(log4jLevel);
    }

    /**
     * Converts a picard log level to a log4j log level.
     * @param picardLevel Picard {@link Log.LogLevel} to convert to a Log4J {@link Level}.
     * @return The {@link Level} that corresponds to the given {@code picardLevel}.
     */
    public static Level levelToLog4jLevel(Log.LogLevel picardLevel) {
        return loggingLevelNamespaceMap.get(picardLevel);
    }


    /**
     * Set the java.util.logging level since some of our dependencies write messages using it.
     *
     * Taken with modifications from a post by michel.iamit
     * http://stackoverflow.com/users/369060/michel-iamit
     * from http://stackoverflow.com/questions/470430/java-util-logging-logger-doesnt-respect-java-util-logging-level
     */
    private static void setJavaUtilLoggingLevel(final Log.LogLevel verbosity) {
        Logger topLogger = java.util.logging.Logger.getLogger("");

        Handler consoleHandler = null;
        for (Handler handler : topLogger.getHandlers()) {
            if (handler instanceof ConsoleHandler) {
                consoleHandler = handler;
                break;
            }
        }

        if (consoleHandler == null) {
            consoleHandler = new ConsoleHandler();
            topLogger.addHandler(consoleHandler);
        }
        consoleHandler.setLevel(levelToJavaUtilLevel(verbosity));
    }

    private static java.util.logging.Level levelToJavaUtilLevel(Log.LogLevel picardLevel) {
        return javaUtilLevelNamespaceMap.get(picardLevel);
    }


    /**
     * Propagate the verbosity level to Picard, log4j, the java built in logger, and Kryo's MinLog
     */
    public static void setLoggingLevel(final Log.LogLevel verbosity) {

        // Call the Picard API to establish the logging level used by Picard
        Log.setGlobalLogLevel(verbosity);

        // set the Log4JLoggingLevel
        setLog4JLoggingLevel(verbosity);

        // set the java.util.logging Level
        setJavaUtilLoggingLevel(verbosity);

        // set the esotericsoft MinLog level, this is used by kryo
        setMinLogLoggingLevel(verbosity);
    }

    private static void setLog4JLoggingLevel(Log.LogLevel verbosity) {
        // Now establish the logging level used by log4j by propagating the requested
        // logging level to all loggers associated with our logging configuration.
        final LoggerContext loggerContext = (LoggerContext) LogManager.getContext(false);
        final Configuration loggerContextConfig = loggerContext.getConfiguration();
        final String contextClassName = LoggingUtils.class.getName();
        final LoggerConfig loggerConfig = loggerContextConfig.getLoggerConfig(contextClassName);

        loggerConfig.setLevel(levelToLog4jLevel(verbosity));
        loggerContext.updateLoggers();
    }

    /**
     * set the logging level for {@link com.esotericsoftware.minlog.Log}, the logger used by Kryo
     */
    private static void setMinLogLoggingLevel(Log.LogLevel verbosity) {
        switch (verbosity) {
            case DEBUG: com.esotericsoftware.minlog.Log.DEBUG(); break;
            case INFO: com.esotericsoftware.minlog.Log.INFO(); break;
            case WARNING: com.esotericsoftware.minlog.Log.WARN(); break;
            case ERROR: com.esotericsoftware.minlog.Log.ERROR(); break;
            default:
                throw new GATKException("This log level is not implemented properly: " + verbosity);
        }
    }
}
