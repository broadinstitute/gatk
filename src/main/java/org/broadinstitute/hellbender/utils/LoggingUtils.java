package org.broadinstitute.hellbender.utils;

import com.google.common.collect.BiMap;
import com.google.common.collect.EnumHashBiMap;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.core.Layout;
import org.apache.logging.log4j.core.LoggerContext;
import org.apache.logging.log4j.core.appender.FileAppender;
import org.apache.logging.log4j.core.config.Configuration;
import org.apache.logging.log4j.core.config.LoggerConfig;
import org.apache.logging.log4j.core.layout.PatternLayout;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.*;
import java.util.logging.*;

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

    private static final String PATTERN_STRING = "%-5p %d{HH:mm:ss,SSS} %C{1} - %m %n";

    /**
     * Class that adds the ability for MinLog to write to file
     */
    static class FileMinLogger extends com.esotericsoftware.minlog.Log.Logger {
        static private BufferedWriter writer;

        public FileMinLogger(final String path) {
            try {
                writer = new BufferedWriter(new FileWriter(path));
            } catch (final IOException ex) {
                System.err.println("Could not construct FileMinLogger : " + ex.getMessage());
                ex.printStackTrace();
            }
        }

        @Override
        protected void print(final String message) {
            try {
                writer.write(message);
            } catch (final IOException ex) {
                System.err.println("Could not write " + message + " : " + ex.getMessage());
            }
        }

        protected void close() {
            try {
                writer.close();
            } catch (final IOException ex) {
                System.err.println("Could not close writer : " + ex.getMessage());
            }
        }
    }

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

    // Map between the logging level used throughout Hellbender code (which is the Picard Log.LogLevel enum),
    // and the Java core log Level values.
    private static BiMap<Log.LogLevel, java.util.logging.Level> javaUtilLevelNamespaceMap;
    static {
        javaUtilLevelNamespaceMap = EnumHashBiMap.create(Log.LogLevel.class);
        javaUtilLevelNamespaceMap.put(Log.LogLevel.ERROR, java.util.logging.Level.SEVERE);
        javaUtilLevelNamespaceMap.put(Log.LogLevel.WARNING, java.util.logging.Level.WARNING);
        javaUtilLevelNamespaceMap.put(Log.LogLevel.INFO, java.util.logging.Level.INFO);
        javaUtilLevelNamespaceMap.put(Log.LogLevel.DEBUG, java.util.logging.Level.FINEST);
    }

    /**
     * Set {@link htsjdk.samtools.util.Log} level
     */
    private static void setHtsjdkLoggingLevel(final Log.LogLevel verbosity){
        Log.setGlobalLogLevel(verbosity);
    }

    // Package-private for unit test access
    static Log.LogLevel levelFromLog4jLevel(Level log4jLevel) {
        return loggingLevelNamespaceMap.inverse().get(log4jLevel);
    }

    private static Level levelToLog4jLevel(Log.LogLevel picardLevel) {
        return loggingLevelNamespaceMap.get(picardLevel);
    }

    /**
     * Set the {@link java.util.logging} level since some of our dependencies write messages using it.
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

        for (Handler handler : topLogger.getHandlers()) {
            handler.setLevel(levelToJavaUtilLevel(verbosity));
        }
    }

    private static java.util.logging.Level levelToJavaUtilLevel(Log.LogLevel picardLevel) {
        return javaUtilLevelNamespaceMap.get(picardLevel);
    }

    /**
     * set the logging level for {@link org.apache.logging.log4j.core.Logger}
     */
    private static void setLog4JLoggingLevel(final Log.LogLevel verbosity) {
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

    /**
     * Propagate the verbosity level to HTSJDK, log4j, the java built in logger, and Kryo's MinLog
     */
    public static void setLoggingLevel(final Log.LogLevel verbosity) {

        // set the HTSJDK logging level, used by Picard
        setHtsjdkLoggingLevel(verbosity);

        // set the Log4JLoggingLevel
        setLog4JLoggingLevel(verbosity);

        // set the java.util.logging Level
        setJavaUtilLoggingLevel(verbosity);

        // set the Esotericsoft MinLog level, this is used by kryo
        setMinLogLoggingLevel(verbosity);
    }

    /**
     * Set the {@link htsjdk.samtools.util.Log} to log to file.
     */
    private static void setHtsjdkLoggingFile(final String path) {
        try {
            Log.setGlobalPrintStream(new PrintStream(new FileOutputStream(path, true)));
        } catch (final FileNotFoundException ex) {
            System.err.println("Could not send log to " + path +  " for htsjdk.samtools.util.Log : " + ex.getMessage());
        }
    }

    /**
     * Set the {@link org.apache.logging.log4j.core.Logger} to log to file.
     */
    private static void setLog4JLoggingFile(final String path) {
        final LoggerContext loggerContext = (LoggerContext) LogManager.getContext(false);
        final Configuration loggerContextConfig = loggerContext.getConfiguration();
        final Layout<?> layout = PatternLayout.createLayout(PATTERN_STRING, null, loggerContextConfig,
                null, null, false, false, null, null);
        final FileAppender appender = FileAppender.createAppender(path, "true", "", path, "",
                "", "", "", layout, null, "", path, null);
        appender.start();
        loggerContextConfig.addLoggerAppender((org.apache.logging.log4j.core.Logger) LogManager.getRootLogger(), appender);
    }

    /**
     * Set the {@link java.util.logging} to log to file.
     * Taken from
     * https://stackoverflow.com/questions/15758685/how-to-write-logs-in-text-file-when-using-java-util-logging-logger
     */
    private static void setJavaUtilLoggingFile(final String path) {
        final Logger topLogger = java.util.logging.Logger.getLogger("");

        try {
            final FileHandler fh = new FileHandler(path);
            topLogger.addHandler(fh);
            final SimpleFormatter formatter = new SimpleFormatter();
            fh.setFormatter(formatter);
        } catch (final IOException ex) {
            System.err.println("Could not send log to " + path +  " for java.util.logging.Logger : " + ex.getMessage());
        }
    }

    /**
     * Set the {@link com.esotericsoftware.minlog.Log} to log to file
     *
     */
    private static void setMinLogLoggingFile(final String path) {
        com.esotericsoftware.minlog.Log.setLogger(new FileMinLogger(path));
    }

    /**
     * Send logging output to file if file name is not null.
     *
     * @param path  File path for log output
     */
    public static void setLoggingFile(final String path) {
        if (path != null) {
            // set HTSJDK log to output to file
            setHtsjdkLoggingFile(path);

            // set log4j log to output to file
            setLog4JLoggingFile(path);

            // set the java.util.logging output to file
            setJavaUtilLoggingFile(path);

            // set the Esotericsoft MinLog output to file, this is used by kryo
            setMinLogLoggingFile(path);
        }
    }
}
