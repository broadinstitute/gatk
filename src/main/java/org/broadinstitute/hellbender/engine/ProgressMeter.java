package org.broadinstitute.hellbender.engine;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.function.LongSupplier;

/**
 * A basic progress meter to print out the number of records processed (and other metrics) during a traversal
 * at a configurable time interval.
 *
 * Clients set the update interval at construction, which controls how many seconds must elapse
 * before printing an update. Then call {@link #start} at traversal start, {@link #update(Locatable)}
 * after processing each record from the primary input, and {@link #stop} at traversal end to print
 * summary statistics.
 *
 * All output is made at INFO level via log4j.
 */
public final class ProgressMeter {
    protected static final Logger logger = LogManager.getLogger(ProgressMeter.class);

    /**
     * By default, we output a line to the logger after this many seconds have elapsed
     */
    public static final double DEFAULT_SECONDS_BETWEEN_UPDATES = 10.0;

    /**
     * We check the current time every time we process this many records
     * (to cut down on the number of system calls)
     */
    public static final long RECORDS_BETWEEN_TIME_CHECKS = 1000l;

    /**
     * By default, we use this function to get the current time
     */
    public static final LongSupplier DEFAULT_TIME_FUNCTION = System::currentTimeMillis;

    /**
     * Number of milliseconds in a second
     */
    public static final long MILLISECONDS_PER_SECOND = 1000l;

    /**
     * Number of milliseconds in a minute
     */
    public static final long MILLISECONDS_PER_MINUTE = MILLISECONDS_PER_SECOND * 60l;

    /**
     * We output a line to the logger after this many seconds have elapsed
     */
    private final double secondsBetweenUpdates;

    /**
     * Total records processed
     */
    private long numRecordsProcessed = 0l;

    /**
     * Our start timestamp in milliseconds as returned by our {@link #timeFunction}
     */
    private long startTimeMs = 0l;

    /**
     * Current timestamp in milliseconds as returned by our {@link #timeFunction}.
     *
     * Updated only every {@link #RECORDS_BETWEEN_TIME_CHECKS} records to cut down
     * on system calls.
     */
    private long currentTimeMs = 0l;

    /**
     * Timestamp in milliseconds as returned by our {@link #timeFunction} of the last time
     * we outputted a progress line to the logger
     */
    private long lastPrintTimeMs = 0l;

    /**
     * The genomic location of the most recently processed record, or null if the most recent record had no location.
     * Updated only when we actually output a line to the logger.
     */
    private Locatable currentLocus = null;

    /**
     * The number of times we've outputted a status line to the logger via {@link #printProgress}.
     * We keep track of this only for unit-testing purposes.
     */
    private long numLoggerUpdates = 0l;

    /**
     * Function that returns the current time in milliseconds (defaults to {@link #DEFAULT_TIME_FUNCTION}).
     * Should only be customized for unit testing purposes.
     */
    private LongSupplier timeFunction;

    /**
     * Create a progress meter with the default update interval of {@link #DEFAULT_SECONDS_BETWEEN_UPDATES} seconds
     * and the default time function {@link #DEFAULT_TIME_FUNCTION}.
     */
    public ProgressMeter() {
        this(DEFAULT_SECONDS_BETWEEN_UPDATES);
    }

    /**
     * Create a progress meter with a custom update interval and the default time function {@link #DEFAULT_TIME_FUNCTION}
     *
     * @param secondsBetweenUpdates number of seconds that should elapse before outputting a line to the logger
     */
    public ProgressMeter( final double secondsBetweenUpdates ) {
        this(secondsBetweenUpdates, DEFAULT_TIME_FUNCTION);
    }

    /**
     * Create a progress meter with a custom update interval and a custom function for getting the current
     * time in milliseconds.
     *
     * Providing your own time function is only useful in unit tests -- in normal usage
     * clients should call one of the other constructors.
     *
     * @param secondsBetweenUpdates number of seconds that should elapse before outputting a line to the logger
     * @param timeFunction function that returns the current time in milliseconds.
     */
    @VisibleForTesting
    ProgressMeter( final double secondsBetweenUpdates, final LongSupplier timeFunction ) {
        Utils.nonNull(timeFunction);
        if ( secondsBetweenUpdates <= 0.0 ) {
            throw new IllegalArgumentException("secondsBetweenUpdates must be > 0.0");
        }

        this.secondsBetweenUpdates = secondsBetweenUpdates;
        this.timeFunction = timeFunction;
    }

    /**
     * Start the progress meter and produce preliminary output such as column headings.
     */
    public void start() {
        logger.info("Starting traversal");
        printHeader();

        startTimeMs = timeFunction.getAsLong();
        currentTimeMs = startTimeMs;
        lastPrintTimeMs = startTimeMs;
        numRecordsProcessed = 0l;
        numLoggerUpdates = 0l;
        currentLocus = null;
    }

    /**
     * Signal to the progress meter that an additional record has been processed. Will output
     * statistics to the logger roughly every {@link #secondsBetweenUpdates} seconds.
     *
     * @param currentLocus the genomic location of the record just processed or null if the most recent record had no location.
     */
    public void update( final Locatable currentLocus ) {
        ++numRecordsProcessed;
        if ( numRecordsProcessed % RECORDS_BETWEEN_TIME_CHECKS == 0 ) {
            currentTimeMs = timeFunction.getAsLong();
            this.currentLocus = currentLocus;

            if ( secondsSinceLastPrint() >= secondsBetweenUpdates ) {
                printProgress();
                lastPrintTimeMs = currentTimeMs;
            }
        }
    }

    /**
     * Stop the progress meter and output summary statistics to the logger
     */
    public void stop() {
        currentTimeMs = timeFunction.getAsLong();
        logger.info(String.format("Traversal complete. Processed %d total records in %.1f minutes.", numRecordsProcessed, elapsedTimeInMinutes()));
    }

    /**
     * Print column headings labelling the output from {@link #printProgress}
     */
    private void printHeader() {
        logger.info(String.format("%20s  %15s  %20s  %15s",
                                  "Current Locus", "Elapsed Minutes", "Records Processed", "Records/Minute"));
    }

    /**
     * Output traversal statistics to the logger.
     */
    private void printProgress() {
        ++numLoggerUpdates;
        logger.info(String.format("%20s  %15.1f  %20d  %15.1f",
                                  currentLocusString(), elapsedTimeInMinutes(), numRecordsProcessed, processingRate()));
    }

    /**
     * @return the total minutes elapsed since we called {@link #start}
     *
     * This is only accurate at set polling intervals and should not be
     * called directly except in tests.
     */
    @VisibleForTesting
    double elapsedTimeInMinutes() {
        return (currentTimeMs - startTimeMs) / (double)MILLISECONDS_PER_MINUTE;
    }

    /**
     * @return the number of seconds that have elapsed since our last progress output to the logger
     *
     * This is only accurate at set polling intervals and should not be
     * called directly except in tests.
     */
    @VisibleForTesting
    double secondsSinceLastPrint() {
        return (currentTimeMs - lastPrintTimeMs) / (double)MILLISECONDS_PER_SECOND;
    }

    /**
     * @return number of records we're processing per minute, on average
     *
     * This is only accurate at set polling intervals and should not be
     * called directly except in tests.
     */
    @VisibleForTesting
    double processingRate() {
        return numRecordsProcessed / elapsedTimeInMinutes();
    }

    /**
     * @return number of times we've outputted a progress line to the logger (for unit testing purposes)
     */
    @VisibleForTesting
    long numLoggerUpdates() {
        return numLoggerUpdates;
    }

    /**
     * @return genomic location of the most recent record formatted for output to the logger
     */
    private String currentLocusString() {
        return currentLocus != null ? currentLocus.getContig() + ":" + currentLocus.getStart() :
                                      "unmapped";
    }
}
