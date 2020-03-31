package org.broadinstitute.hellbender.engine.progressmeter;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledExecutorService;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;

/**
 * A basic progress meter to print out the number of records processed (and other metrics) during a traversal
 * at a configurable time interval.
 *
 * Clients set the update interval at construction, which controls how many seconds must elapse
 * before printing an update. Then call {@link #start} at traversal start, {@link #update(T)}
 * after processing each record from the primary input, and {@link #stop} at traversal end to print
 * summary statistics.
 *
 * Note that {@link #start} must only be called once, before any {@link #update(T)}.
 * Note no {@link #update(T)} must be called after {@link #stop}.
 *
 * All output is made at INFO level via log4j.
 */
public abstract class ProgressMeter<T> {
    private static final Logger logger = LogManager.getLogger(ProgressMeter.class);

    /**
     * By default, we output a line to the logger after this many seconds have elapsed
     */
    public static final double DEFAULT_SECONDS_BETWEEN_UPDATES = 10.0;

    /**
     * Number of milliseconds in a second
     */
    public static final long MILLISECONDS_PER_SECOND = 1000L;

    /**
     * Number of milliseconds in a minute
     */
    public static final long MILLISECONDS_PER_MINUTE = MILLISECONDS_PER_SECOND * 60L;

    /**
     * Default label for records in logger messages. For display purposes only.
     */
    public static final String DEFAULT_RECORD_LABEL = "records";

    /**
     * We output a line to the logger after approximately this many milliseconds have elapsed
     */
    private long millisecondsBetweenUpdates;

    /**
     * Total records processed
     */
    private long numRecordsProcessed = 0L;

    /**
     * Our start timestamp in milliseconds as returned by our {@link #getTime()}
     */
    private long startTimeMs = 0L;

    /**
     * Current timestamp in milliseconds as returned by our {@link #getTime()}.
     */
    private long currentTimeMs = 0L;

    /**
     * Timestamp in milliseconds as returned by {@link #getTime()} of the last time
     * we outputted a progress line to the logger
     */
    private long lastPrintTimeMs = 0L;

    /**
     * The most recently processed record, or null if the most recent update was a null.
     */
    private T currentRecord = null;

    /**
     * The number of times we've outputted a status line to the logger via {@link #printProgress}.
     * We keep track of this only for unit-testing purposes.
     */
    private long numLoggerUpdates = 0L;

    /**
     * Keeps track of whether the progress meter has ever been started.
     */
    private boolean started;

    /**
     * Keeps track of whether the progress meter has been stopped.
     */
    private boolean stopped;

    /**
     * Label for records in logger messages. For display purposes only.
     */
    private String recordLabel = DEFAULT_RECORD_LABEL;

    /**
     * Timer to run updates in the background
     */
    private final ScheduledExecutorService scheduler;

    /**
     * Create a progress meter with the default update interval of {@link #DEFAULT_SECONDS_BETWEEN_UPDATES} seconds.
     */
    public ProgressMeter() {
        this(DEFAULT_SECONDS_BETWEEN_UPDATES);
    }

    /**
     * Create a progress meter with a custom update interval.
     *
     * @param secondsBetweenUpdates number of seconds that should elapse before outputting a line to the logger
     */
    public ProgressMeter( final double secondsBetweenUpdates ) {
        Utils.validateArg(secondsBetweenUpdates > 0, "secondsBetweenUpdates must be > 0.0");
        this.started = false;
        this.stopped = false;
        this.millisecondsBetweenUpdates = (long)(secondsBetweenUpdates*(double)MILLISECONDS_PER_SECOND);
        Utils.validate(millisecondsBetweenUpdates > 0, "millisecondsBetweenUpdates must be > 0");

        this.scheduler = Executors.newScheduledThreadPool(1,
                new ThreadFactoryBuilder().setDaemon(true).setNameFormat("Progress Meter").build());
    }

    /**
     * Change the label used for records in logger messages. Default label is {@link #DEFAULT_RECORD_LABEL}
     *
     * @param label Label to use for records in logger messages. Not null.
     */
    public void setRecordLabel( final String label ) {
        Utils.nonNull(label);
        this.recordLabel = label;
    }

    /**
     * Start the progress meter and produce preliminary output such as column headings.
     * @throws IllegalStateException if the meter has been started before or has been stopped already
     */
    public synchronized void start() {
        Utils.validate( !started, "the progress meter has been started already");
        Utils.validate( !stopped, "the progress meter has been stopped already");
        started = true;
        logger.info("Starting traversal");
        printHeader();

        startTimeMs = getTime();
        currentTimeMs = startTimeMs;
        lastPrintTimeMs = startTimeMs;
        numRecordsProcessed = 0L;
        numLoggerUpdates = 0L;
        currentRecord = null;

        scheduler.scheduleAtFixedRate(this::printProgress, millisecondsBetweenUpdates, millisecondsBetweenUpdates, TimeUnit.MILLISECONDS);
    }

    private long getTime() {
        return System.currentTimeMillis();
    }

    /**
     * Signal to the progress meter that an additional record has been processed. Will output
     * statistics to the logger roughly every {@link #millisecondsBetweenUpdates} milliseconds.
     *
     * @param currentLocus the genomic location of the record just processed or null if the most recent record had no location.
     * @throws IllegalStateException if the meter has not been started yet or has been stopped already
     */
    public synchronized void update( final T currentLocus ) {
        Utils.validate(started, "the progress meter has not been started yet");
        Utils.validate( !stopped, "the progress meter has been stopped already");
        ++numRecordsProcessed;
        this.currentRecord = currentLocus;
    }

    /**
     * Stop the progress meter and output summary statistics to the logger
     * @throws IllegalStateException if the meter has not been started yet or has been stopped already
     */
    public synchronized void stop() {
        Utils.validate(started, "the progress meter has not been started yet");
        Utils.validate( !stopped, "the progress meter has been stopped already");
        this.stopped = true;
        currentTimeMs = getTime();
        // Output progress a final time at the end
        printProgress();
        scheduler.shutdown();
        logger.info(String.format("Traversal complete. Processed %d total %s in %.1f minutes.", numRecordsProcessed, recordLabel, elapsedTimeInMinutes()));
    }

    /**
     * Print column headings labelling the output from {@link #printProgress}
     */
    private void printHeader() {
        logger.info(String.format("%20s  %15s  %20s  %15s",
                                  "Current Locus", "Elapsed Minutes",
                                  StringUtils.capitalize(recordLabel) + " Processed",
                                  StringUtils.capitalize(recordLabel) + "/Minute"));
    }

    /**
     * Output traversal statistics to the logger.
     */
    private synchronized void printProgress() {
        currentTimeMs = getTime();
        lastPrintTimeMs = currentTimeMs;
        ++numLoggerUpdates;
        logger.info(String.format("%20s  %15.1f  %20d  %15.1f",
                                  formatRecord(currentRecord), elapsedTimeInMinutes(), numRecordsProcessed, processingRate()));
    }

    /**
     * Format the current record into a string for printing.  This should ideally give the user some information about
     * the state of progress.
     *
     * @param currentRecord the most recent update the ProgressMeter has recieved.
     * @return a String summarizing the state of progress from this given record.
     */
    protected abstract String formatRecord(final T currentRecord);

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
     */
    private double processingRate() {
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
     * @return number of records that have been processed by the progress meter (for unit testing purposes)
     */
    @VisibleForTesting
    long getNumRecordsProcessed(){
        return numRecordsProcessed;
    }

    /**
     * Returns whether the meter has been started. It returns false before the call to {@link #start} and true forever after.
     */
    public boolean started() {
        return started;
    }

    /**
     * Returns whether the meter has been stopped. It returns false before the call to {@link #stop} and true forever after.
     */
    public boolean stopped() {
        return stopped;
    }

}
