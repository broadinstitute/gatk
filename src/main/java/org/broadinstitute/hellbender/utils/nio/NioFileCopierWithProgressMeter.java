package org.broadinstitute.hellbender.utils.nio;

import org.apache.commons.collections4.queue.CircularFifoQueue;
import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

import javax.xml.bind.DatatypeConverter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.time.Duration;
import java.time.Instant;
import java.time.ZoneId;
import java.util.Queue;

/**
 * Class to copy a file using {@link java.nio}.
 * Operates using paths.
 *
 * INSTANCES OF THIS CLASS ARE NOT THREAD-SAFE!
 *
 * Created by jonn on 8/27/18.
 */
public class NioFileCopierWithProgressMeter {

    //==================================================================================================================
    // Standard logger:
    private static final Logger logger = LogManager.getLogger(NioFileCopierWithProgressMeter.class);

    //==================================================================================================================
    // Public Static Members:

    //==================================================================================================================
    // Private Static Members:

    protected static final int    BUFFER_SIZE_BYTES                          = 1024 * 1024;
    protected static final double DEFAULT_PROGRESS_DISPLAY_PERCENT_INCREMENT = 0.25;

    protected static final long SECOND_IN_MS = 1000;
    protected static final long MINUTE_IN_MS = SECOND_IN_MS * 60;
    protected static final long HOUR_IN_MS   = MINUTE_IN_MS * 60;
    protected static final long DAY_IN_MS    = HOUR_IN_MS * 24;

    protected static final int KB_TO_BYTES      = 1024;
    protected static final int MS_TO_SEC        = 1000;
    protected static final int NANOS_TO_MILLIS  = 1000000;
    protected static final int NANOS_TO_SECONDS = 1000000000;

    protected static final int COPY_SPEED_HISTORY_SIZE = 10;

    protected static final boolean OVERWRITE_EXISTING_DEFAULT   = false;
    protected static final Verbosity VERBOSITY_DEFAULT          = Verbosity.MODERATE;

    //==================================================================================================================
    // Private Members:

    // Data variables:
    protected final Path source;
    protected final Path dest;

    protected long srcFileSize;
    protected int srcFileSizeNumDigits;

    protected String        checksum         = "";
    protected MessageDigest messageDigest    = null;
    protected String        expectedChecksum = "";

    // Flag defaults:
    protected boolean   overwriteExisting              = OVERWRITE_EXISTING_DEFAULT;
    protected Verbosity verbosity                      = Verbosity.MODERATE;
    protected boolean   formatTimeRemainingAsTimestamp = true;

    // Copy buffer:
    protected final byte          copyBuffer[]                = new byte[ BUFFER_SIZE_BYTES ];

    // Progress variables:
    protected double progressPercentDisplayIncrement             = DEFAULT_PROGRESS_DISPLAY_PERCENT_INCREMENT;
    protected final Queue<Double> downloadBytesPerMilliSecond = new CircularFifoQueue<>(COPY_SPEED_HISTORY_SIZE);
    protected       boolean       copyComplete                = false;
    protected       long          totalBytesRead              = 0;
    protected       long          progressBytesRead           = 0;
    protected       double        lastProgressValue           = 0;
    protected long lastProgressTime_ns;

    //==================================================================================================================
    // Constructors:

    /**
     * {@link NioFileCopierWithProgressMeter} uses a factory pattern.
     * This internal constructor is to be used by the class itself.
     * @param source The {@link Path} to the source file for the copy.
     * @param dest The {@link Path} to the destination file for the copy.
     * @param overwriteExisting If {@code true} will overwrite an existing file in the location specified by {@code dest}.
     * @param verbosity {@link Verbosity} of the progress progress log over the duration of the copy.
     */
    protected NioFileCopierWithProgressMeter(final Path source, final Path dest, final boolean overwriteExisting, final Verbosity verbosity) {
        this.source = source.toAbsolutePath();
        this.dest = dest.toAbsolutePath();
        this.overwriteExisting = overwriteExisting;
        this.verbosity = verbosity;
    }

    //==================================================================================================================
    // Static Methods:

    /**
     * Create an {@link NioFileCopierWithProgressMeter}.
     * By default the resulting {@link NioFileCopierWithProgressMeter} will not overwrite the destination if anything already exists there.
     * @param source The {@link Path} to the source file for the copy.
     * @param dest The {@link Path} to the destination file for the copy.
     * @return An {@link NioFileCopierWithProgressMeter} initialized to copy the file located at {@code source} to the location specified by {@code dest}.
     */
    public static NioFileCopierWithProgressMeter create(final Path source, final Path dest) {
        return create(source, dest, OVERWRITE_EXISTING_DEFAULT);
    }

    /**
     * Create an {@link NioFileCopierWithProgressMeter}.
     * Will periodically display progress of copying files.
     * @param source The {@link Path} to the source file for the copy.
     * @param dest The {@link Path} to the destination file for the copy.
     * @param overwriteExisting If {@code true} will overwrite an existing file in the location specified by {@code dest}.
     * @return An {@link NioFileCopierWithProgressMeter} initialized to copy the file located at {@code source} to the location specified by {@code dest}.
     */
    public static NioFileCopierWithProgressMeter create(final Path source, final Path dest, final boolean overwriteExisting) {
        return new NioFileCopierWithProgressMeter(source, dest, overwriteExisting, VERBOSITY_DEFAULT);
    }

    /**
     * Create an {@link NioFileCopierWithProgressMeter}.
     * @param source The {@link Path} to the source file for the copy.
     * @param dest The {@link Path} to the destination file for the copy.
     * @param overwriteExisting If {@code true} will overwrite an existing file in the location specified by {@code dest}.
     * @param verbosity {@link Verbosity} of the progress progress log over the duration of the copy.
     * @return An {@link NioFileCopierWithProgressMeter} initialized to copy the file located at {@code source} to the location specified by {@code dest}.
     */
    public static NioFileCopierWithProgressMeter create(final Path source, final Path dest, final boolean overwriteExisting, final Verbosity verbosity) {
        return new NioFileCopierWithProgressMeter(source, dest, overwriteExisting, verbosity);
    }

    //==================================================================================================================
    // Getters / Setters:

    /**
     * @return A copy of the {@link Path} used as the source for this {@link NioFileCopierWithProgressMeter}.
     */
    public Path getSource() {
        return source;
    }

    /**
     * @return A copy of the {@link Path} used as the destination for this {@link NioFileCopierWithProgressMeter}.
     */
    public Path getDest() {
        return dest;
    }

    /**
     * @return {@code true} iff the {@link #source} has been copied to the {@link #dest} in this {@link NioFileCopierWithProgressMeter}.  {@code false} otherwise.
     */
    public boolean isCopyComplete() {
        return copyComplete;
    }

    /**
     * @return {@code true} iff this {@link NioFileCopierWithProgressMeter} will overwrite {@link #dest} when copying.  {@code false} otherwise.
     */
    public boolean isOverwriteExisting() {
        return overwriteExisting;
    }

    /**
     * Sets whether this {@link NioFileCopierWithProgressMeter} will overwrite {@link #dest} when copying.
     */
    public NioFileCopierWithProgressMeter setOverwriteExisting(final boolean overwriteExisting) {
        this.overwriteExisting = overwriteExisting;

        return this;
    }

    /**
     * Sets the {@link #messageDigest} and expected checksum to be used by this {@link NioFileCopierWithProgressMeter} to validate the copied file.
     * NOTE: Setting these values will automatically cause the file to be checked for integrity upon copy completion.
     * @param algorithm {@link String} specifying the checksum algorithm to be used to calculate the checksum of the copied file.
     * @param expectedChecksum Expected value of the checksum calculated by the given {@code messageDigest} for the copied file.
     */
    public NioFileCopierWithProgressMeter setChecksumAlgorithmAndExpectedChecksum(final String algorithm,
                                                                                  final String expectedChecksum) {
        try {
            this.messageDigest = MessageDigest.getInstance(algorithm);
        }
        catch ( final NoSuchAlgorithmException ex ) {
            throw new IllegalArgumentException("Provided checksum algorithm does not exist: " + algorithm, ex);
        }
        this.expectedChecksum = expectedChecksum;

        return this;
    }

    /**
     * Sets the logger to log the time remaining in timestamp format ala 'DD:HH:MM:ss.SS'.
     */
    public NioFileCopierWithProgressMeter setFormatTimeRemainingAsTimestamp() {
        formatTimeRemainingAsTimestamp = true;

        return this;
    }

    /**
     * Sets the logger to log the time remaining in word format ala 'D days, H hours, M minutes, s seconds'.
     */
    public NioFileCopierWithProgressMeter setFormatTimeRemainingAsWords() {
        formatTimeRemainingAsTimestamp = false;

        return this;
    }

    /**
     * @return The {@link Verbosity} at which this {@link NioFileCopierWithProgressMeter} will log copy progress.
     */
    public Verbosity getVerbosity() {
        return verbosity;
    }

    /**
     * Sets the progress meter {@link #verbosity} of this {@link NioFileCopierWithProgressMeter}.
     */
    public NioFileCopierWithProgressMeter setVerbosity(final Verbosity verbosity) {
        this.verbosity = verbosity;

        return this;
    }

    //==================================================================================================================
    // Instance Methods:

    protected void updateMessageDigest(final byte[] copyBuffer, final int startIndex, final int endIndex) {
        if ( messageDigest != null ) {
            messageDigest.update(copyBuffer, startIndex, endIndex - startIndex);
        }
    }

    protected void calculateChecksumFromMessageDigest() {
        if ( messageDigest != null ) {
            checksum = DatatypeConverter.printHexBinary(messageDigest.digest());
        }
    }

    protected boolean isSilent() {
        return verbosity == Verbosity.SILENT;
    }

    protected String formatMillisecondsTime(final long time_ms) {
        if ( formatTimeRemainingAsTimestamp ) {
            return new AsTimeTimeFormatter(time_ms).format();
        }
        else {
            return new AsWordsTimeFormatter(time_ms).format();
        }
    }

    protected void logProgress(final double progressValue, final long totalBytesRead, final double bytesPerMillisecond) {
        if ( verbosity == Verbosity.VERBOSE ) {
            logProgressVerbose(progressValue, totalBytesRead, bytesPerMillisecond);
        }
        else if (verbosity.isAbove(Verbosity.MINIMAL)) {
            logProgressSimple(progressValue, totalBytesRead, bytesPerMillisecond);
        }
    }


    protected Duration getRemainingDuration(final long totalBytesRead, final double bytesPerMillisecond) {
        final long     remainingFileSize_bytes  = srcFileSize - totalBytesRead;
        final double   estTimeRemaining_ms      = remainingFileSize_bytes / bytesPerMillisecond;
        return Duration.ofMillis((long) estTimeRemaining_ms);
    }

    protected void logProgressSimple(final double progressValue, final long totalBytesRead, final double bytesPerMillisecond) {

        // Get the remaining time estimate:
        final Duration estTimeRemainingDuration = getRemainingDuration(totalBytesRead, bytesPerMillisecond);

        logger.info(
                String.format("    Transfer: % 2.2f%% complete.  Est. time remaining: %s (@%3.02f kbps)",
                        progressValue,
                        formatMillisecondsTime(estTimeRemainingDuration.toMillis()),
                        bytesPerMillisecond / KB_TO_BYTES * MS_TO_SEC
                )
        );
    }

    protected void logProgressVerbose(final double progressValue, final long totalBytesRead, final double bytesPerMillisecond) {

        // Get the remaining time estimate:
        final Duration estTimeRemainingDuration = getRemainingDuration(totalBytesRead, bytesPerMillisecond);

        final Instant endTime = Instant.now().plus(estTimeRemainingDuration);

        logger.info(
                String.format("    Transfer: % 2.2f%% complete (%" + srcFileSizeNumDigits + "d bytes; %6s).  Est. time remaining: %s (Complete time: %s) (@%3.02f kbps)",
                        progressValue,
                        totalBytesRead,
                        FileUtils.byteCountToDisplaySize(totalBytesRead),
                        formatMillisecondsTime(estTimeRemainingDuration.toMillis()),
                        endTime.atZone(ZoneId.systemDefault()).toLocalDateTime().toString(),
                        bytesPerMillisecond / KB_TO_BYTES * MS_TO_SEC
                )
        );
    }

    private void initializeCopyProgressTime(final long startTime_ns) {

        // Track the time it takes to download each chunk:
        lastProgressTime_ns = startTime_ns;
    }

    protected void updateCopyProgress(final int bytesRead) {

        // Only bother with all this if we're logging in the first place:
        if ( !isSilent() ) {

            // Update our progress counters:
            totalBytesRead += bytesRead;
            progressBytesRead += bytesRead;

            // Get our progress percentage:
            final double rawProgressValuePercent = ((double) totalBytesRead / (double) srcFileSize) * 100.0;

            // Round our progress to nearest PROGRESS_DISPLAY_PERCENT_INCREMENT:
            final double progressValue = progressPercentDisplayIncrement * (Math.floor(Math.abs(rawProgressValuePercent / progressPercentDisplayIncrement)));

            // Output our progress if we're ready for it:
            if ( progressValue != lastProgressValue ) {

                // Update our time:
                final long   currentProgressTime_ns = System.nanoTime();
                final long   dt_ms                  = (currentProgressTime_ns - lastProgressTime_ns) / NANOS_TO_MILLIS;
                final double bytesPerMs             = ((double) progressBytesRead) / ((double) dt_ms);
                lastProgressTime_ns = currentProgressTime_ns;

                // Add the bytes/ms to our queue:
                downloadBytesPerMilliSecond.add(bytesPerMs);

                // Log our progress so far:
                logProgress(progressValue, totalBytesRead, downloadBytesPerMilliSecond.stream().mapToDouble(x -> x).average().orElse(dt_ms));

                // Get ready for the next call:
                lastProgressValue = progressValue;
                progressBytesRead = 0;
            }
        }
    }

    protected void determineProgessDisplayIncrement(final long fileSize) {
        // Simple checks on input file size to make sure we don't overwhelm or underwhelm the user with updates:
        // TODO: Refactor class to have download and logger in separate threads and display on every percentage complete or delta-time.

        final long SIZE_STEP = 1024;
        final long KB        = 1024;
        final long MB        = KB * SIZE_STEP;
        final long GB        = MB * SIZE_STEP;

        // 100Gb or larger:
        if ( fileSize >= (100*GB) ) {
            progressPercentDisplayIncrement = 0.1;
        }
        // 10Gb or larger:
        else if ( fileSize >= (10*GB) ) {
            progressPercentDisplayIncrement = 0.25;
        }
        // 5Gb or larger:
        else if ( fileSize >= (5*GB) ) {
            progressPercentDisplayIncrement = 0.5;
        }
        // 1Gb or larger:
        else if ( fileSize >= GB ) {
            progressPercentDisplayIncrement = 1;
        }
        // 100Mb or larger:
        else if ( fileSize >= (MB*100) ) {
            progressPercentDisplayIncrement = 5;
        }
        // 1Mb or larger:
        else if ( fileSize >= MB ) {
            progressPercentDisplayIncrement = 10;
        }
        // Less than 1Mb
        else {
            progressPercentDisplayIncrement = 25;
        }
    }

    protected void doCopy() {

        try ( final InputStream inputStream = Files.newInputStream(getSource());
              final OutputStream outputStream = Files.newOutputStream(getDest()) ){

            // Get the file size of our source file:
            srcFileSize = Files.size(getSource());
            srcFileSizeNumDigits = (int)Math.ceil(Math.log10(srcFileSize));

            determineProgessDisplayIncrement(srcFileSize);

            if ( verbosity.isAbove(Verbosity.SILENT) ) {
                logger.info("Initiating copy from " + getSource().toUri().toString() + " to " + getDest().toUri().toString());
                logger.info("File size: " + srcFileSize + " bytes (" + FileUtils.byteCountToDisplaySize(srcFileSize) + ").");
                logger.info("Please wait.  This could take a while...");
            }

            // Perform the copy:
            while (true) {

                // Read from our input:
                final int bytesRead = inputStream.read(copyBuffer);
                if ( bytesRead == -1 ) {
                    break;
                }

                // Write to our output:
                outputStream.write(copyBuffer, 0, bytesRead);

                // Update the message digest so we can calculate the file checksum on the fly:
                updateMessageDigest(copyBuffer, 0, bytesRead);

                // Update our internal progress meter:
                updateCopyProgress(bytesRead);
            }

            // Calculate the checksum from the message digest:
            calculateChecksumFromMessageDigest();
        }
        catch (final IOException ex) {
            throw new UserException("Could not copy file: " + source.toUri().toString() + " -> " + dest.toUri().toString(), ex);
        }
    }

    /**
     * Initiate the copy from {@link #source} to {@link #dest}.
     */
    public NioFileCopierWithProgressMeterResults initiateCopy() {

        // Make sure we haven't copied the file already:
        if (copyComplete) {
            throw new GATKException("Attempted multiple file copies.  NioFileCopierWithProgressMeter can copy a file only once!");
        }

        // Do a quick existence check for safety:
        if ( Files.exists(getDest()) ) {
            if ( !isOverwriteExisting() ) {
                throw new UserException.CouldNotCreateOutputFile(getDest().toUri().toString(), "Download aborted!  Output data sources file already exists!");
            }
            else if ( verbosity.isAbove(Verbosity.SILENT) ) {
                logger.warn("Destination already exists.  Overwriting file at location: " + getDest().toUri().toString());
            }
        }

        // Keep track of the start time so we can estimate completion time:
        final long startTime_ns = System.nanoTime();

        // Initialize our internal progress meter:
        initializeCopyProgressTime(startTime_ns);

        // Now copy from our source to our destination:
        doCopy();

        // Let the world know the glory that is a complete file copy:
        if ( verbosity.isAbove(Verbosity.SILENT) ) {
            logger.info(String.format("Download Complete! - Total elapsed time: %ds", ((System.nanoTime() - startTime_ns) / NANOS_TO_SECONDS)));
        }

        // Make sure we don't copy this file more than once:
        copyComplete = true;

        return new NioFileCopierWithProgressMeterResults(
                source,
                dest,
                srcFileSize,
                messageDigest != null,
                checksum,
                messageDigest == null ? "" : messageDigest.getAlgorithm(),
                expectedChecksum
        );
    }

    //==================================================================================================================
    // Helper Data Types:

    /**
     * An enum to allow for verbosity of logging progress of an {@link NioFileCopierWithProgressMeter}.
     */
    public enum Verbosity {
        /**
         * Output no logging messages whatsoever.
         */
        SILENT(0),
        /**
         * Output logging messages at the start and end of the copy, but no progress during.
         */
        MINIMAL(1),
        /**
         * Output basic progress information during the copy.
         */
        MODERATE(2),
        /**
         * Output verbose progress information during the copy.
         */
        VERBOSE(3);

        final private int sev;

        Verbosity(final int sev) { this.sev = sev; }

        public boolean isAbove(final Verbosity other) {
            return this.sev > other.sev;
        }
    }

    /**
     * An interface that defines a method to use to calculate a checksum on an {@link InputStream}.
     * Used to verify file contents are correct and have not been corrupted in-transit.
     */
    public interface ChecksumCalculator {
        String calculateChecksumOnInputStream(InputStream data) throws IOException;
    }

    /**
     * Simple class to keep track of time information and format it.
     */
    private abstract class SimpleTimeFormatter {

        final long rawTime_ms;
        final long days;
        final long hours;
        final long minutes;
        final long seconds;
        final long millis;

        SimpleTimeFormatter(final long time_ms) {
            rawTime_ms = time_ms;

            long remainder = time_ms;

            days = formatTimeHelper(remainder, DAY_IN_MS);
            remainder -= days * DAY_IN_MS;

            hours = formatTimeHelper(remainder, HOUR_IN_MS);
            remainder -= hours * HOUR_IN_MS;

            minutes = formatTimeHelper(remainder, MINUTE_IN_MS);
            remainder -= minutes * MINUTE_IN_MS;

            seconds = formatTimeHelper(remainder, SECOND_IN_MS);
            remainder -= seconds * SECOND_IN_MS;

            millis = remainder;
        }

        private long formatTimeHelper(final long duration, final long conversionFactor ) {
            final long outTime;
            if ( duration > conversionFactor ) {
                outTime = Math.floorDiv(duration, conversionFactor);
            }
            else {
                outTime = 0;
            }

            return outTime;
        }

        protected String sHelper(final long value) {
            return (value == 1 ? "" : "s");
        }

        public abstract String format();
    }

    private class AsWordsTimeFormatter extends SimpleTimeFormatter {

        AsWordsTimeFormatter(final long time_ms){
            super(time_ms);
        }

        public String format() {
            if ( days > 0 ) {
                return String.format("%d day" + sHelper(days) + ", %02d hour" + sHelper(hours) + ", %02d minute" + sHelper(minutes) + ", %2d.%03d seconds", days, hours, minutes, seconds, millis);
            }
            if ( hours > 0 ) {
                return String.format("%02d hour" + sHelper(hours) + ", %02d minute" + sHelper(minutes) + ", %02d.%03d seconds", hours, minutes, seconds, millis);
            }
            if ( minutes > 0 ) {
                return String.format("%02d minute" + sHelper(minutes) + ", %02d.%03d seconds", minutes, seconds, millis);
            }
            if ( seconds > 0 ) {
                return String.format("%02d.%03d seconds", seconds, millis);
            }
            return String.format("0.%03d seconds", millis);
        }
    }
    private class AsTimeTimeFormatter extends SimpleTimeFormatter {

        AsTimeTimeFormatter(final long time_ms){
            super(time_ms);
        }
        public String format() {

            if ( days > 0 ) {
                return String.format("%d:%02d:%02d:%02d.%03d", days, hours, minutes, seconds, millis);
            }
            if ( hours > 0 ) {
                return String.format("%02d:%02d:%02d.%03d", hours, minutes, seconds, millis);
            }
            if ( minutes > 0 ) {
                return String.format("%02d:%02d.%03d", minutes, seconds, millis);
            }
            if ( seconds > 0 ) {
                return String.format("%02d.%03d", seconds, millis);
            }
            return String.format("00.%03d", millis);
        }
    }
}
