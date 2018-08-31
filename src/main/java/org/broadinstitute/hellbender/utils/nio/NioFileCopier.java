package org.broadinstitute.hellbender.utils.nio;

import org.apache.commons.collections4.queue.CircularFifoQueue;
import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.time.Duration;
import java.time.Instant;
import java.time.ZoneId;
import java.util.HashMap;
import java.util.Map;
import java.util.Queue;

/**
 * Class to copy a file using {@link java.nio}.
 * Operates using paths.
 *
 * INSTANCES OF THIS CLASS ARE NOT THREAD-SAFE!
 *
 * Created by jonn on 8/27/18.
 */
public class NioFileCopier {

    //==================================================================================================================
    // Standard logger:
    private static final Logger logger = LogManager.getLogger(NioFileCopier.class);

    //==================================================================================================================
    // Public Static Members:

    //==================================================================================================================
    // Private Static Members:

    private static final int    BUFFER_SIZE_BYTES                  = 1024 * 1024;
    private static final double PROGRESS_DISPLAY_PERCENT_INCREMENT = 0.25;

    private static final long SECOND_IN_MS = 1000;
    private static final long MINUTE_IN_MS = SECOND_IN_MS * 60;
    private static final long HOUR_IN_MS   = MINUTE_IN_MS * 60;
    private static final long DAY_IN_MS    = HOUR_IN_MS * 24;

    private static final int KB_TO_BYTES = 1024;
    private static final int MS_TO_SEC = 1000;
    private static final int NANOS_TO_MILLIS  = 1000000;
    private static final int NANOS_TO_SECONDS = 1000000000;

    private static final int COPY_SPEED_HISTORY_SIZE = 10;

    private static final boolean OVERWRITE_EXISTING_DEFAULT = false;
    private static final boolean DO_LOG_PROGRESS_DEFAULT    = true;
    private static final boolean DO_LOG_VERBOSE_DEFAULT     = false;

    //==================================================================================================================
    // Private Members:

    // Data variables:
    protected final Path source;
    protected final Path dest;

    private long srcFileSize;
    private int srcFileSizeNumDigits;

    private Map<ChecksumCalculator, String> checksumMap              = new HashMap<>();
    private String                          latestChecksum           = "";
    private ChecksumCalculator              latestChecksumCalculator = null;

    // Flag defaults:
    private boolean overwriteExisting = OVERWRITE_EXISTING_DEFAULT;
    private boolean doLogProgress     = DO_LOG_PROGRESS_DEFAULT;
    private boolean logVerbose        = DO_LOG_VERBOSE_DEFAULT;
    private boolean silentCopy        = false;

    // Copy buffer:
    private final byte          copyBuffer[]                = new byte[ BUFFER_SIZE_BYTES ];

    // Progress variables:
    private final Queue<Double> downloadBytesPerMilliSecond = new CircularFifoQueue<>(COPY_SPEED_HISTORY_SIZE);
    private       boolean       copyComplete                = false;
    private       long          totalBytesRead              = 0;
    private       long          progressBytesRead           = 0;
    private       double        lastProgressValue           = 0;
    private long lastProgressTime_ns;

    //==================================================================================================================
    // Constructors:

    /**
     * {@link NioFileCopier} uses a factory pattern.
     * This internal constructor is to be used by the class itself.
     * @param source The {@link Path} to the source file for the copy.
     * @param dest The {@link Path} to the destination file for the copy.
     * @param overwriteExisting If {@code true} will overwrite an existing file in the location specified by {@code dest}.
     * @param doLogProgress If {@code true} will log progress over the duration of the copy.
     */
    private NioFileCopier(final Path source, final Path dest, final boolean overwriteExisting, final boolean doLogProgress) {
        this.source = source;
        this.dest = dest;
        this.overwriteExisting = overwriteExisting;
        this.doLogProgress = doLogProgress;
    }

    //==================================================================================================================
    // Static Methods:

    /**
     * Create an {@link NioFileCopier}.
     * By default the resulting {@link NioFileCopier} will not overwrite the destination if anything already exists there.
     * @param source The {@link Path} to the source file for the copy.
     * @param dest The {@link Path} to the destination file for the copy.
     * @return An {@link NioFileCopier} initialized to copy the file located at {@code source} to the location specified by {@code dest}.
     */
    public static NioFileCopier create(final Path source, final Path dest) {
        return create(source, dest, OVERWRITE_EXISTING_DEFAULT);
    }

    /**
     * Create an {@link NioFileCopier}.
     * Will periodically display progress of copying files.
     * @param source The {@link Path} to the source file for the copy.
     * @param dest The {@link Path} to the destination file for the copy.
     * @param overwriteExisting If {@code true} will overwrite an existing file in the location specified by {@code dest}.
     * @return An {@link NioFileCopier} initialized to copy the file located at {@code source} to the location specified by {@code dest}.
     */
    public static NioFileCopier create(final Path source, final Path dest, final boolean overwriteExisting) {
        return new NioFileCopier(source, dest, overwriteExisting, DO_LOG_PROGRESS_DEFAULT);
    }

    /**
     * Create an {@link NioFileCopier}.
     * @param source The {@link Path} to the source file for the copy.
     * @param dest The {@link Path} to the destination file for the copy.
     * @param overwriteExisting If {@code true} will overwrite an existing file in the location specified by {@code dest}.
     * @param doLogProgress If {@code true} will log progress over the duration of the copy.
     * @return An {@link NioFileCopier} initialized to copy the file located at {@code source} to the location specified by {@code dest}.
     */
    public static NioFileCopier create(final Path source, final Path dest, final boolean overwriteExisting, final boolean doLogProgress) {
        return new NioFileCopier(source, dest, overwriteExisting, doLogProgress);
    }

    //==================================================================================================================
    // Getters / Setters:

    /**
     * @return A copy of the {@link Path} used as the source for this {@link NioFileCopier}.
     */
    public Path getSource() {
        return Paths.get(source.toUri());
    }

    /**
     * @return A copy of the {@link Path} used as the destination for this {@link NioFileCopier}.
     */
    public Path getDest() {
        return Paths.get(dest.toUri());
    }

    /**
     * @return {@code true} iff the {@link #source} has been copied to the {@link #dest} in this {@link NioFileCopier}.  {@code false} otherwise.
     */
    public boolean isCopyComplete() {
        return copyComplete;
    }

    /**
     * @return {@code true} iff this {@link NioFileCopier} will overwrite {@link #dest} when copying.  {@code false} otherwise.
     */
    public boolean isOverwriteExisting() {
        return overwriteExisting;
    }

    /**
     * Sets whether this {@link NioFileCopier} will overwrite {@link #dest} when copying.
     */
    public NioFileCopier setOverwriteExisting(final boolean overwriteExisting) {
        this.overwriteExisting = overwriteExisting;

        return this;
    }

    /**
     * @return {@code true} iff this {@link NioFileCopier} will log intermediate progress when copying.  {@code false} otherwise.
     */
    public boolean isDoLogProgress() {
        return doLogProgress;
    }

    /**
     * Sets whether this {@link NioFileCopier} will log intermediate progress when copying.
     */
    public NioFileCopier setDoLogProgress(final boolean doLogProgress) {
        this.doLogProgress = doLogProgress;

        return this;
    }

    /**
     * {@code true} if this {@link NioFileCopier} will verbosely log intermediate progress when copying.
     * This is independent of whether intermediate progress logging is enabled overall.
     * Intermediate progress logging is mediated by {@link #doLogProgress}.
     * @return {@code true} iff this {@link NioFileCopier} will verbosely log intermediate progress when copying.  {@code false} otherwise.
     */
    public boolean isLogVerbose() {
        return logVerbose;
    }

    /**
     * Sets whether this {@link NioFileCopier} will verbosely log intermediate progress when copying.
     * This is independent of whether intermediate progress logging is enabled overall.
     * Intermediate progress logging is mediated by {@link #doLogProgress}.
     */
    public NioFileCopier setLogVerbose(final boolean logVerbose) {
        this.logVerbose = logVerbose;

        return this;
    }

    /**
     * @return {@code true} iff this {@link NioFileCopier} will copy completely silently, with no logging output whatsoever.  {@code false} otherwise.
     */
    public boolean isSilentCopy() {
        return silentCopy;
    }

    /**
     * Sets whether this {@link NioFileCopier} will copy completely silently, with no logging output whatsoever.
     */
    public NioFileCopier setSilentCopy(final boolean silentCopy) {
        this.silentCopy = silentCopy;

        return this;
    }

    /**
     * Get the last checksum calculated against {@link #dest}.
     * If no checksum has been calculated, will return an empty {@link String}.
     * @return The last checksum calculated against {@link #dest}.
     */
    public String getLatestChecksum() {
        return latestChecksum;
    }

    /**
     * Get the last {@link ChecksumCalculator} used to calculate the checksum of {@link #dest}.
     * If no checksum has been calculated, will return {@code null}.
     * @return The last {@link ChecksumCalculator} used to calculate the checksum of {@link #dest}.
     */
    public ChecksumCalculator getLatestChecksumCalculator() {
        return latestChecksumCalculator;
    }

    //==================================================================================================================
    // Instance Methods:

    private String formatMillisecondsTime(final long time_ms) {
//        return new AsWordsTimeFormatter(time_ms).format();
        return new AsTimeTimeFormatter(time_ms).format();
    }

    private String sHelper(final long value) {
        return (value == 1 ? "" : "s");
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

    private void logProgress(final double progressValue, final long totalBytesRead, final double bytesPerMillisecond) {
        if ( logVerbose ) {
            logProgressVerbose(progressValue, totalBytesRead, bytesPerMillisecond);
        }
        else {
            logProgressSimple(progressValue, totalBytesRead, bytesPerMillisecond);
        }
    }

    private void logProgressSimple(final double progressValue, final long totalBytesRead, final double bytesPerMillisecond) {

        // Get the remaining time estimate:
        final long     remainingFileSize_bytes  = srcFileSize - totalBytesRead;
        final double   estTimeRemaining_ms      = remainingFileSize_bytes / bytesPerMillisecond;
        final Duration estTimeRemainingDuration = Duration.ofMillis((long) estTimeRemaining_ms);

        logger.info(
                String.format("    Transfer: % 2.2f%% complete.  Est. time remaining: %s (@%3.02f kbps)",
                        progressValue,
                        formatMillisecondsTime(estTimeRemainingDuration.toMillis()),
                        bytesPerMillisecond / KB_TO_BYTES * MS_TO_SEC
                )
        );
    }

    private void logProgressVerbose(final double progressValue, final long totalBytesRead, final double bytesPerMillisecond) {

        // Get the remaining time estimate:
        final long     remainingFileSize_bytes  = srcFileSize - totalBytesRead;
        final double   estTimeRemaining_ms      = remainingFileSize_bytes / bytesPerMillisecond;
        final Duration estTimeRemainingDuration = Duration.ofMillis((long) estTimeRemaining_ms);

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

    private void initializeCopyProgress(final long startTime_ns) {

        // Track the time it takes to download each chunk:
        lastProgressTime_ns = startTime_ns;

        // Keep track of the total number of bytes read and our progress value:
        totalBytesRead = 0;
        progressBytesRead = 0;
        lastProgressValue = 0;
    }

    private void updateCopyProgress(final int bytesRead) {

        // Only bother with all this if we're logging in the first place:
        if ( (!isSilentCopy()) && isDoLogProgress() ) {

            // Update our progress counters:
            totalBytesRead += bytesRead;
            progressBytesRead += bytesRead;

            // Get our progress percentage:
            final double rawProgressValuePercent = ((double) totalBytesRead / (double) srcFileSize) * 100.0;

            // Round our progress to nearest PROGRESS_DISPLAY_PERCENT_INCREMENT:
            final double progressValue = PROGRESS_DISPLAY_PERCENT_INCREMENT * (Math.floor(Math.abs(rawProgressValuePercent / PROGRESS_DISPLAY_PERCENT_INCREMENT)));

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

    private void copyLoopWorker() {

        try ( final InputStream inputStream = Files.newInputStream(getSource());
              final OutputStream outputStream = Files.newOutputStream(getDest()) ){

            // Get the file size of our source file:
            srcFileSize = Files.size(getSource());
            srcFileSizeNumDigits = (int)Math.ceil(Math.log10(srcFileSize));

            if ( !isSilentCopy() ) {
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

                // Update our internal progress meter:
                updateCopyProgress(bytesRead);
            }
        }
        catch (final IOException ex) {
            throw new GATKException("Could not copy file: " + source.toUri().toString() + " -> " + dest.toUri().toString(), ex);
        }
    }

    /**
     * Initiate the copy from {@link #source} to {@link #dest}.
     */
    public void initiateCopy() {

        // Make sure we haven't copied the file already:
        if (copyComplete) {
            throw new GATKException("Attempted multiple file copies.  NioFileCopier can copy a file only once!");
        }

        // Do a quick existence check for safety:
        if ( Files.exists(getDest()) ) {
            if ( !isOverwriteExisting() ) {
                throw new UserException.CouldNotCreateOutputFile(getDest().toUri().toString(), "Download aborted!  Output data sources file already exists!");
            }
            else if ( !isSilentCopy() ) {
                logger.warn("Destination already exists.  Overwriting file at location: " + getDest().toUri().toString());
            }
        }

        // Keep track of the start time so we can estimate completion time:
        final long startTime_ns = System.nanoTime();

        // Initialize our internal progress meter:
        initializeCopyProgress(startTime_ns);

        // Now copy from our source to our destination:
        copyLoopWorker();

        // Let the world know the glory that is a complete file copy:
        if ( !isSilentCopy() ) {
            logger.info(String.format("Download Complete! - Total elapsed time: %ds", ((System.nanoTime() - startTime_ns) / NANOS_TO_SECONDS)));
        }

        // Make sure we don't copy this file more than once:
        copyComplete = true;
    }

    /**
     * Validate the integrity of the downloaded file against the given {@code expectedChecksum} using the given {@link ChecksumCalculator}.
     * Must be called after the copy is complete (i.e. when {@link #copyComplete} = {@code true}.
     * Will not recalculate the checksum if it has already been calculated for the given {@link ChecksumCalculator}.
     * @param expectedChecksum A {@link String} containing the expected value of the checksum for {@link #source}.
     * @param checksumCalculator A method to calculate the checksum for the downloaded file, {@link #dest}.
     * @return {@code true} iff the given {@code expectedChecksum} matches the calculated checksum for {@link #dest} using {@code checksumCalculator}.  {@code false} otherwise.
     */
    public boolean validateIntegrity(final String expectedChecksum, final ChecksumCalculator checksumCalculator) {
        return validateIntegrity(expectedChecksum, checksumCalculator, false);
    }

    /**
     * Validate the integrity of the downloaded file against the given {@code expectedChecksum} using the given {@link ChecksumCalculator}.
     * Must be called after the copy is complete (i.e. when {@link #copyComplete} = {@code true}.
     * @param expectedChecksum A {@link String} containing the expected value of the checksum for {@link #source}.
     * @param checksumCalculator A method to calculate the checksum for the downloaded file, {@link #dest}.
     * @param forceRecalculate If {@code true} will recalculate the checksum for {@link #dest} even if it has already been calculated for the given {@link ChecksumCalculator}.
     * @return {@code true} iff the given {@code expectedChecksum} matches the calculated checksum for {@link #dest} using {@code checksumCalculator}.  {@code false} otherwise.
     */
    public boolean validateIntegrity(final String expectedChecksum, final ChecksumCalculator checksumCalculator, final boolean forceRecalculate) {

        // Can we actually do this yet?
        if ( !isCopyComplete() ) {
            throw new GATKException("Can only validate file integrity after the file copy is complete!");
        }

        final String actualChecksum;

        // If we have calculated this before, we should just return the cached version:
        if ( (!checksumMap.containsKey(checksumCalculator)) || forceRecalculate ) {

            if ( !Files.exists(dest) ) {
                throw new GATKException("File no longer exists.  Cannot calculate checksum for: " + dest.toUri());
            }

            // Calculate the checksum of dest:
            try ( final InputStream dataStream = Files.newInputStream(dest, StandardOpenOption.READ) ) {
                if ( !isSilentCopy() ) {
                    logger.info("Calculating checksum...");
                }

                actualChecksum = checksumCalculator.calculateChecksumOnInputStream(dataStream).trim().toLowerCase();

                if ( !isSilentCopy() ) {
                    logger.info("Calculation complete!");
                }
            }
            catch ( final IOException ex ) {
                throw new GATKException("Could not read destination file to calculate hash: " + dest.toUri().toString(), ex);
            }

            // Store our checksum for later:
            checksumMap.put(checksumCalculator, actualChecksum);
            latestChecksum = actualChecksum;
            latestChecksumCalculator = checksumCalculator;
        }
        else {
            // Grab the checksum from the cached map:
            actualChecksum = checksumMap.get(checksumCalculator);

            // Update the latest checksum we have:
            latestChecksum = actualChecksum;
            latestChecksumCalculator = checksumCalculator;
        }

        // Verify the checksums are the same:
        if ( !expectedChecksum.equals(actualChecksum) ) {
            if ( !isSilentCopy() ) {
                logger.warn("Warning: destination file is corrupt!  Unexpected checksum: " + actualChecksum + " != " + expectedChecksum);
            }
            return false;
        }
        else {
            if ( !isSilentCopy() ) {
                logger.info("Destination file is valid.");
            }
            return true;
        }
    }

    //==================================================================================================================
    // Helper Data Types:

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
