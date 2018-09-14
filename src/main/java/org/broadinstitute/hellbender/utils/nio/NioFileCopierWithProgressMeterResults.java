package org.broadinstitute.hellbender.utils.nio;

import java.nio.file.Path;

/**
 * An object to hold the results of a copy operation performed by {@link NioFileCopierWithProgressMeterResults}.
 * Designed to be a simple data object that holds information, not to be instantiated by a user directly.
 * Created by jonn on 9/4/18.
 */
public class NioFileCopierWithProgressMeterResults {

    //==================================================================================================================
    // Private Members:

    private final Path    src;
    private final Path    dest;
    private final long    size;
    private final boolean wasValidationRequested;
    private final String  checksum;
    private final String  checksumAlgorithm;
    private final String  expectedChecksum;


    //==================================================================================================================
    // Constructors:

    /**
     * Package-private constructor.
     * Designed never to be instantiated directly.
     * Should only be created by {@link NioFileCopierWithProgressMeter}.
     */
    NioFileCopierWithProgressMeterResults(
            final Path    src,
            final Path    dest,
            final long    size,
            final boolean checksumWasCalculated,
            final String  checksum,
            final String  checksumAlgorithm,
            final String  expectedChecksum) {

        this.src                   = src;
        this.dest                  = dest;
        this.size                  = size;
        this.wasValidationRequested = checksumWasCalculated;
        this.checksum              = checksum;
        this.checksumAlgorithm     = checksumAlgorithm;
        this.expectedChecksum      = expectedChecksum;
    }

    //==================================================================================================================
    // Instance Methods:

    /**
     * @return The {@link Path} of the source file to copy.
     */
    public Path getSource() {
        return src;
    }

    /**
     * @return The {@link Path} to the resulting destination of the file after copy.
     */
    public Path getDestination() {
        return dest;
    }

    /**
     * @return The size of the file downloaded.
     */
    public long getSize() {
        return size;
    }

    /**
     * Get whether a user requested that a checksum be calculated on the downloaded file.
     * @return {@code true} iff a checksum was calculated by during the transfer of the file.  {@code false} otherwise.
     */
    public boolean wasValidationRequested() {
        return wasValidationRequested;
    }

    /**
     * Get the checksum calculated against {@link #dest}.
     * If no checksum has been calculated, will return an empty {@link String}.
     * @return The last checksum calculated against {@link #dest}.
     */
    public String getChecksum() {
        return checksum;
    }

    /**
     * Get the checksum algorithm used to calculate the checksum of the downloaded file.
     * If no checksum algorithm was used, will return an empty {@link String}.
     * @return The checksum algorithm used to calculate the checksum of the downloaded file.
     */
    public String getChecksumAlgorithm() {
        return checksumAlgorithm;
    }

    /**
     * Get the expected checksum value for {@link #dest} as set by the user in {@link NioFileCopierWithProgressMeter#setChecksumAlgorithmAndExpectedChecksum(String, String)}.
     * If no checksum has been set, will return an empty {@link String}.
     * @return The expected checksum for {@link #dest} provided by the user in {@link NioFileCopierWithProgressMeter#setChecksumAlgorithmAndExpectedChecksum(String, String)}.
     */
    public String getExpectedChecksum() {
        return expectedChecksum;
    }

    /**
     * Returns whether the downloaded file is valid by comparing the {@link #checksum} and {@link #expectedChecksum} values.
     * Will perform a case-insensitive comparison between the values of {@link #checksum} and {@link #expectedChecksum}.
     * NOTE: If no checksum algorithm and expected value were specified in the {@link NioFileCopierWithProgressMeterResults},
     * then this method will return {@code false}.
     * @return {@code true} if the checksum from the downloaded file matches the expected checksum previously given by the user.
     */
    public boolean isDestFileValid() {
        if ( wasValidationRequested  ) {
            return checksum.toLowerCase().equals(expectedChecksum.toLowerCase());
        }
        else {
            return false;
        }
    }

}
