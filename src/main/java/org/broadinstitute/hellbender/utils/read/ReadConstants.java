package org.broadinstitute.hellbender.utils.read;

/**
 * Constants for use with the GATKRead interface
 */
public final class ReadConstants {

    /**
     * Value used to represent the absence of a defined start/end position in a read
     */
    public static final int UNSET_POSITION = 0;

    /**
     * Value used to represent the absence of a mapping quality in a read
     */
    public static final int NO_MAPPING_QUALITY = 0;

    /**
     * Value used for String representation of a non-existent sequence
     */
    public static final String NULL_SEQUENCE_STRING = "*";
}
