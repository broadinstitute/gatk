package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.ValidationStringency;

/**
 * Constants for use with the GATKRead interface
 */
public final class ReadConstants {
    private ReadConstants(){}

    /**
     * Value used as the default validation stringency for all read input
     */
    public static ValidationStringency DEFAULT_READ_VALIDATION_STRINGENCY = ValidationStringency.SILENT;

    /**
     * Value used to represent the absence of a defined start/end position in a read
     */
    public static final int UNSET_POSITION = 0;

    /**
     * Value used to represent the absence of a defined contig in a read
     */
    public static final String UNSET_CONTIG = "*";

    /**
     * Value used to represent the absence of a mapping quality in a read
     */
    public static final int NO_MAPPING_QUALITY = 0;

    /**
     * Value used for String representation of a non-existent sequence
     */
    public static final String NULL_SEQUENCE_STRING = "*";
}
