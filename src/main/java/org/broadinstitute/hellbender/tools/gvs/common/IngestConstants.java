package org.broadinstitute.hellbender.tools.gvs.common;

public class IngestConstants {

    public static final char SEPARATOR = '\t';
    public static final String FILETYPE = ".tsv";
    public static final String sampleMetadataFilePrefix = "sample_"; // used for arrays
    public static final String sampleInfoFilePrefix = "sample_info_";
    public static final int partitionPerTable = 4000;
    public static final int MAX_REFERENCE_BLOCK_BASES = 1000;
    public static final int MAX_DELETION_SIZE = 1000;
}
