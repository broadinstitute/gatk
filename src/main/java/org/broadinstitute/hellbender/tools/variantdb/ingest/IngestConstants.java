package org.broadinstitute.hellbender.tools.variantdb.ingest;

public class IngestConstants {

    public static final char SEPARATOR = '\t';
    public static final String FILETYPE = ".tsv";
    public static final String metadataFilePrefix = "metadata_";
    public static final String metadataDirectoryName = "metadata"; // TODO remove
    public static final int partitionPerTable = 4000;

}
