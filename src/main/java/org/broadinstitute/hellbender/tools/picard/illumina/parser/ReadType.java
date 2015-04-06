package org.broadinstitute.hellbender.tools.picard.illumina.parser;

/**
 * A read type describes a stretch of cycles in an ReadStructure
 * (e.g. Assume we have a paired end/barcoded run with the 76 template cycles followed by 8 barcode cycles followed by
 * another 76 template reads, the run would be represented by the ReadStructure 76T8B76T)
 * Note: Currently SKIP is unused by IlluminaBasecallsToSam, ExtractIlluminaBarcodes, and IlluminaDataProvider
 */
public enum ReadType {
    T, B, S;

    public static final ReadType Template = T;
    public static final ReadType Barcode = B;
    public static final ReadType Skip = S;
}
