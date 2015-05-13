package org.broadinstitute.hellbender.utils.reference;

/**
 * A collection of static methods for dealing with references.
 */
public final class ReferenceUtils {

    // Private so that no one will instantiate this class.
    private ReferenceUtils() {}

    /**
     * Given a fasta filename, return the name of the corresponding index file.
     * (This also works if the file is in gs://)
     */
    public static String getFastaIndexFileName(String fastaFilename) {
        return fastaFilename + ".fai";
    }

    /**
     * Given a fasta filename, return the name of the corresponding dictionary file.
     * (This also works if the file is in gs://)
     */
    public static String getFastaDictionaryFileName(String fastaFilename) {
        int lastDot = fastaFilename.lastIndexOf('.');
        return fastaFilename.substring(0, lastDot) + ".dict";
    }


}
