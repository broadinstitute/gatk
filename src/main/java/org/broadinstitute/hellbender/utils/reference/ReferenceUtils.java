package org.broadinstitute.hellbender.utils.reference;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.BufferedLineReader;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

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

    /**
     * Given a fasta dictionary file, returns its sequence dictionary
     *
     * @param fastaDictionaryFile fasta dictionary file
     * @return the SAMSequenceDictionary from fastaDictionaryFile
     */
    public static SAMSequenceDictionary loadFastaDictionary( final File fastaDictionaryFile ) {
        try ( final FileInputStream fastaDictionaryStream = new FileInputStream(fastaDictionaryFile) ) {
            return loadFastaDictionary(fastaDictionaryStream);
        }
        catch ( IOException e ) {
            throw new UserException.CouldNotReadInputFile("Error loading fasta dictionary file " + fastaDictionaryFile, e);
        }
        catch ( UserException.MalformedFile e ) {
            throw new UserException.MalformedFile(
                    "Could not read sequence dictionary from given fasta file " +
                            fastaDictionaryFile
            );
        }
    }

    /**
     * Given an InputStream connected to a fasta dictionary, returns its sequence dictionary
     *
     * Note: does not close the InputStream it's passed
     *
     * @param fastaDictionaryStream InputStream connected to a fasta dictionary
     * @return the SAMSequenceDictionary from the fastaDictionaryStream
     */
    public static SAMSequenceDictionary loadFastaDictionary( final InputStream fastaDictionaryStream ) {
        // Don't close the reader when we're done, since we don't want to close the client's InputStream for them
        final BufferedLineReader reader = new BufferedLineReader(fastaDictionaryStream);

        final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
        final SAMFileHeader header = codec.decode(reader, fastaDictionaryStream.toString());

        // Make sure we have a valid sequence dictionary before continuing:
        if (header.getSequenceDictionary() == null || header.getSequenceDictionary().isEmpty()) {
            throw new UserException.MalformedFile (
                    "Could not read sequence dictionary from given fasta stream " +
                            fastaDictionaryStream
            );
        }

        return header.getSequenceDictionary();
    }
}
