package org.broadinstitute.hellbender.utils.codecs.gtf;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.readers.LineIterator;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/**
 * Codec to decode data in GTF format from ENSEMBL.
 * According to ENSEMBL, GTF files downloaded from them conform to GFF version 2 (http://gmod.org/wiki/GFF2).
 */
final public class EnsemblGtfCodec extends AbstractGtfCodec {

    private static final Logger logger = LogManager.getLogger(EnsemblGtfCodec.class);

    //==================================================================================================================
    // Public Static Members:

    public static String GTF_FILE_TYPE_STRING = "ENSEMBL";

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Private Members:

    private final List<String> header         = new ArrayList<>();
    private       int          currentLineNum = 1;
    private       String       version        = "";

    //==================================================================================================================
    // Constructors:

    public EnsemblGtfCodec() {
        super();
    }

    //==================================================================================================================
    // Override Methods:

    @Override
    String getGtfFileType() {
        return GTF_FILE_TYPE_STRING;
    }

    @Override
    String getLineComment() {
        return "#!";
    }

    @Override
    int getCurrentLineNumber() {
        return currentLineNum;
    }

    @Override
    List<String> getHeader() {
        return header;
    }

    @Override
    boolean passesFileNameCheck(final String inputFilePath) {
        try {
            final Path p = IOUtil.getPath(inputFilePath);

            return p.getFileName().toString().toLowerCase().endsWith("." + GTF_FILE_EXTENSION);
        }
        catch (final FileNotFoundException ex) {
            logger.warn("File does not exist! - " + inputFilePath + " - returning name check as failure.");
        }
        catch (final IOException ex) {
            logger.warn("Caught IOException on file: " + inputFilePath + " - returning name check as failure.");
        }

        return false;
    }

    @Override
    List<String> readActualHeader(final LineIterator reader) {

        // Read in the header lines:
        ingestHeaderLines(reader);

        // Validate our header:
        validateHeader(header, true);

        // Set our line number to be the line of the first actual Feature:
        currentLineNum = HEADER_NUM_LINES + 1;

        return header;
    }

    /**
     * Get the version information from the header.
     */
    private String getVersionFromHeader() {
        // header version is of the form:
        //     #!genome-version ASM584v2
        // So we get the stuff after the space:
        return header.get(1).split("[ \t]")[1];
    }

    @Override
    boolean validateFeatureSubtype(final GencodeGtfFeature feature) {
        return validateEnsemblGtfFeature( feature );
    }

    @Override
    void incrementLineNumber() {
        ++currentLineNum;
    }

    @Override
    String getUcscVersionNumber() {
        return getVersionFromHeader();
    }

    //==================================================================================================================
    // Static Methods:

    /**
     * Validates a given {@link GencodeGtfFeature} against a given version of the ENSEMBL GTF file spec.
     * This method ensures that all required fields are defined, but does not interrogate their values.
     * @param feature A {@link GencodeGtfFeature} to validate.  MUST NOT BE {@code null}.
     * @return True if {@code feature} contains all required fields for the given GENCODE GTF version, {@code gtfVersion}
     */
    static boolean validateEnsemblGtfFeature(final GencodeGtfFeature feature) {

        final GencodeGtfFeature.FeatureType featureType = feature.getFeatureType();

        if ( featureType != GencodeGtfFeature.FeatureType.GENE) {
            if (feature.getTranscriptId() == null) {
                return false;
            }
            if (feature.getTranscriptType() == null) {
                return false;
            }
            if (feature.getTranscriptName() == null) {
                return false;
            }
        }

        return true;
    }

    //==================================================================================================================
    // Instance Methods:

    /**
     * Check if the given header of a tentative ENSEMBL GTF file is, in fact, the header to such a file.
     * @param header Header lines to check for conformity to ENSEMBL GTF specifications.
     * @param throwIfInvalid If true, will throw a {@link UserException.MalformedFile} if the header is invalid.
     * @return true if the given {@code header} is that of a ENSEMBL GTF file; false otherwise.
     */
    @VisibleForTesting
    boolean validateHeader(final List<String> header, final boolean throwIfInvalid) {
        if ( header.size() != HEADER_NUM_LINES) {
            if ( throwIfInvalid ) {
                throw new UserException.MalformedFile(
                        "ENSEMBL GTF Header is of unexpected length: " +
                                header.size() + " != " + HEADER_NUM_LINES);
            }
            else {
                return false;
            }
        }

        // Check the normal commented fields:
        return checkHeaderLineStartsWith(header, 0, "genome-build") &&
               checkHeaderLineStartsWith(header, 1, "genome-version") &&
               checkHeaderLineStartsWith(header, 2, "genome-date") &&
               checkHeaderLineStartsWith(header, 3, "genome-build-accession") &&
               checkHeaderLineStartsWith(header, 4, "genebuild-last-updated");
    }

    //==================================================================================================================
    // Helper Data Types:

}
