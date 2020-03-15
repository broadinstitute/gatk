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
import java.util.*;

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

    private static String            VERSION_FIELD    = "genome-version";
    private static String            DEFAULT_VERSION  = "ENSEMBL_DEFAULT_VERSION";
    private static final Set<String> COMMENT_PREFIXES = Collections.unmodifiableSet(new LinkedHashSet<>(Arrays.asList("#!", "##")));

    //==================================================================================================================
    // Private Members:

    private final        List<String> header          = new ArrayList<>();
    private              int          currentLineNum  = 1;
    private              String       version         = null;

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
    String getDefaultLineComment() {
        return COMMENT_PREFIXES.iterator().next();
    }

    @Override
    Set<String> getAllLineComments() {
        return COMMENT_PREFIXES;
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
        currentLineNum = header.size() + 1;

        // Set up our version number:
        populateVersionNumber();

        return header;
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
        return version;
    }

    @Override
    /**
     * {@inheritDoc}
     *
     * Because ENSEMBL GTF files are strictly a superset of GENCODE GTF files, we need to do some extra checks here to
     * make sure that this file can NOT be decoded by {@link GencodeGtfCodec} but can still be decoded by this
     * {@link EnsemblGtfCodec}.
     */
    public boolean canDecode(final String inputFilePath) {

        // Create a GencodeGtfCodec so we can see if it will decode the input file.
        final GencodeGtfCodec gencodeGtfCodec = new GencodeGtfCodec();
        if ( gencodeGtfCodec.canDecode(inputFilePath) ) {
            // Uh oh!  We can decode this as GENCODE.
            // So we should NOT decode this as ENSEMBL.
            return false;
        }

        return super.canDecode(inputFilePath);
    }

    //==================================================================================================================
    // Static Methods:

    /**
     * Validates a given {@link GencodeGtfFeature} against a given version of the ENSEMBL GTF file spec.
     * This method ensures that all required fields are defined, but does not interrogate their values.
     * @param feature A {@link GencodeGtfFeature} to validate.  MUST NOT BE {@code null}.
     * @return True if {@code feature} contains all required fields for the given GENCODE GTF version, {@code gtfVersion}
     */
    private static boolean validateEnsemblGtfFeature(final GencodeGtfFeature feature) {

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

    private void populateVersionNumber() {

        // If `genome-version` was specified in the header, we should use that.
        // Otherwise we can return a placeholder.

        String ver = DEFAULT_VERSION;

        // Attempt to get the version from the header:
        for ( final String line : header ) {
            for ( final String comment : getAllLineComments() ) {
                if ( line.startsWith(comment + VERSION_FIELD) ) {
                    ver = line.replaceFirst(comment + VERSION_FIELD + "\\s*", "").trim();
                }
            }
        }

        version = ver;
    }

    /**
     * Check if the given header of a tentative ENSEMBL GTF file is, in fact, the header to such a file.
     * @param header Header lines to check for conformity to ENSEMBL GTF specifications.
     * @param throwIfInvalid If true, will throw a {@link UserException.MalformedFile} if the header is invalid.
     * @return true if the given {@code header} is that of a ENSEMBL GTF file; false otherwise.
     */
    @VisibleForTesting
    boolean validateHeader(final List<String> header, final boolean throwIfInvalid) {
        // As it turns out, the ENSEMBL GTF header is pretty loosy-goosy.
        // No fields are required, and therefore it could actually be empty.

        // Rather than attempting to validate the file, here we just
        // assert that all header lines begin with a comment (they should already).
        int lineNum = 1;
        for (final String line : header) {

            if ( !isLineCommented(line) ) {
                if ( throwIfInvalid ) {
                    throw new UserException.MalformedFile("ENSEMBL GTF Header line " + lineNum + " is not commented: " + line);
                }
                else {
                    return false;
                }
            }

            ++lineNum;
        }

        return true;
    }

    //==================================================================================================================
    // Helper Data Types:

}
