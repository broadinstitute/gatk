package org.broadinstitute.hellbender.utils.codecs.gtf;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.AbstractFeatureCodec;
import htsjdk.tribble.readers.LineIterator;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Path;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * {@link htsjdk.tribble.Tribble} Codec to read data from a GENCODE GTF file.
 *
 * GENCODE GTF Files are defined here: https://www.gencodegenes.org/pages/data_format.html
 *
 * This codec will scan through a GENCODE GTF file and return {@link GencodeGtfFeature} objects.
 * {@link GencodeGtfFeature} objects contain fields that have sub-features.  All features are
 * grouped by gene (this is the natural formatting of a GENCODE GTF file).
 *
 * All fields exist in the Abstract {@link GencodeGtfFeature}.  The subclasses contain representations of the logical
 * data hierarchy that reflect how the data were presented in the feature file itself (to preserve the natural
 * grouping by gene).
 * The {@link GencodeGtfFeature} logical data hierarchy (NOT the class hierarchy) is as follows
 * (with | representing a "has a" relationship)
 *
 * +--> {@link GencodeGtfGeneFeature}
 *    |
 *    +--> {@link GencodeGtfTranscriptFeature}
 *       |
 *       +--> {@link GencodeGtfSelenocysteineFeature}
 *       +--> {@link GencodeGtfUTRFeature}
 *       +--> {@link GencodeGtfExonFeature}
 *          |
 *          +--> {@link GencodeGtfCDSFeature}
 *          +--> {@link GencodeGtfStartCodonFeature}
 *          +--> {@link GencodeGtfStopCodonFeature}
 *
 * {@link htsjdk.tribble.Tribble} indexing has been tested and works as expected.
 * Does not support {@link htsjdk.tribble.index.tabix.TabixIndex} indexing.
 *
 * Unlike many other {@link htsjdk.tribble.Tribble} codecs, this one scans multiple input file lines to produce
 * a single feature.  This is due to how GENCODE GTF files are structured (essentially grouped by contig and gene).
 * For this reason, {@link GencodeGtfCodec} inherits from {@link AbstractFeatureCodec}, as opposed to {@link htsjdk.tribble.AsciiFeatureCodec}
 * (i.e. {@link htsjdk.tribble.AsciiFeatureCodec}s read a single line at a time, and {@link AbstractFeatureCodec} do not have that explicit purpose).
 *
 * Created by jonn on 7/21/17.
 */
final public class GencodeGtfCodec extends AbstractGtfCodec {

    private static final Logger logger = LogManager.getLogger(GencodeGtfCodec.class);

    // ============================================================================================================

    public static final String GENCODE_GTF_FILE_PREFIX = "gencode";
    public static final String GTF_FILE_TYPE_STRING = "GENCODE";

    // ============================================================================================================

    private static final int GENCODE_GTF_MIN_VERSION_NUM_INCLUSIVE = 19;

    /**
     * Maximum version of gencode that will not generate a warning.  This parser will still attempt to parse versions above this number, but a warning about potential errors will appear.
     */
    private static final int GENCODE_GTF_MAX_VERSION_NUM_INCLUSIVE = 34;

    private int currentLineNum = 1;
    private final List<String> header = new ArrayList<>();
    private static final int HEADER_NUM_LINES = 5;

    private static final Pattern VERSION_PATTERN = Pattern.compile("version (\\d+)");
    private int versionNumber;

    private static final String commentPrefix = "##";

    // ============================================================================================================

    /**
     * Gets the UCSC version corresponding to the given gencode version.
     * Version equivalences obtained here:
     *
     *  https://genome.ucsc.edu/FAQ/FAQreleases.html
     *  https://www.gencodegenes.org/human/releases.html
     *
     * @param gencodeVersion The gencode version to convert to UCSC version.
     * @return The UCSC version in a {@link String} corresponding to the given gencode version.
     */
    private static String getUcscVersionFromGencodeVersion(final int gencodeVersion) {
        if (gencodeVersion < GENCODE_GTF_MIN_VERSION_NUM_INCLUSIVE) {
            throw new GATKException("Gencode version is too far out of date.  Cannot decode: " + gencodeVersion);
        }

        if ( gencodeVersion < 25 ) {
            return "hg19";
        }
        else {
            return "hg38";
        }
    }

    // ============================================================================================================

    public GencodeGtfCodec() {
        super();
    }

    // ============================================================================================================

    @Override
    int getCurrentLineNumber() {
        return currentLineNum;
    }

    @Override
    List<String> getHeader() {
        return header;
    }

    @Override
    List<String> readActualHeader(final LineIterator reader) {

        // Clear our version number too:
        versionNumber = -1;

        // Read in the header lines:
        ingestHeaderLines(reader);

        // Validate our header:
        validateHeader(header, true);

        // Set our version number:
        setVersionNumber();

        // Set our line number to be the line of the first actual Feature:
        currentLineNum = header.size() + 1;

        return header;
    }

    /**
     * Sets {@link #versionNumber} to the number corresponding to the value in the header.
     */
    private void setVersionNumber() {
        try {
            final Matcher versionMatcher = VERSION_PATTERN.matcher(header.get(0));
            if (!versionMatcher.find() ) {
                throw new UserException.MalformedFile("Cannot find version number from Gencode GTF header.");
            }
            versionNumber = Integer.valueOf(versionMatcher.group(1));
        }
        catch (final NumberFormatException ex) {
            throw new UserException("Could not read version number from header", ex);
        }
    }

    /**
     * Validates a given {@link GencodeGtfFeature} against a given version of the GENCODE GTF file spec.
     * This method ensures that all required fields are defined, but does not interrogate their values.
     * @param feature A {@link GencodeGtfFeature} to validate.  MUST NOT BE {@code null}.
     * @param gtfVersion The GENCODE GTF version against which to validate {@code feature}
     * @return True if {@code feature} contains all required fields for the given GENCODE GTF version, {@code gtfVersion}
     */
    private static boolean validateGencodeGtfFeature(final GencodeGtfFeature feature, final int gtfVersion) {

        if (gtfVersion < GencodeGtfCodec.GENCODE_GTF_MIN_VERSION_NUM_INCLUSIVE) {
            throw new GATKException("Invalid version number for validation: " + gtfVersion +
                    " must be above: " + GencodeGtfCodec.GENCODE_GTF_MIN_VERSION_NUM_INCLUSIVE);
        }

        final Level logLevel = Level.FATAL;

        if (feature.getGeneName() == null) {
            logger.log(logLevel, "Feature gene name is null.");
            return false;
        }

        final GencodeGtfFeature.FeatureType featureType = feature.getFeatureType();

        if ( gtfVersion < 26 ) {
            if (feature.getGeneStatus() == null) {
                logger.log(logLevel, "Gencode version < 26 and feature gene status is null.");
                return false;
            }
            if (feature.getTranscriptStatus() == null) {
                logger.log(logLevel, "Gencode version < 26 and feature transcript status is null.");
                return false;
            }
        }

        if ( (featureType != GencodeGtfFeature.FeatureType.GENE) ||
                (gtfVersion < 21) ) {
            if (feature.getTranscriptId() == null) {
                logger.log(logLevel, "Gencode version < 21 and feature transcript ID is null.");
                return false;
            }
            if (feature.getTranscriptType() == null) {
                logger.log(logLevel, "Gencode version < 21 and feature transcript type is null.");
                return false;
            }
            if (feature.getTranscriptName() == null) {
                logger.log(logLevel, "Gencode version < 21 and feature transcript name is null.");
                return false;
            }
        }

        // Gencode can only have 2 feature types:
        if (!feature.getAnnotationSource().equals(GencodeGtfFeature.ANNOTATION_SOURCE_ENSEMBL) &&
                !feature.getAnnotationSource().equals(GencodeGtfFeature.ANNOTATION_SOURCE_HAVANA) ) {
            logger.log(logLevel, "Gencode data came from an unexpected source: " + feature.getAnnotationSource());
            return false;
        }

        return true;
    }

    @Override
    boolean passesFileNameCheck(final String inputFilePath) {
        try {
            final Path p = IOUtil.getPath(inputFilePath);

            return p.getFileName().toString().toLowerCase().startsWith(GENCODE_GTF_FILE_PREFIX) &&
                    p.getFileName().toString().toLowerCase().endsWith("." + GTF_FILE_EXTENSION);
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
    String getDefaultLineComment() {
        return commentPrefix;
    }

    @Override
    Set<String> getAllLineComments() {
        return Collections.unmodifiableSet(new HashSet<>(Collections.singletonList(commentPrefix)));
    }

    @Override
    String  getGtfFileType() {
        return GTF_FILE_TYPE_STRING;
    }

    @Override
    boolean validateFeatureSubtype(final GencodeGtfFeature feature ) {
        return validateGencodeGtfFeature( feature, versionNumber );
    }

    @Override
    void incrementLineNumber() {
        ++currentLineNum;
    }

    @Override
    String getUcscVersionNumber() {
        return getUcscVersionFromGencodeVersion(versionNumber);
    }

    // ============================================================================================================

    /**
     * Check if the given header of a tentative GENCODE GTF file is, in fact, the header to such a file.
     * Will also return true if the file is a general GTF file (i.e. a GTF file that was not created and
     * maintained by GENCODE).
     * @param header Header lines to check for conformity to GENCODE GTF specifications.
     * @param throwIfInvalid If true, will throw a {@link UserException.MalformedFile} if the header is invalid.
     * @return true if the given {@code header} is that of a GENCODE GTF file; false otherwise.
     */
    @VisibleForTesting
    boolean validateHeader(final List<String> header, final boolean throwIfInvalid) {
        if ( header.size() != HEADER_NUM_LINES) {
            if ( throwIfInvalid ) {
                throw new UserException.MalformedFile(
                        "GENCODE GTF Header is of unexpected length: " +
                                header.size() + " != " + HEADER_NUM_LINES);
            }
            else {
                return false;
            }
        }

        // Check the normal commented fields:
        if ( !checkHeaderLineStartsWith(header,0, "description:") ) {
            return false;
        }

        if ( !header.get(0).contains("version") ) {
            if ( throwIfInvalid ) {
                throw new UserException.MalformedFile(
                        "GENCODE GTF Header line 1 does not contain version specification: " +
                                header.get(0));
            }
            else {
                return false;
            }
        }

        // Grab the version from the file and make sure it's within the acceptable range:
        final Matcher versionMatcher = VERSION_PATTERN.matcher(header.get(0));
        if ( !versionMatcher.find() ) {
            if ( throwIfInvalid ) {
                throw new UserException.MalformedFile(
                        "GENCODE GTF Header line 1 does not contain a recognizable version number: " +
                                header.get(0));
            }
            else {
                return false;
            }
        }

        try {
            final int versionNumber = Integer.valueOf(versionMatcher.group(1));
            if (versionNumber < GENCODE_GTF_MIN_VERSION_NUM_INCLUSIVE) {
                final String message = "GENCODE GTF Header line 1 has an out-of-date (< v" + GENCODE_GTF_MIN_VERSION_NUM_INCLUSIVE + " version number (" +
                        versionNumber + "): " + header.get(0);
                if (throwIfInvalid) {
                    throw new UserException.MalformedFile(message);
                } else {
                    logger.warn(message + "   Continuing, but errors may occur.");
                }
            }

            if (versionNumber > GENCODE_GTF_MAX_VERSION_NUM_INCLUSIVE) {
                logger.warn("GENCODE GTF Header line 1 has a version number that is above maximum tested version (v " + GENCODE_GTF_MAX_VERSION_NUM_INCLUSIVE + ") (given: " +
                        versionNumber + "): " + header.get(0) + "   Continuing, but errors may occur.");
            }
        }
        catch (final NumberFormatException ex) {
            if ( throwIfInvalid ) {
                throw new UserException("Could not create number value for version: " + versionMatcher.group(1), ex);
            }
            else {
                return false;
            }
        }

        return checkHeaderLineStartsWith(header, 1, "provider: GENCODE") &&
                checkHeaderLineStartsWith(header, 2, "contact: gencode") &&
                checkHeaderLineStartsWith(header, 3, "format: gtf") &&
                checkHeaderLineStartsWith(header, 4, "date:");
    }
}
