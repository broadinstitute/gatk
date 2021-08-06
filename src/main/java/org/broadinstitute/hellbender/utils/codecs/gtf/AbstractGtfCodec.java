package org.broadinstitute.hellbender.utils.codecs.gtf;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.LocationAware;
import htsjdk.tribble.AbstractFeatureCodec;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.readers.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

public abstract class AbstractGtfCodec extends AbstractFeatureCodec<GencodeGtfFeature, LineIterator> {

    static final Logger logger = LogManager.getLogger(AbstractGtfCodec.class);

    //==================================================================================================================
    // Public Static Members:
    public static final String GTF_FILE_EXTENSION = "gtf";

    //==================================================================================================================
    // Private/Protected Static Members:
    static final String FIELD_DELIMITER = "\t";

    static final int NUM_COLUMNS = 9;
    static final int FEATURE_TYPE_FIELD_INDEX = 2;

    // A number much larger than the number of header lines we expect (to future proof):
    static final int MAX_NUM_HEADER_LINES_TO_CHECK = 100;

    //==================================================================================================================
    // Private Members:


    //==================================================================================================================
    // Constructors:
    protected AbstractGtfCodec() {
        super(GencodeGtfFeature.class);
    }

    //==================================================================================================================
    // Override Methods:

    @Override
    public boolean canDecode(final String inputFilePath) {
        boolean canDecode;
        try {
            // Simple file and name checks to start with:
            final Path p = IOUtil.getPath(inputFilePath);

            canDecode = passesFileNameCheck(inputFilePath);

            if (canDecode) {

                // Crack open the file and look at the top of it:
                try ( final BufferedReader br = new BufferedReader(new InputStreamReader(Files.newInputStream(p))) ) {

                    // The first few lines compose the header of a valid GTF File:
                    final List<String> headerLines = new ArrayList<>(MAX_NUM_HEADER_LINES_TO_CHECK);

                    // Only check a few lines at the top:
                    for (int i = 0; i < MAX_NUM_HEADER_LINES_TO_CHECK; ++i) {
                        final String line = br.readLine();
                        if ( line == null ) {
                            break;
                        }

                        // When lines are no longer commented we're done with the header:
                        if ( isLineCommented(line) ) {
                            headerLines.add(line);
                        }
                        else {
                            break;
                        }
                    }

                    // Validate our header:
                    canDecode = validateHeader(headerLines);
                }

            }
        }
        catch (final FileNotFoundException ex) {
            logger.warn("File does not exist! - " + inputFilePath + " - returning can decode as failure.");
            canDecode = false;
        }
        catch (final IOException ex) {
            logger.warn("Caught IOException on file: " + inputFilePath + " - returning can decode as failure.");
            canDecode = false;
        }

        return canDecode;
    }

    @Override
    public GencodeGtfFeature decode(final LineIterator lineIterator) {

        GencodeGtfFeature decodedFeature = null;

        // Create some caches for our data (as we need to group it):
        GencodeGtfGeneFeature gene = null;
        GencodeGtfTranscriptFeature transcript = null;
        final List<GencodeGtfExonFeature> exonStore = new ArrayList<>();
        final List<GencodeGtfFeature> leafFeatureStore = new ArrayList<>();

        boolean needToFlushRecords = false;

        // Accumulate lines until we have a full gene and all of its internal features:
        while ( lineIterator.hasNext() ) {

            final String line = lineIterator.peek();

            // We must assume we can get header lines.
            // If we get a header line, we return null.
            // This allows indexing to work.
            if ( isLineCommented(line) ) {
                lineIterator.next();
                return null;
            }

            // Split the line into different GTF Fields
            final String[] splitLine = splitGtfLine(line);

            // We need to key off the feature type to collapse our accumulated records:
            final GencodeGtfFeature.FeatureType featureType = GencodeGtfFeature.FeatureType.getEnum( splitLine[FEATURE_TYPE_FIELD_INDEX] );

            // Create a baseline feature to add into our data:
            final GencodeGtfFeature feature = GencodeGtfFeature.create(splitLine, getGtfFileType());

            // Make sure we keep track of the line number for if and when we need to write the file back out:
            feature.setFeatureOrderNumber(getCurrentLineNumber());

            // Set our UCSC version number:
            feature.setUcscGenomeVersion(getUcscVersionNumber());

            // Once we see another gene we take all accumulated records and combine them into the
            // current GencodeGtfFeature.
            // Then we then break out of the loop and return the last full gene object.
            if ((gene != null) && (featureType == GencodeGtfFeature.FeatureType.GENE)) {

                aggregateRecordsIntoGeneFeature(gene, transcript, exonStore, leafFeatureStore);

                // If we found a new gene line, we set our decodedFeature to be
                // the gene we just finished building.
                //
                // We intentionally break here so that we do not call lineIterator.next().
                // This is so that the new gene (i.e. the one that triggered us to be in this if statement)
                // remains intact for the next call to decode.
                decodedFeature = gene;

                needToFlushRecords = false;

                break;
            }
            // Once we see a transcript we aggregate our data into our current gene object and
            // set the current transcript object to the new transcript we just read.
            // Then we continue reading from the line iterator.
            else if ((transcript != null) && (featureType == GencodeGtfFeature.FeatureType.TRANSCRIPT)) {

                aggregateRecordsIntoGeneFeature(gene, transcript, exonStore, leafFeatureStore);

                transcript = (GencodeGtfTranscriptFeature) feature;
                incrementLineNumber();

                needToFlushRecords = true;
            }
            else {
                // We have not reached the end of this set of gene / transcript records.
                // We must cache these records together so we can create a meaningful data hierarchy from them all.
                // Records are stored in their Feature form, not string form.

                // Add the feature to the correct storage unit for easy assembly later:
                switch (featureType) {
                    case GENE:
                        gene = (GencodeGtfGeneFeature)feature;
                        break;
                    case TRANSCRIPT:
                        transcript = (GencodeGtfTranscriptFeature)feature;
                        break;
                    case EXON:
                        exonStore.add((GencodeGtfExonFeature)feature);
                        break;
                    default:
                        leafFeatureStore.add(feature);
                        break;
                }

                needToFlushRecords = false;
                incrementLineNumber();
            }

            // Increment our iterator here so we don't accidentally miss any features from the following gene
            lineIterator.next();
        }

        // For the last record in the file, we need to do one final check to make sure that we don't miss it.
        // This is because there will not be a subsequent `gene` line to read:
        if ( (gene != null) && (needToFlushRecords || (!exonStore.isEmpty()) || (!leafFeatureStore.isEmpty())) ) {

            aggregateRecordsIntoGeneFeature(gene, transcript, exonStore, leafFeatureStore);
            decodedFeature = gene;
        }

        // If we have other records left over we should probably yell a lot,
        // as this is bad.
        //
        // However, this should never actually happen.
        //
        if ( (!exonStore.isEmpty()) || (!leafFeatureStore.isEmpty()) ) {

            if (!exonStore.isEmpty()) {
                logger.error("Gene Feature Aggregation: Exon store not empty: " + exonStore.toString());
            }

            if (!leafFeatureStore.isEmpty()) {
                logger.error("Gene Feature Aggregation: leaf feature store not empty: " + leafFeatureStore.toString());
            }

            final String msg = "Aggregated data left over after parsing complete: Exons: " + exonStore.size() + " ; LeafFeatures: " + leafFeatureStore.size();
            throw new GATKException.ShouldNeverReachHereException(msg);
        }

        // Now we validate our feature before returning it:
        if ( ! validateFeature(decodedFeature) ) {
            throw new UserException.MalformedFile("Decoded feature is not valid: " + decodedFeature);
        }

        return decodedFeature;
    }

    // ============================================================================================================
    // Trivial override methods that are pulled form AsciiFeatureCodec
    // This was done to ensure that this was a reasonable Codec class (with good interfaces for reading features).

    @Override
    public void close(final LineIterator lineIterator) {
        CloserUtil.close(lineIterator);
    }

    @Override
    public boolean isDone(final LineIterator lineIterator) {
        return !lineIterator.hasNext();
    }

    @Override
    public LineIterator makeSourceFromStream(final InputStream bufferedInputStream) {
        return new LineIteratorImpl(new SynchronousLineReader(bufferedInputStream));
    }

    @Override
    public FeatureCodecHeader readHeader(final LineIterator lineIterator) {
        return new FeatureCodecHeader(readActualHeader(lineIterator), FeatureCodecHeader.NO_HEADER_END);
    }

    @Override
    public LocationAware makeIndexableSourceFromStream(final InputStream bufferedInputStream) {
        return new AsciiLineReaderIterator(AsciiLineReader.from(bufferedInputStream));
    }

    //==================================================================================================================
    // Static Methods:

    /**
     * Aggregates the given feature sets into a single gene feature.
     *
     * The given gene is updated using modifiers.
     * {@code exonStore} and {@code leafFeatureStore} are cleared of all data.
     *
     * @param gene {@link GencodeGtfGeneFeature} into which to aggregate features.
     * @param transcript {@link GencodeGtfTranscriptFeature} to insert into {@code gene}
     * @param exonStore {@link List} of {@link GencodeGtfExonFeature}s to insert into corresponding {@link GencodeGtfTranscriptFeature} {@code transcript}
     * @param leafFeatureStore {@link List} of {@link GencodeGtfFeature}s to insert into corresponding {@link GencodeGtfExonFeature} objects in {@code exonStore}
     */
    private static void aggregateRecordsIntoGeneFeature(final GencodeGtfGeneFeature gene,
                                                final GencodeGtfTranscriptFeature transcript,
                                                final List< GencodeGtfExonFeature > exonStore,
                                                final List< GencodeGtfFeature > leafFeatureStore ) {

        // OK, we go through the record and consolidate the sub parts of the record.
        // We must consolidate these records through grouping by genomic position.

        // Loop through the Exons and put the correct leaf features into each:
        for ( final GencodeGtfExonFeature exon : exonStore ) {
            for ( final Iterator<GencodeGtfFeature> iterator = leafFeatureStore.iterator(); iterator.hasNext(); ) {

                final GencodeGtfFeature feature = iterator.next();

                // Features that are within the extents of an exon belong in that exon:
                if ( exon.contains(feature) ) {

                    final GencodeGtfFeature.FeatureType featureType = feature.getFeatureType();

                    // Add the feature to the correct place in the exon:
                    switch (featureType) {
                        case CDS:
                            exon.setCds((GencodeGtfCDSFeature) feature);
                            break;
                        case START_CODON:
                            exon.setStartCodon((GencodeGtfStartCodonFeature) feature);
                            break;
                        case STOP_CODON:
                            exon.setStopCodon((GencodeGtfStopCodonFeature) feature);
                            break;
                        case UTR:
                            transcript.addUtr((GencodeGtfUTRFeature) feature);
                            break;
                        case SELENOCYSTEINE:
                            transcript.addSelenocysteine((GencodeGtfSelenocysteineFeature) feature);
                            break;
                        case FIVE_PRIME_UTR:
                            transcript.addFivePrimeUtr((GencodeGtfFivePrimeUtrFeature) feature);
                            break;
                        case THREE_PRIME_UTR:
                            transcript.addThreePrimeUtr((GencodeGtfThreePrimeUtrFeature) feature);
                            break;
                        default:
                            throw new UserException.MalformedFile(
                                    "Found unexpected Feature Type in GENCODE GTF File (line " +
                                            feature.getFeatureOrderNumber() + "): " +
                                            featureType.toString()
                            );
                    }

                    // We have used this iterator item.
                    // We should remove it now so we don't keep going through the list each exon.
                    iterator.remove();
                }
            }

            // Now insert this exon into the transcript:
            transcript.addExon(exon);
        }

        // Add in the transcript:
        gene.addTranscript(transcript);

        // Clear the input data:
        exonStore.clear();
        leafFeatureStore.clear();
    }

    //==================================================================================================================
    // Instance Methods:

    /**
     * Check if the given header of a tentative GTF file is, in fact, the header to such a file.
     * @param header Header lines to check for conformity to GTF specifications.
     * @return true if the given {@code header} is that of a GTF file; false otherwise.
     */
    @VisibleForTesting
    boolean validateHeader(final List<String> header) {
        return validateHeader(header, false);
    }

    /**
     * Validate all fields in the given {@link GencodeGtfFeature}.
     * @param feature {@link GencodeGtfFeature} to validate.
     * @return {@code true} IFF the given {@code feature} / {@code versionNumber} combination is valid.  {@code false} otherwise.
     */
    private boolean validateFeature(final GencodeGtfFeature feature) {
        return validateBaseGtfFeatureFields(feature) && validateFeatureSubtype(feature);
    }

    /**
     * Validates a given {@link GencodeGtfFeature} for basic GTF criteria.
     * This method ensures that required fields are defined, but does not interrogate their values.
     * @param feature A {@link GencodeGtfFeature} to validate.
     * @return True if {@code feature} contains all required base GTF fields.
     */
    private static boolean validateBaseGtfFeatureFields(final GencodeGtfFeature feature) {

        if ( feature == null ) {
            return false;
        }

        final GencodeGtfFeature.FeatureType featureType = feature.getFeatureType();

        if (feature.getChromosomeName() == null) {
            return false;
        }
        if (feature.getAnnotationSource() == null) {
            return false;
        }
        if (feature.getFeatureType() == null) {
            return false;
        }
        if (feature.getGenomicStrand() == null) {
            return false;
        }
        if (feature.getGenomicPhase() == null) {
            return false;
        }

        if (feature.getGeneId() == null) {
            return false;
        }
        if (feature.getGeneType() == null) {
            return false;
        }
        if (feature.getGeneName() == null) {
            return false;
        }

        if ( (featureType != GencodeGtfFeature.FeatureType.GENE) &&
             (featureType != GencodeGtfFeature.FeatureType.TRANSCRIPT) &&
             (featureType != GencodeGtfFeature.FeatureType.SELENOCYSTEINE) &&
                (featureType != GencodeGtfFeature.FeatureType.FIVE_PRIME_UTR) &&
                (featureType != GencodeGtfFeature.FeatureType.THREE_PRIME_UTR) ) {

            if (feature.getExonNumber() == GencodeGtfFeature.NO_EXON_NUMBER) {
                return false;
            }
            return feature.getExonId() != null;
        }

        return true;
    }

    /**
     * Split the given line in a GTF file into fields.
     * Throws a {@link UserException} if the file is not valid.
     * @param line {@link String} containing one line of a GTF file to split.
     * @return A {@link String[]} with each entry containing a field from the GTF line.
     */
    private String[] splitGtfLine(final String line) {
        // Split the line into different GTF Fields
        // Note that we're using -1 as the limit so that empty tokens will still be counted
        // (as opposed to discarded).
        final String[] splitLine = line.split(FIELD_DELIMITER, -1);

        // Ensure the file is at least trivially well-formed:
        if (splitLine.length != NUM_COLUMNS) {
            throw new UserException.MalformedFile("Found an invalid number of columns in the given GTF file on line "
                    + getCurrentLineNumber() + " - Given: " + splitLine.length + " Expected: " + NUM_COLUMNS + " : " + line);
        }
        return splitLine;
    }

    /**
     * Read in lines from the given {@link LineIterator} and put them in the header file.
     * Will read until the lines no longer start with comments.
     * @param reader {@link LineIterator} a reader pointing at the top of a GTF file.
     */
    void ingestHeaderLines(final LineIterator reader) {
        while ( reader.hasNext() ) {
            final String line = reader.peek();

            // The file will start with commented out lines.
            // Grab them until there are no more commented out lines.
            if ( isLineCommented(line) ) {

                getHeader().add(line);
                reader.next();
            }
            else {
                break;
            }
        }
    }

    /**
     * Checks that the given header line number starts with the given text.
     * Uses {@link #getAllLineComments()} to check for comments at the line starts.
     * @param header A {@link List<String>} containing a header to validate.
     * @param lineNum Line number in the header to check.
     * @param startingText {@link String} containing text that the line should start with
     * @return {@code true} IFF the header line number {@code lineNum} starts with {@code startingText}; {@code false} otherwise.
     */
    boolean checkHeaderLineStartsWith(final List<String> header, final int lineNum, final String startingText) {
        return checkHeaderLineStartsWith(header, lineNum, startingText, false);
    }

    /**
     * Checks that the given header line number starts with the given text.
     * Uses {@link #getAllLineComments()} to check for comments at the line starts.
     * @param header A {@link List<String>} containing a header to validate.
     * @param lineNum Line number in the header to check.
     * @param startingText {@link String} containing text that the line should start with
     * @param throwIfInvalid If {@code true} will throw a {@link UserException} instead of returning false.
     * @return {@code true} IFF the header line number {@code lineNum} starts with {@code startingText}; {@code false} otherwise.
     */
    boolean checkHeaderLineStartsWith(final List<String> header, final int lineNum, final String startingText, final boolean throwIfInvalid ) {
        final boolean foundGoodLine = isLineCommented(header.get(lineNum), startingText);

        if (!foundGoodLine) {
            if ( throwIfInvalid ) {
                throw new UserException.MalformedFile(
                        getGtfFileType() + " GTF Header line " + (lineNum + 1) + " does not contain expected information (" +
                                getDefaultLineComment() + startingText + "): " + header.get(lineNum));
            }
            else {
                return false;
            }
        }

        return true;
    }

    /**
     * Checks whether the given line is commented out or not.
     * @param line A {@link String} representing a line in the GTF file.
     * @return {@code true} iff the line starts with a comment delimiter.  {@code false} otherwise.
     */
    boolean isLineCommented(final String line) {
        return isLineCommented(line, "");
    }

    /**
     * Checks whether the given line is commented out or not and if the content of the line begins with the given
     * {@code linePrefix}.
     * @param line A {@link String} representing a line in the GTF file.
     * @param linePrefix A {@link String} containing content that the line must start with in addition to the comment prefix.
     * @return {@code true} iff the line starts with a comment delimiter.  {@code false} otherwise.
     */
    boolean isLineCommented(final String line, final String linePrefix) {
        boolean isCommented = false;
        for ( final String commentPrefix : getAllLineComments() ) {
            if ( line.startsWith(commentPrefix + linePrefix) ) {
                isCommented = true;
                break;
            }
        }

        return isCommented;
    }

    /** @return The current line number for this {@link AbstractGtfCodec}. */
    abstract int getCurrentLineNumber();

    /**
     * Increments the current line number for this {@link AbstractGtfCodec} by one.
     *
     * This method is needed because of how Tribble initializes codecs.
     * Each codec needs its own internal state for what line number it is on.  This parent class cannot store the state,
     * or the child codecs will not work.
     */
    abstract void incrementLineNumber();

    /** @return The header AbstractGtfCodec. */
    abstract List<String> getHeader();

    /**
     * @return The default {@link String} a line begins with to indicate that line is commented out.
     */
    abstract String getDefaultLineComment();

    /**
     * @return A {@link Set<String>} containing all possible line prefixes that indicate the line is commented out.
     */
    abstract Set<String> getAllLineComments();

    /**
     * @return The type of GTF file in this {@link AbstractGtfCodec}.
     */
    abstract String getGtfFileType();

    /**
     * @param inputFilePath A {@link String} containing the path to a potential GTF file.
     * @return {@code true} IFF the given {@code inputFilePath} is a valid name for this {@link AbstractGtfCodec}.
     */
    abstract boolean passesFileNameCheck(final String inputFilePath);

    /**
     * Check if the given header of a tentative GTF file is, in fact, the header to such a file.
     * @param header Header lines to check for conformity to GTF specifications.
     * @param throwIfInvalid If true, will throw a {@link UserException.MalformedFile} if the header is invalid.
     * @return true if the given {@code header} is that of a GTF file; false otherwise.
     */
    @VisibleForTesting
    abstract boolean validateHeader(final List<String> header, final boolean throwIfInvalid);

    /**
     * Validate the given {@link GencodeGtfFeature} according to what type of feature it is.
     * @param feature {@link GencodeGtfFeature} to validate.
     * @return {@code true} IFF the given {@code feature} / {@code versionNumber} combination is valid.  {@code false} otherwise.
     */
    abstract boolean validateFeatureSubtype(final GencodeGtfFeature feature);

    /**
     * Get the UCSC Version Number of this {@link AbstractGtfCodec}.
     * @return The {@link String} representation of the UCSC version number corresponding to the backing data source behind this {@link AbstractGtfCodec}.
     */
    abstract String getUcscVersionNumber();

    /**
     * Read the {@code header} from the given {@link LineIterator} for the GTF File.
     * Will also validate this {@code header} for correctness before returning it.
     * Throws a {@link UserException.MalformedFile} if the header is malformed.
     *
     * This must be called before {@link #decode(LineIterator)}
     *
     * @param reader The {@link LineIterator} from which to read the header.
     * @return The header as read from the {@code reader}
     */
    abstract List<String> readActualHeader(final LineIterator reader);

    //==================================================================================================================
    // Helper Data Types:

}
