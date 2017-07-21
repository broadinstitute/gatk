package org.broadinstitute.hellbender.utils.codecs.GENCODE;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.LocationAware;
import htsjdk.tribble.AbstractFeatureCodec;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.readers.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * {@link htsjdk.tribble.Tribble} Codec to read data from a GENCODE GTF file.
 *
 * GENCODE GTF Files are defined here: https://www.gencodegenes.org/data_format.html
 *
 * Created by jonn on 7/21/17.
 */
final public class GencodeGtfCodec extends AbstractFeatureCodec<GencodeGtfFeature, LineIterator> {

    protected final Logger logger = LogManager.getLogger(this.getClass());

    private static final String COMMENT_START = "##";

    static final int NUM_COLUMNS = 9;

    private long currentLineNum = 1;
    protected List<String> header = new ArrayList<>();

    // ============================================================================================================

    public GencodeGtfCodec() {
        super(GencodeGtfFeature.class);
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
    public FeatureCodecHeader readHeader(final LineIterator lineIterator) throws IOException {
        return new FeatureCodecHeader(readActualHeader(lineIterator), FeatureCodecHeader.NO_HEADER_END);
    }

    @Override
    public LocationAware makeIndexableSourceFromStream(final InputStream bufferedInputStream) {
        final PositionalBufferedStream pbs;
        if (bufferedInputStream instanceof PositionalBufferedStream) {
            pbs = (PositionalBufferedStream) bufferedInputStream;
        } else {
            pbs = new PositionalBufferedStream(bufferedInputStream);
        }
        return new AsciiLineReaderIterator(new AsciiLineReader(pbs));
    }

    // ============================================================================================================

    @Override
    public GencodeGtfFeature decode(final LineIterator lineIterator) {

        GencodeGtfFeature decodedFeature = null;

        // Create some caches for our data (as we need to group it):
        GencodeGtfGeneFeature gene = null;
        GencodeGtfTranscriptFeature transcript = null;
        final ArrayList< GencodeGtfExonFeature > exonStore = new ArrayList<>();
        final ArrayList< GencodeGtfFeature > leafFeatureStore = new ArrayList<>();

        // Accumulate lines until we have a full gene and all of its internal features:
        while ( lineIterator.hasNext() ) {

            String line = lineIterator.peek();

            // We must assume we can get header lines.
            // If we get a header line, we return null.
            // This allows indexing to work.
            if ( line.startsWith(COMMENT_START) ) {
                lineIterator.next();
                return null;
            }

            String[] splitLine = line.split("\t", -1);

            // Ensure the file is at least trivially well-formed:
            if (splitLine.length != NUM_COLUMNS) {
                throw new UserException.MalformedFile("Found an invalid number of columns in the given GENCODE file on line "
                        + currentLineNum + " - Given: " + splitLine.length + " Expected: " + NUM_COLUMNS + " : " + line);
            }

            // We need to key off the feature type to collapse our accumulated records:
            final String featureType = splitLine[2];

            // Create a baseline feature to add into our data:
            GencodeGtfFeature feature = GencodeGtfFeature.create(splitLine);

            // Make sure we keep track of the line number for if and when we need to write the file back out:
            feature.setFeatureOrderNumber(currentLineNum);

            // Once we see another gene, we take all accumulated records and combine them into a
            // GencodeGtfFeature.
            if ((gene != null) && (transcript != null) && (featureType.equals("gene") || featureType.equals("transcript") )) {

                aggregateRecordsIntoGeneFeature(gene, transcript, exonStore, leafFeatureStore);

                if ( featureType.equals("gene") ) {
                    decodedFeature = gene;
                    break;
                }
                else if ( featureType.equals("transcript") ) {
                    transcript = (GencodeGtfTranscriptFeature) feature;
                    ++currentLineNum;
                }
            }
            else {
                // We have not reached the end of this set of gene / transcript records.
                // We must cache these records together so we can create a meaningful data hierarchy from them all.
                // Records are stored in their Feature form, not string form.

                // Add the feature to the correct storage unit for easy assembly later:
                switch (featureType) {
                    case "gene":
                        gene = (GencodeGtfGeneFeature)feature;
                        break;
                    case "transcript":
                        transcript = (GencodeGtfTranscriptFeature)feature;
                        break;
                    case "exon":
                        exonStore.add((GencodeGtfExonFeature)feature);
                        break;
                    default:
                        leafFeatureStore.add(feature);
                        break;
                }

                ++currentLineNum;
            }

            // Increment our iterator here so we don't accidentally miss any features from the following gene
            lineIterator.next();
        }

        // Do we have some records leftover that we need to aggregate into a feature:
        if ( (gene != null) &&
                (exonStore.size() != 0) || (leafFeatureStore.size() != 0) ) {

            aggregateRecordsIntoGeneFeature(gene, transcript, exonStore, leafFeatureStore);
            decodedFeature = gene;
        }

        // If we have other records left over we should probably yell a lot,
        // as this is bad.
        //
        // However, this should never actually happen.
        //
        if ( (exonStore.size() != 0) || (leafFeatureStore.size() != 0) ) {

            if (exonStore.size() != 0) {
                logger.error("Gene Feature Aggregation: Exon store not empty: " + exonStore.toString());
            }

            if (leafFeatureStore.size() != 0) {
                logger.error("Gene Feature Aggregation: leaf feature store not empty: " + leafFeatureStore.toString());
            }

            String msg = "Aggregated data left over after parsing complete: Exons: " + exonStore.size() + " ; LeafFeatures: " + leafFeatureStore.size();
            throw new RuntimeException(msg);
        }

        return decodedFeature;
    }

    /**
     * Read the {@code header} from the given {@link LineIterator} for the GENCODE GTF File.
     * Will also validate this {@code header} for correctness before returning it.
     * Throws a {@link UserException.MalformedFile} if the header is malformed.
     * @param reader The {@link LineIterator} from which to read the header.
     * @return The header as read from the {@code reader}
     */
    private List<String> readActualHeader(LineIterator reader) {

        boolean isFirst = true;

        while ( reader.hasNext() ) {
            String line = reader.peek();

            // The file will start with commented out lines.
            // Grab them until there are no more commented out lines.
            if ( line.startsWith(COMMENT_START) ) {
                header.add(line);
                reader.next();
                isFirst = false;
            }
            else if ( isFirst ) {
                throw new UserException.MalformedFile("GENCODE file does not have a header!");
            }
            else {
                // Validate our header:
                if ( !validateHeader(header.toArray(new String[0])) ) {
                    throw new UserException.MalformedFile("Invalid GENCODE GTF File - Header does not conform to GENCODE GTF Specifications!");
                }

                break;
            }
        }

        // Set our line number to be the line of the first actual Feature:
        currentLineNum = 6;

        return header;
    }

    @Override
    public boolean canDecode(String path) {

        // Simple file and name checks to start with:
        File f = new File(path);
        boolean canDecode = f.getName().toLowerCase().startsWith("gencode") && f.getName().toLowerCase().endsWith(".gtf");

        if ( canDecode ) {

            String localPath = path;
            if ( path.startsWith("file://") ) {
                localPath = path.substring(7);
            }

            // Crack open the file and look at the top of it:
            try {
                try (BufferedReader br = new BufferedReader(new FileReader(localPath))) {
                    // Read the first 5 lines.
                    // They compose the header of a valid GTF File.
                    String[] stringArray = new String[5];

                    for (int i = 0; i < stringArray.length; ++i){
                        stringArray[i] = br.readLine();
                    }

                    // Validate our header:
                    canDecode = validateHeader(stringArray);
                }
            }
            catch (FileNotFoundException ex) {
                logger.warn("File does not exist! - " + path + " - returning can decode as failure.");
                canDecode = false;
            }
            catch (IOException ex) {
                logger.warn("Caught IOException on file: " + path + " - returning can decode as failure.");
                canDecode = false;
            }
        }
        else {
            logger.warn("Given file name does not conform to GENCODE GTF standards: " + path);
        }

        return canDecode;
    }

    // ============================================================================================================

    /**
     * Check if the given header of a tentative GENCODE GTF file is, in fact, the header to such a file.
     * @param header Header lines to check for conformity to GENCODE GTF specifications.
     * @return true if the given {@code header} is that of a GENCODE GTF file; false otherwise.
     */
    static boolean validateHeader(final String[] header) {
        boolean isValid = false;

        if ( header.length == 5 ) {
            // Check the normal commented fields:
            isValid = header[0].startsWith("##description:");
            isValid = isValid && header[1].startsWith("##provider: GENCODE");

            isValid = isValid && header[2].startsWith("##contact: gencode");
            isValid = isValid && header[2].endsWith("@sanger.ac.uk");

            isValid = isValid && header[3].startsWith("##format: gtf");
            isValid = isValid && header[4].startsWith("##date:");
        }

        return isValid;
    }

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
    private void aggregateRecordsIntoGeneFeature(final GencodeGtfGeneFeature gene,
                                                                  final GencodeGtfTranscriptFeature transcript,
                                                                  final ArrayList< GencodeGtfExonFeature > exonStore,
                                                                  final ArrayList< GencodeGtfFeature > leafFeatureStore ) {

        // OK, we go through the record and consolidate the sub parts of the record.
        // We must consolidate these records through grouping by genomic position.

        // Loop through the Exons and put the correct leaf features into each:
        for ( GencodeGtfExonFeature exon : exonStore ) {
            for (Iterator<GencodeGtfFeature> iterator = leafFeatureStore.iterator(); iterator.hasNext(); ) {

                GencodeGtfFeature feature = iterator.next();

                // Features that are within the extents of an exon belong in that exon:
                if ( exon.contains(feature) ) {

                    GencodeGtfFeature.FeatureType featureType = feature.getFeatureType();

                    // Add the feature to the correct place in the exon:
                    switch (featureType) {
                        case cds:
                            exon.setCds((GencodeGtfCDSFeature) feature);
                            break;
                        case start_codon:
                            exon.setStartCodon((GencodeGtfStartCodonFeature) feature);
                            break;
                        case stop_codon:
                            exon.setStopCodon((GencodeGtfStopCodonFeature) feature);
                            break;
                        case utr:
                            transcript.addUtr((GencodeGtfUTRFeature) feature);
                            break;
                        case selenocysteine:
                            transcript.addSelenocysteine(((GencodeGtfSelenocysteineFeature) feature));
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


}
