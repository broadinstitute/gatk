package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.LocationAware;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.readers.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetTableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.URI;
import java.util.Collections;
import java.util.Iterator;

import static htsjdk.tribble.bed.BEDCodec.BED_EXTENSION;

/**
 * Target table file codec.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class TargetCodec extends AsciiFeatureCodec<Target> {

    private LineIteratorReader sourceReader;
    private TargetTableReader sourceTargetTableReader;

    @SuppressWarnings("unused") // used thru reflection.
    public TargetCodec() {
        super(Target.class);
    }

    @Override
    public Target decode(final String line) {
        return sourceTargetTableReader.readRecord(line);
    }

    /**
     * TargetCodec uses a com.opencsv.CSVReader, which uses a buffered reader. Since just initializing the
     * CSVReader causes it to buffer a large portion of the underlying input stream, it can't be used directly
     * on an input stream that is being indexed.
     *
     * To avoid this, we override makeIndexableSourceFromStream and instantiate and return a marker class that is a
     * subclass of AsciiLineReaderIterator, and which is recognized by the readActualHeader method (where the actual
     * TargetTableReader/CSVReader state is established). The marker class signals the codec to manually extract
     * everything up to the end of the header from the stream, and instantiate the CSVReader on a separate stream
     * containing only the (header) portion of the input. This allows the CSV reader to correctly establish it's
     * internal state based on the header, without consuming the rest of the stream.
     *
     * The rest of the indexing process is driven by the indexer, which pulls from the input directly and passes
     * candidate features to the codec's decode metho), so the features are then consumed and indexed with the correct
     * stream offsets.
     *
     * @param inputStream stream to be indexed
     * @return An AsciiLineReaderIterator implementation.
     */
    @Override
    public LocationAware makeIndexableSourceFromStream(final InputStream inputStream) {
        return new IndexableAsciiLineReaderIterator(AsciiLineReader.from(inputStream));
    }

    // Marker class to indicate to readActualHeader that we're being indexed
    private class IndexableAsciiLineReaderIterator extends AsciiLineReaderIterator {
        public IndexableAsciiLineReaderIterator(final AsciiLineReader asciiLineReader) {
            super(asciiLineReader);
        }
    }

    /**
     * Read the actual header. This method instantiates the underlying TargetTableReader on either the entirety
     * of the passed in LineIterator (in the case where we're not indexing) or on the subset of the input lines
     * representing the header lines (in the case that were indexing).
     * @param reader
     * @return column header object
     */
    @Override
    public Object readActualHeader(final LineIterator reader) {
        // If we're being indexed, create a shim Reader that consumes only the header portion of the input
        // stream, leaving the rest of the stream underlying the LineIterator available for the indexer
        sourceReader = reader instanceof IndexableAsciiLineReaderIterator ?
                makeIndexableReaderShim(reader) :
                new LineIteratorReader(reader);

        try {
            sourceTargetTableReader = new TargetTableReader(sourceReader);
        } catch (final IOException ex) {
            throw new GATKException("cannot read target table", ex);
        }
        return sourceTargetTableReader.columns();
    }

    // Return a LineIteratorReader backed by only the header line from the underlying reader.
    // This method needs to be certain not to consume any input from the stream beyond the headerline
    // (that is, it shoud consume no features).
    //
    // NOTE: This implementation discards leading comments, which is based on the assumption
    // that TableReader.processComment is a no-op.
    private LineIteratorReader makeIndexableReaderShim(final LineIterator reader) {
        // find the first non-comment line (which must be the header)
        String line = null;
        while (reader.hasNext()) {
            line = reader.peek();
            if (line != null && line.startsWith(TableUtils.COMMENT_PREFIX)) {
                reader.next();
            } else {
                break;
            }
        }

        // return a line reader iterator backed by the single string representing the header line
        final String sourceLine = line;
        return new LineIteratorReader(
                new LineIteratorImpl(
                    new LineReader() {
                        private Iterator<String> iterator = Collections.singletonList(sourceLine).iterator();

                        @Override
                        public String readLine() throws IOException {
                            return iterator.hasNext() ? iterator.next() : null;
                        }

                        @Override
                        public void close() {
                        }
                    }
                )
        );
}

    @Override
    public boolean canDecode(final String path) {
        File file;
        try {
            // Use the URI constructor so that we can handle file:// URIs
            final URI uri = new URI(path);
            file = uri.isAbsolute() ? new File(uri) : new File(path);
        }
        catch ( Exception e ) {  // Contract for canDecode() mandates that all exception be trapped
            return false;
        }
        
        if (!file.canRead() || !file.isFile()) {
            return false;
        }
        try {
            new TargetTableReader(file).close();
        } catch (final IOException|RuntimeException ex) {
            // IOException correspond to low-level IO errors
            // whereas RuntimeException would be caused by a formatting error in the file.
            return false;
        }

        //disallow .bed extension
        final String toDecode = AbstractFeatureReader.hasBlockCompressedExtension(path) ?
                path.substring(0, path.lastIndexOf(".")) :
                path;
        return !toDecode.toLowerCase().endsWith(BED_EXTENSION);
    }
}
