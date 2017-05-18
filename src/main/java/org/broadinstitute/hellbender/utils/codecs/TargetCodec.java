package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetTableReader;

import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;

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

    @Override
    public Object readActualHeader(final LineIterator reader) {
        sourceReader = new LineIteratorReader(reader);
        try {
            sourceTargetTableReader = new TargetTableReader(sourceReader);
        } catch (final IOException ex) {
            throw new GATKException("cannot read target table", ex);
        }
        return sourceTargetTableReader.columns();
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
