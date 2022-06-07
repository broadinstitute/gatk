package org.broadinstitute.hellbender.utils.io;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.LocationAware;
import htsjdk.samtools.util.PositionalOutputStream;
import htsjdk.tribble.Feature;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexCreator;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.FeatureSink;

import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Path;
import java.util.function.Function;

/**
 * Class for output streams that encode Tribble {@link Feature}s. Supports block-compressed output, which is
 * detected by the output file path extension, in which case a Tabix index is also generated.
 */
public class FeatureOutputStream <F extends Feature> implements FeatureSink<F> {

    private static final String NEWLINE_CHARACTER = "\n";

    private final Function<F, String> encoder;
    private final OutputStream outputStream;
    private final LocationAware locationAware; // same object as outputStream, above
    // This is necessary because there is no LocationAwareOutputStream class for
    // PositionalOutputStream and BlockCompressedOutputStream to extend.  Sadly.
    // Do not wrap the BlockCompressedOutputStream in a PositionalOutputStream. The
    // PositionalOutputStream just counts bytes output, but that's not a valid virtual file offset
    // for a bgzip-compressed file.

    private final IndexCreator indexCreator;
    private final Path featurePath;

    /**
     * @param file file to write to
     * @param tabixFormat column descriptions for the tabix index
     * @param encoder functor to transform feature into a line of text
     */
    public FeatureOutputStream( final GATKPath file,
                                final TabixFormat tabixFormat,
                                final Function<F, String> encoder,
                                final SAMSequenceDictionary dict,
                                final int compressionLevel ) {
        Utils.nonNull(file);
        Utils.nonNull(tabixFormat);
        Utils.nonNull(encoder);
        Utils.nonNull(dict);
        this.encoder = encoder;
        if (IOUtil.hasBlockCompressedExtension(file.toPath())) {
            final BlockCompressedOutputStream bcos =
                    new BlockCompressedOutputStream(file.toString(), compressionLevel);
            outputStream = bcos;
            locationAware = bcos;
            indexCreator = new TabixIndexCreator(dict, tabixFormat);
        } else {
            final PositionalOutputStream pos = new PositionalOutputStream(file.getOutputStream());
            outputStream = pos;
            locationAware = pos;
            indexCreator = null;
        }
        featurePath = file.toPath();
    }

    /**
     * Writes header to file. Should be called before adding features
     *
     * @param header header text (without final newline character), cannot be null
     */
    public void writeHeader(final String header) {
        Utils.nonNull(header);
        try {
            outputStream.write((header + NEWLINE_CHARACTER).getBytes());
        } catch (final IOException e) {
            throw new GATKException("Error writing header", e);
        }
    }

    /**
     * Adds new feature and writes to file
     *
     * @param feature feature to write, cannot be null
     */
    @Override
    public void write(final F feature) {
        Utils.nonNull(feature);
        if (indexCreator != null) {
            indexCreator.addFeature(feature, locationAware.getPosition());
        }
        try {
            outputStream.write((encoder.apply(feature) + NEWLINE_CHARACTER).getBytes());
        } catch (final IOException e) {
            throw new GATKException("Error writing record", e);
        }
    }

    /**
     * Closes the underlying stream and indexer
     */
    @Override
    public void close() {
        try {
            if (indexCreator != null) {
                final Index index = indexCreator.finalizeIndex(locationAware.getPosition());
                index.writeBasedOnFeaturePath(featurePath);
            }
            outputStream.close();
        } catch (final IOException e) {
            throw new GATKException("Error closing output", e);
        }
    }
}
