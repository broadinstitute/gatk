package org.broadinstitute.hellbender.utils.io;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.PositionalOutputStream;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexCreator;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.FeatureSink;
import org.broadinstitute.hellbender.utils.codecs.FeaturesHeader;

import java.io.IOException;
import java.nio.file.Path;
import java.util.function.Function;

/**
 * Class for output streams that encode Tribble {@link Feature}s. Supports block-compressed output, which is
 * detected by the output file path extension, in which case a Tabix index is also generated.
 */
public class FeatureOutputStream <F extends Feature> implements FeatureSink<F> {

    private static final String NEWLINE_CHARACTER = "\n";

    private final Function<F, String> encoder;
    private final PositionalOutputStream outputStream;
    private final IndexCreator indexCreator;
    private final Path featurePath;

    /**
     * @param file file to write to
     * @param tabixFormat column descriptions for the tabix index
     * @param header metaData for the features
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
            outputStream = new PositionalOutputStream(
                            new BlockCompressedOutputStream(file.toString(), compressionLevel));
            indexCreator = new TabixIndexCreator(dict, tabixFormat);
        } else {
            outputStream = new PositionalOutputStream(file.getOutputStream());
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
            indexCreator.addFeature(feature, outputStream.getPosition());
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
                final Index index = indexCreator.finalizeIndex(outputStream.getPosition());
                index.writeBasedOnFeaturePath(featurePath);
            }
            outputStream.close();
        } catch (final IOException e) {
            throw new GATKException("Error closing output", e);
        }
    }
}
