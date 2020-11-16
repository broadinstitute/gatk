package org.broadinstitute.hellbender.utils.io;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.PositionalOutputStream;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexCreator;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.IOException;
import java.nio.file.Path;
import java.util.function.Function;

/**
 * Encodes Tribble {@link Feature}s to a block compressed file with Tabix index, given a feature encoding function
 */
public final class TabixIndexedFeatureOutputStream<F extends Feature> implements FeatureOutputStream<F> {

    private static final String NEWLINE_CHARACTER = "\n";

    private final PositionalOutputStream outputStream;
    private final IndexCreator indexCreator;
    private final Function<F, String> encoder;
    private final Path featurePath;

    /**
     * @param file file to write to
     * @param codec feature codec
     * @param encoder encoding function mapping a feature to its text record (whithout a newline character)
     * @param dictionary sequence dictionary
     * @param compressionLevel block compression level
     */
    public TabixIndexedFeatureOutputStream(final GATKPath file, final FeatureCodec<F, Object> codec,
                                           final Function<F, String> encoder, final SAMSequenceDictionary dictionary,
                                           final int compressionLevel) {
        Utils.nonNull(file);
        Utils.nonNull(codec);
        Utils.nonNull(encoder);
        Utils.nonNull(dictionary);
        Utils.validateArg(IOUtil.hasBlockCompressedExtension(file.toPath()), "File must must have a gzip extension");
        featurePath = file.toPath();
        outputStream = new PositionalOutputStream(new BlockCompressedOutputStream(file.toString(), compressionLevel));
        indexCreator = new TabixIndexCreator(dictionary, codec.getTabixFormat());
        this.encoder = encoder;
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
    public void add(final F feature) {
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
            final Index index = indexCreator.finalizeIndex(outputStream.getPosition());
            index.writeBasedOnFeaturePath(featurePath);
            outputStream.close();
        } catch (final IOException e) {
            throw new GATKException("Error closing output", e);
        }
    }
}
