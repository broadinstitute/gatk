package org.broadinstitute.hellbender.utils.io;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import org.broadinstitute.hellbender.engine.GATKPath;

import java.util.function.Function;

/**
 * Convenient for creating the appropriate type of {@link FeatureOutputStream} for a given output path
 */
public final class FeatureOutputStreamFactory {
    /**
     * Creates block compressed output stream with index if possible, otherwise plain text
     */
    public <F extends Feature> FeatureOutputStream<F> create(final GATKPath path,
                                                             final FeatureCodec<? extends Feature, ?> codec,
                                                             final Function<F, String> encoder,
                                                             final SAMSequenceDictionary dictionary,
                                                             final int compressionLevel) {
        if (IOUtil.hasBlockCompressedExtension(path.toPath())) {
            return new TabixIndexedFeatureOutputStream<>(path, codec, encoder, dictionary, compressionLevel);
        }
        return new UncompressedFeatureOutputStream<>(path, encoder);
    }
}
