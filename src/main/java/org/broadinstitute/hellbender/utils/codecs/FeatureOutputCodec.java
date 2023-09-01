package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.engine.GATKPath;

import java.io.IOException;
import java.util.List;

/**
 * A FeatureOutputCodec can encode Features into some type of FeatureSink.
 *
 * It probably makes sense to have most concrete implementations on this interface implement
 * FeatureCodec as well so that the encoding and decoding happen in the same class.  But this
 * doesn't extend FeatureCodec because it's not strictly necessary.
 *
 * (Note that a FeatureCodec isn't a codec at all: It's just a dec.  It doesn't know how to encode.)
 */
public interface FeatureOutputCodec<F extends Feature, S extends FeatureSink<F>> {
    /** File formats are typically recognized by their extensions.  Is this our kind of file? */
    boolean canDecode( String path );

    /** Each codec knows how to handle one subtype of Feature, and this is the one. */
    Class<F> getFeatureType();

    /** Create a consumer of features. */
    S makeSink( GATKPath path, SAMSequenceDictionary dict, List<String> sampleNames, int compressionLevel );

    /** Push a feature into the consumer. */
    void encode( F feature, S sink ) throws IOException;

    /**
     * Get an object that wraps the usual FeatureSink that post-processes the interval-ordered
     * stream of features.
     * Typically, this handles multiple features on the same interval:
     * It might impose additional ordering criteria, uniqueness criteria, or record-merging behavior
     * according to the requirements of each feature type.
     */
    FeatureSink<F> makeSortMerger( GATKPath path, SAMSequenceDictionary dict, List<String> sampleNames, int compressionLevel );
}
