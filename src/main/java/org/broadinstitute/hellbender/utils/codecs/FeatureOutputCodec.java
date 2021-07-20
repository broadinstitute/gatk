package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.engine.GATKPath;

import java.io.IOException;
import java.util.List;

public interface FeatureOutputCodec<F extends Feature, S extends FeatureSink<F>> {
    boolean canDecode( final String path );
    Class<F> getFeatureType();
    S makeSink( GATKPath path, SAMSequenceDictionary dict, List<String> sampleNames, int compressionLevel );
    void encode( F feature, S sink  ) throws IOException;
}
