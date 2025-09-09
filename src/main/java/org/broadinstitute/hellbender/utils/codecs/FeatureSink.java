package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.tribble.Feature;

public interface FeatureSink<F extends Feature> extends AutoCloseable {
    void write( F feature );
    @Override
    void close();
}
