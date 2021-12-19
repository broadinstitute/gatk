package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.tribble.Feature;

public interface FeatureSink<F extends Feature> {
    void write( F feature );
    void close();
}
