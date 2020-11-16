package org.broadinstitute.hellbender.utils.io;

import htsjdk.tribble.Feature;

import java.io.Closeable;

/**
 * Interface for output streams that encode Tribble {@link Feature}s
 */
public interface FeatureOutputStream<F extends Feature> extends Closeable {
    void writeHeader(final String header);
    void add(final F feature);
    void close();
}
