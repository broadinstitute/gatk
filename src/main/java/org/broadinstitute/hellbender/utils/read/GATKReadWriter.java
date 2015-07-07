package org.broadinstitute.hellbender.utils.read;

import java.io.Closeable;

/**
 * Interface for classes that are able to write GATKReads to some output destination.
 */
public interface GATKReadWriter extends Closeable {
    void addRead(final GATKRead read);
}
