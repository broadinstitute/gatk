package org.broadinstitute.hellbender.tools.gvs.ingest;

import java.io.Closeable;
import java.io.IOException;

public interface RefRangesWriter extends Closeable {
    void write(long location, long sampleId, int length, String state) throws IOException;
    void writeCompressed(long packedData, long sampleId) throws IOException;

    default void commitData() {
        // no-op
    }
}
