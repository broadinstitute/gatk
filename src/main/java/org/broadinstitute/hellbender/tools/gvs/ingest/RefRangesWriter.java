package org.broadinstitute.hellbender.tools.gvs.ingest;

import java.io.Closeable;
import java.io.IOException;

public abstract class RefRangesWriter implements Closeable {
    public abstract void write(long location, long sampleId, int length, String state) throws IOException;

    public void commitData() {
        // no-op
    }
}
