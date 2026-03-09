package org.broadinstitute.hellbender.tools.gvs.ingest;

import java.io.Closeable;
import java.io.IOException;

public interface SamplePloidyWriter extends Closeable {
    void write(Long sampleId, Long chromosome, Integer ploidy) throws IOException;

    default void commitData() {
        // no-op
    }
}
