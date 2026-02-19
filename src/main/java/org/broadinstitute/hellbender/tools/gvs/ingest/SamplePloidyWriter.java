package org.broadinstitute.hellbender.tools.gvs.ingest;

import java.io.Closeable;
import java.io.IOException;

public abstract class SamplePloidyWriter implements Closeable {
    public abstract void write(Long sampleId, Long chromosome, Integer ploidy) throws IOException;

    public void commitData() {
        // no-op
    }
}
