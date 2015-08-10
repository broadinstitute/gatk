package org.broadinstitute.hellbender.engine.dataflow.datasources;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;

public interface ReferenceSource {
    public ReferenceBases getReferenceBases( final SimpleInterval interval ) throws IOException;
}
