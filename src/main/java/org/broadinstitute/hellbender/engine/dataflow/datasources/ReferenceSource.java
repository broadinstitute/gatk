package org.broadinstitute.hellbender.engine.dataflow.datasources;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;

public interface ReferenceSource {
    public ReferenceBases getReferenceBases(SimpleInterval interval, PipelineOptions pipelineOptions) throws IOException;
}
