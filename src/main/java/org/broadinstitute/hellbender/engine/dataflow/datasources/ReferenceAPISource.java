package org.broadinstitute.hellbender.engine.dataflow.datasources;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;
import java.io.Serializable;

public class ReferenceAPISource implements ReferenceSource, Serializable {
    private static final long serialVersionUID = 1L;

    private RefAPIMetadata refAPIMetadata;

    public ReferenceAPISource(RefAPIMetadata refAPIMetadata) {
        this.refAPIMetadata = refAPIMetadata;
    }

    @Override
    public ReferenceBases getReferenceBases(SimpleInterval interval, PipelineOptions pipelineOptions) throws IOException {
        RefAPISource refAPISource = RefAPISource.getRefAPISource();
        return refAPISource.getReferenceBases(pipelineOptions, refAPIMetadata, interval);
    }
}
