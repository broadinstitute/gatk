package org.broadinstitute.hellbender.engine.dataflow.datasources;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;
import java.io.Serializable;

public class ReferenceDataflowSource implements Serializable {
    private static final long serialVersionUID = 1L;

    public ReferenceDataflowSource(String referenceName, PipelineOptions pipelineOptions) {
        // if reference name ends with .fasta and hdfs - use Hadoop, if not hdfs - use local, otherwise - use GCS
    }

    public ReferenceBases getReferenceBases( final SimpleInterval interval ) throws IOException {
        return null;
    }
}
