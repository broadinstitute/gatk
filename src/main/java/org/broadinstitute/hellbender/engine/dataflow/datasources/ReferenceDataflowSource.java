package org.broadinstitute.hellbender.engine.dataflow.datasources;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.dataflow.sdk.transforms.SerializableFunction;
import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;
import java.io.Serializable;
import java.util.Map;

public class ReferenceDataflowSource implements ReferenceSource, Serializable {
    private static final long serialVersionUID = 1L;

    private ReferenceSource referenceSource;
    private SerializableFunction<GATKRead, SimpleInterval> referenceWindowFunction;

    public ReferenceDataflowSource(String referenceName, SerializableFunction<GATKRead, SimpleInterval> referenceWindowFunction,
                                   PipelineOptions pipelineOptions) {
        if (referenceName.endsWith(".fasta")) {
            if (BucketUtils.isHadoopUrl(referenceName)) {
                referenceSource = new ReferenceHadoopSource(referenceName);
            } else {
                referenceSource = new ReferenceFileSource(referenceName);
            }
        } else { // use GCS
            Map<String, String> referenceNameToIdTable = RefAPISource.buildReferenceNameToIdTable(pipelineOptions, referenceName);
            RefAPIMetadata refAPIMetadata = new RefAPIMetadata(referenceName, referenceNameToIdTable, referenceWindowFunction);
            referenceSource = new ReferenceAPISource(refAPIMetadata);
        }
        this.referenceWindowFunction = referenceWindowFunction;
    }

    @VisibleForTesting
    public ReferenceDataflowSource(RefAPIMetadata refAPIMetadata) {
        this.referenceSource = new ReferenceAPISource(refAPIMetadata);
        this.referenceWindowFunction = refAPIMetadata.getReferenceWindowFunction();
    }

    public SerializableFunction<GATKRead, SimpleInterval> getReferenceWindowFunction() {
        return referenceWindowFunction;
    }

    public ReferenceBases getReferenceBases(final SimpleInterval interval, PipelineOptions pipelineOptions) throws IOException {
        return referenceSource.getReferenceBases(interval, pipelineOptions);
    }
}
