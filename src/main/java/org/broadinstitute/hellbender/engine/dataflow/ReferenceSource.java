package org.broadinstitute.hellbender.engine.dataflow;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

public class ReferenceSource {

    public ReferenceSource( final String referenceName, final PipelineOptions options) {

    }

    public static int getShardIDForInterval( final Locatable l ) {
        return 0;
    }


    public ReferenceBases getReferenceBases( final SimpleInterval interval ) {
        return null;
    }
}
