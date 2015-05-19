package org.broadinstitute.hellbender.engine.dataflow;

import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.utils.read.Read;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

/**
 * Created by davidada on 5/19/15.
 */
public class APIToRefBasesKeyedByReadDataflowTransform extends PTransform<PCollection<Read>, PCollection<KV<Read, ReferenceBases>>> {
    @Override
    public PCollection<KV<Read, ReferenceBases>> apply(PCollection<Read> input) {
        PCollection<KV<ReferenceShard, Read>> shardAndRead = input.apply(new GroupReadsForRef());
        PCollection<KV<ReferenceBases, Iterable<Read>>> GroupedReads = shardAndRead.apply(new RefBasesFromAPIDataflowTransform());
        return GroupedReads.apply(new RefBasesToReadsDataflowTransform());
    }
}
