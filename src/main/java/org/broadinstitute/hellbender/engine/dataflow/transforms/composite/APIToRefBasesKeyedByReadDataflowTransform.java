package org.broadinstitute.hellbender.engine.dataflow.transforms.composite;

import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.transforms.GroupReadsForRef;
import org.broadinstitute.hellbender.engine.dataflow.transforms.RefBasesFromAPI;
import org.broadinstitute.hellbender.engine.dataflow.transforms.GroupReadWithRefBases;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceShard;
import org.broadinstitute.hellbender.utils.read.Read;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

/**
 * Created by davidada on 5/19/15.
 */
public class APIToRefBasesKeyedByReadDataflowTransform extends PTransform<PCollection<Read>, PCollection<KV<Read, ReferenceBases>>> {
    private final String refName;

    public APIToRefBasesKeyedByReadDataflowTransform(String refName) {

        this.refName = refName;
    }
    @Override
    public PCollection<KV<Read, ReferenceBases>> apply(PCollection<Read> input) {
        PCollection<KV<ReferenceShard, Read>> shardAndRead = input.apply(new GroupReadsForRef());
        PCollection<KV<ReferenceBases, Iterable<Read>>> GroupedReads = shardAndRead.apply(new RefBasesFromAPI(refName));
        return GroupedReads.apply(new GroupReadWithRefBases());
    }
}
