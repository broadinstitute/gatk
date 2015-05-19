package org.broadinstitute.hellbender.engine.dataflow.transforms.composite;

import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.ReadContextData;
import org.broadinstitute.hellbender.utils.read.Read;

public class AddContextDataToReads extends PTransform<PCollection<Read>, PCollection<KV<Read, ReadContextData>>> {

    @Override
    public PCollection<KV<Read, ReadContextData>> apply( PCollection<Read> input ) {
        return super.apply(input);
    }

}
