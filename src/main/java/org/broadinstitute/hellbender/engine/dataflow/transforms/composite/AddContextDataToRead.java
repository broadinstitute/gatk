package org.broadinstitute.hellbender.engine.dataflow.transforms.composite;

import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadContextData;
import org.broadinstitute.hellbender.utils.read.Read;

import java.util.List;

public class AddContextDataToRead extends PTransform<PCollection<Read>, PCollection<KV<Read, ReadContextData>>> {

    private final String referenceName;
    private final List<String> variantSources;

    public AddContextDataToRead( final String referenceName, final List<String> variantSources ) {
        this.referenceName = referenceName;
        this.variantSources = variantSources;
    }

    @Override
    public PCollection<KV<Read, ReadContextData>> apply( PCollection<Read> input ) {
        return super.apply(input);
    }

}
