package org.broadinstitute.hellbender.engine.dataflow.transforms.composite;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadContextData;
import org.broadinstitute.hellbender.utils.read.Read;
import org.broadinstitute.hellbender.utils.variant.Variant;

import java.util.List;

public class AddContextDataToRead extends PTransform<PCollection<Read>, PCollection<KV<Read, ReadContextData>>> {

    private final String referenceName;
    private final List<String> variantSources;
    private final Pipeline pipeline;

    public AddContextDataToRead( final String referenceName, final List<String> variantSources, final Pipeline pipeline ) {
        this.referenceName = referenceName;
        this.variantSources = variantSources;
        this.pipeline = pipeline;
    }

    @Override
    public PCollection<KV<Read, ReadContextData>> apply( PCollection<Read> input ) {
        // First, get the KV<Read, Iterable<Variant>>
        PCollection<KV<Read, Iterable<Variant>>> keyedVariants = input.apply(new LoadAndKeyVariants(variantSources, pipeline));


        return super.apply(input);
    }

}
