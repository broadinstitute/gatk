package org.broadinstitute.hellbender.engine.dataflow.transforms.composite;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.transforms.*;
import com.google.cloud.dataflow.sdk.transforms.join.CoGbkResult;
import com.google.cloud.dataflow.sdk.transforms.join.CoGroupByKey;
import com.google.cloud.dataflow.sdk.transforms.join.KeyedPCollectionTuple;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.TupleTag;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadContextData;
import org.broadinstitute.hellbender.engine.dataflow.datasources.VariantShard;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.Read;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
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
        // Then, get the KV<Read, ReferenceBases>
        PCollection<KV<Read, Iterable<Variant>>> keyedVariants = input.apply(new LoadAndKeyVariants(variantSources, pipeline));
        PCollection<KV<Read, ReferenceBases>> keyedReference = input.apply(new APIToRefBasesKeyedByReadDataflowTransform(referenceName));
        // Now group by Read and put into the ___ class.

        // GroupBy Read
        // NOTE: Assuming that both have been deduped.
        final TupleTag<Iterable<Variant>> variantTag = new TupleTag<>();
        final TupleTag<ReferenceBases> referenceTag = new TupleTag<>();
        PCollection<KV<Read, CoGbkResult>> coGbkInput = KeyedPCollectionTuple
                .of(variantTag, keyedVariants)
                .and(referenceTag, keyedReference).apply(CoGroupByKey.<Read>create());

        // GroupBy Read
        return coGbkInput.apply(ParDo.of(
                new DoFn<KV<Read, CoGbkResult>, KV<Read, ReadContextData>>() {
                    @Override
                    public void processElement(ProcessContext c) throws Exception {
                        Iterable<Variant> v = c.element().getValue().getOnly(variantTag);
                        ReferenceBases refs = c.element().getValue().getOnly(referenceTag);
                        c.output(KV.of(c.element().getKey(), new ReadContextData(refs, v)));
                    }
                }));
    }

}
