package org.broadinstitute.hellbender.engine.dataflow.transforms;


import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.join.CoGbkResult;
import com.google.cloud.dataflow.sdk.transforms.join.CoGroupByKey;
import com.google.cloud.dataflow.sdk.transforms.join.KeyedPCollectionTuple;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionTuple;
import com.google.cloud.dataflow.sdk.values.TupleTag;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.Read;
import org.broadinstitute.hellbender.utils.variant.Variant;

public class GroupReadWithOverlappingData extends PTransform<KeyedPCollectionTuple<SimpleInterval>, PCollection<KV<Read, Iterable<Variant>>>> {

    private TupleTag<Read> readTag;
    private TupleTag<Variant> variantTag;

    @Override
    public PCollection<KV<Read, Iterable<Variant>>> apply( KeyedPCollectionTuple<SimpleInterval> input ) {
        PCollection<KV<SimpleInterval, CoGbkResult>> dataByShard = input.apply(CoGroupByKey.<SimpleInterval>create().withName("Group by shard"));

        dataByShard.apply(ParDo.of(new DoFn<KV<SimpleInterval, CoGbkResult>, KV<Read, Iterable<Variant>>>() {

            @Override
            public void processElement( ProcessContext c ) throws Exception {
                final SimpleInterval shard = c.element().getKey();
                final CoGbkResult shardData = c.element().getValue();

                final Iterable<Read> readsInShard = shardData.getAll(readTag);
                final Iterable<Variant> variantsInShard = shardData.getAll(variantTag);

                c.output(processShard(readsInShard, variantsInShard));
            }
        }));

        return null;
    }

    public void setTags( final Pair<TupleTag<Read>, TupleTag<Variant>> tags ) {
        this.readTag = tags.getLeft();
        this.variantTag = tags.getRight();
    }

    private KV<Read, Iterable<Variant>> processShard( final Iterable<Read> reads, final Iterable<Variant> variants ) {
        return null;
    }
}
