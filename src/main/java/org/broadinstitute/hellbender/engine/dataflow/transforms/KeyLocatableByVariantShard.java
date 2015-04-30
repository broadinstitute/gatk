package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.engine.dataflow.datasources.VariantShard;

import java.util.List;

public class KeyLocatableByVariantShard extends PTransform<PCollection<? extends Locatable>, PCollection<KV<VariantShard, Locatable>>> {

    @Override
    public PCollection<KV<VariantShard, Locatable>> apply( PCollection<? extends Locatable> input ) {
        return input.apply(ParDo.of(new DoFn<Locatable, KV<VariantShard, Locatable>>() {
            @Override
            public void processElement( ProcessContext c ) throws Exception {
                List<VariantShard> shards = VariantShard.getVariantShardsFromInterval(c.element());
                for (VariantShard shard : shards) {
                    c.output(KV.of(shard, c.element()));
                }
            }
        }).named("KeyLocatableByVariantShard"));
    }
}
