package org.broadinstitute.hellbender.engine.dataflow;

import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.Read;

import java.util.List;

public class KeyLocatableByVariantShardDataflowTransform extends PTransform<PCollection<Locatable>, PCollection<KV<VariantShard, Locatable>>> {

    @Override
    public PCollection<KV<VariantShard, Locatable>> apply( PCollection<Locatable> input ) {
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
