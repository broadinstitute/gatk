package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.datasources.VariantShard;
import org.broadinstitute.hellbender.utils.variant.Variant;

import java.util.List;

/**
 * Created by davidada on 5/15/15.
 */
public class KeyVariantByShard extends PTransform<PCollection<Variant>, PCollection<KV<VariantShard, Variant>>> {
    @Override
    public PCollection<KV<VariantShard, Variant>> apply(PCollection<Variant> input) {
        return input.apply(ParDo.of(new DoFn<Variant, KV<VariantShard, Variant>>() {
            @Override
            public void processElement(ProcessContext c) throws Exception {
                List<VariantShard> shards = VariantShard.getVariantShardsFromInterval(c.element());
                for (VariantShard shard : shards) {
                    c.output(KV.of(shard, c.element()));
                }
            }
        }).named("KeyVariantByShard"));
    }

}
