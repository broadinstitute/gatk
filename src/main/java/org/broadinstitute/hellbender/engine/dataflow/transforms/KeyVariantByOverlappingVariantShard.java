package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.dev.DoFnWLog;
import org.broadinstitute.hellbender.engine.dataflow.datasources.VariantShard;
import org.broadinstitute.hellbender.utils.variant.Variant;

import java.util.List;

/**
 * KeyVariantByOverlappingVariantShard takes a PCollection of variants and keys each variant by the variants shards
 * it overlaps. Note that if a variant span several shards multiple KVs will be produced. This is part of the larger
 * transform to join reads with overlapping variants.
 *
 * |---- shard 0 -----|---- shard 1 -----|---- shard 2 -----|---- shard 3 -----|---- shard 4 -----|
 *           |-------------- variant 1 --------------|
 *  results in
 *  KV<shard 0, variant 1>
 *  KV<shard 1, variant 1>
 *  KV<shard 2, variant 1>
 */
public class KeyVariantByOverlappingVariantShard extends PTransform<PCollection<Variant>, PCollection<KV<VariantShard, Variant>>> {
    private static final long serialVersionUID = 1L;

    @Override
    public PCollection<KV<VariantShard, Variant>> apply( PCollection<Variant> input ) {
        return input.apply(ParDo.of(new DoFnWLog<Variant, KV<VariantShard, Variant>>("KeyVariantByOverlappingVariantShard") {
            private static final long serialVersionUID = 1L;
            @Override
            public void processElement( ProcessContext c ) throws Exception {
                List<VariantShard> shards = VariantShard.getVariantShardsFromInterval(c.element());
                for (VariantShard shard : shards) {
                    c.output(KV.of(shard, c.element()));
                }
            }
        }).named("KeyVariantByVariantShard"));
    }
}
