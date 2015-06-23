package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.datasources.VariantShard;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;

/**
 * KeyReadByVariantShard takes a PCollection of reads and keys each read by the variants shards it overlaps. Note that
 * if a read span several shards multiple KVs will be produced. This part of the larger transform to join reads with
 * overlapping variants.
 *
 * |---- shard 0 -----|---- shard 1 -----|---- shard 2 -----|---- shard 3 -----|---- shard 4 -----|
 *           |-------------- read 1 --------------|
 *  results in
 *  KV<shard 0, read 1>
 *  KV<shard 1, read 1>
 *  KV<shard 2, read 1>
 */
public class KeyReadByVariantShard extends PTransform<PCollection<GATKRead>, PCollection<KV<VariantShard, GATKRead>>> {
    private static final long serialVersionUID = 1L;

    @Override
    public PCollection<KV<VariantShard, GATKRead>> apply( PCollection<GATKRead> input ) {
        return input.apply(ParDo.of(new DoFn<GATKRead, KV<VariantShard, GATKRead>>() {
            private static final long serialVersionUID = 1L;
            @Override
            public void processElement(ProcessContext c) throws Exception {
                List<VariantShard> shards = VariantShard.getVariantShardsFromInterval(c.element());
                for (VariantShard shard : shards) {
                    c.output(KV.of(shard, c.element()));
                }
            }
        }).named("KeyReadByVariantShard"));
    }
}
