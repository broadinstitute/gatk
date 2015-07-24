package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.GroupByKey;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceShard;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * KeyReadByRefShard takes a PCollection of reads and keys each read by the reference shard using start position of
 * the read.
 * This is part of a transform that joins references bases to reads.
 *
 * |---- shard 0 -----|---- shard 1 -----|---- shard 2 -----|---- shard 3 -----|---- shard 4 -----|
 *           |-------------- read 1 --------------|
 *  results in
 *  KV<shard 0, read 1>
 */
public class KeyReadsByRefShard extends PTransform<PCollection<GATKRead>, PCollection<KV<ReferenceShard, Iterable<GATKRead>>>> {
    private static final long serialVersionUID = 1L;

    @Override
    public PCollection<KV<ReferenceShard, Iterable<GATKRead>>> apply(PCollection<GATKRead> input) {
        PCollection<KV<ReferenceShard, GATKRead>> keyReadByReferenceShard = input.apply(ParDo.of(new DoFn<GATKRead, KV<ReferenceShard, GATKRead>>() {
            private static final long serialVersionUID = 1L;
            @Override
            public void processElement(ProcessContext c) throws Exception {
                ReferenceShard shard = ReferenceShard.getShardNumberFromInterval(c.element());
                c.output(KV.of(shard, c.element()));
            }
        }).named("KeyReadByRefShard"));
        return keyReadByReferenceShard.apply(GroupByKey.<ReferenceShard, GATKRead>create());
    }
}

