package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.transforms.*;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefWindowFunctions;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceShard;
import org.broadinstitute.hellbender.utils.SimpleInterval;
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
 *
 * If a custom reference window function is being used to map each read to arbitrary reference bases
 * (and not just the bases that overlap each read), that function should be passed in at construction.
 */
public class KeyReadsByRefShard extends PTransform<PCollection<GATKRead>, PCollection<KV<ReferenceShard, Iterable<GATKRead>>>> {
    private static final long serialVersionUID = 1L;

    private final SerializableFunction<GATKRead, SimpleInterval> referenceWindowFunction;

    public KeyReadsByRefShard() {
        this(RefWindowFunctions.IDENTITY_FUNCTION);
    }

    /**
     * @param referenceWindowFunction custom reference window function used to map each read to arbitrary reference bases
     */
    public KeyReadsByRefShard( final SerializableFunction<GATKRead, SimpleInterval> referenceWindowFunction ) {
        this.referenceWindowFunction = referenceWindowFunction;
    }

    @Override
    public PCollection<KV<ReferenceShard, Iterable<GATKRead>>> apply(PCollection<GATKRead> input) {
        PCollection<KV<ReferenceShard, GATKRead>> keyReadByReferenceShard = input.apply(ParDo.of(new DoFn<GATKRead, KV<ReferenceShard, GATKRead>>() {
            private static final long serialVersionUID = 1L;
            @Override
            public void processElement(ProcessContext c) throws Exception {
                // Apply our reference window function to each read before deciding which reference shard it belongs to.
                ReferenceShard shard = ReferenceShard.getShardNumberFromInterval(referenceWindowFunction.apply(c.element()));
                c.output(KV.of(shard, c.element()));
            }
        }).named("KeyReadByRefShard"));
        return keyReadByReferenceShard.apply(GroupByKey.<ReferenceShard, GATKRead>create());
    }
}

