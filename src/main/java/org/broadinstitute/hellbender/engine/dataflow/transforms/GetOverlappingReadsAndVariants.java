package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.transforms.*;
import com.google.cloud.dataflow.sdk.transforms.join.CoGbkResult;
import com.google.cloud.dataflow.sdk.transforms.join.CoGroupByKey;
import com.google.cloud.dataflow.sdk.transforms.join.KeyedPCollectionTuple;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.TupleTag;
import org.broadinstitute.hellbender.engine.VariantShard;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.Variant;

/**
 * PairReadsAndVariants takes two PCollections (GATKRead and Variant) and returns a PCollection with a
 * (GATKRead,Variant) pair for every read and variant that overlap. This means that a read or variant may be present
 * multiple times in the output. Also, there may be duplicate (GATKRead,Variant) pairs in the output.
 * Currently, all reads must be mapped.
 *
 * The function works by creating a single PCollection of KV where GATKRead is the key and Variant is the value.
 * We do this join by sharding both collections by "variant shard" and checking for overlap on each "shard."
 *
 * step 1: key reads and variants by shard
 * |---- shard 0 -----|---- shard 1 -----|---- shard 2 -----|---- shard 3 -----|---- shard 4 -----|
 *     |---------- read a ---------|               |----- read b ------|
 *   |- variant 1 -|    |- variant 2 -|               |- variant 3 -|
 *
 * step 2: shard read and variant by variant shard
 *                      |---- shard 0 -----|
 *                          |---------- read a ---------|
 *                        |- variant 1 -|
 *
 *
 *                      |---- shard 1 -----|
 *       |---------- read a ---------|
 *                        |- variant 2 -|
 *
 *
 *                      |---- shard 2 -----|
 *                                |----- read b ------|
 *                                   |- variant 3 -|
 *
 *                      |---- shard 3 -----|
 *             |----- read b ------|
 *                |- variant 3 -|
 *
 * step 3: pair reads and variants
 * KV<read a, variant 1> // from shard 0
 * KV<read a, variant 2> // from shard 1
 * KV<read b, variant 3> // from shard 2
 * KV<read b, variant 3> // from shard 3
 *
 */
public class GetOverlappingReadsAndVariants {
    public static PCollection<KV<GATKRead, Variant>> pair(PCollection<GATKRead> pRead, PCollection<Variant> pVariant) {

        PCollection<KV<VariantShard, GATKRead>> vkReads = pRead.apply(new KeyReadsByOverlappingVariantShard());
        PCollection<KV<VariantShard, Variant>> vkVariants =
                pVariant.apply(new KeyVariantByOverlappingVariantShard());

        // GroupBy VariantShard
        final TupleTag<Variant> variantTag = new TupleTag<>();
        final TupleTag<GATKRead> readTag = new TupleTag<>();
        PCollection<KV<VariantShard, CoGbkResult>> coGbkInput = KeyedPCollectionTuple
                .of(variantTag, vkVariants)
                .and(readTag, vkReads).apply(CoGroupByKey.<VariantShard>create());

        // GroupBy Read
        return coGbkInput.apply(ParDo.of(
                new DoFn<KV<VariantShard, CoGbkResult>, KV<GATKRead, Variant>>() {
                    private static final long serialVersionUID = 1L;
                    @Override
                    public void processElement(ProcessContext c) throws Exception {
                        Iterable<Variant> kVariants = c.element().getValue().getAll(variantTag);
                        Iterable<GATKRead> kReads = c.element().getValue().getAll(readTag);
                        // Compute overlap.
                        for (GATKRead r : kReads) {
                            SimpleInterval readInterval = new SimpleInterval(r);
                            for (Variant v : kVariants) {
                                if (readInterval.overlaps(v)) {
                                    c.output(KV.of(r, v));
                                }
                            }
                        }
                    }
                })).setName("PairReadsAndVariants_GroupByRead");

    }
}
