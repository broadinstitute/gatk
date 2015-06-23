package org.broadinstitute.hellbender.engine.dataflow.transforms.composite;

import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.transforms.PairReadsAndVariants;
import org.broadinstitute.hellbender.engine.dataflow.transforms.RemoveDuplicateReadVariantPairs;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.Variant;

/**
 * KeyVariantsByRead takes a collection of Reads and Variants and pairs each variant with every read it overlaps.
 * The final result is an iterable of Variants keyed by GATKRead. The join is accomplished by sharding the reads and
 * variants by VariantShard (chunks of contig) and checking matches per shard.
 *
 * See the diagram in PairReadsAndVariants.java for more details.
 */
public class KeyVariantsByRead {
    public static PCollection<KV<GATKRead, Iterable<Variant>>> Key(PCollection<Variant> pVariants, PCollection<GATKRead> pReads) {
        PCollection<KV<GATKRead, Variant>> readVariants = PairReadsAndVariants.pair(pReads, pVariants);

        // At this point, we ALMOST have what we want, but we need to remove duplicates KV<Read, Variant> pairs.
        // And we need to group by Read. Both of these require having deterministic coding, so we need to switch to
        // UUIDS.
        return readVariants.apply(new RemoveDuplicateReadVariantPairs());
    }
}
