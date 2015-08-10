package org.broadinstitute.hellbender.engine.dataflow.transforms.composite;

import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefAPIMetadata;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefWindowFunctions;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceDataflowSource;
import org.broadinstitute.hellbender.engine.dataflow.transforms.KeyReadsByRefShard;
import org.broadinstitute.hellbender.engine.dataflow.transforms.RefBasesFromAPI;
import org.broadinstitute.hellbender.engine.dataflow.transforms.PairReadWithRefBases;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceShard;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

/**
 * RefBasesForReads queries the Google Genomics API for reference bases overlapping all of the reads.
 *
 * step 1: key reads by reference shards
 *
 * |--------- shard 0 ----------|---------- shard 1 ----------|--------- shard 2 ----------|--------- shard 3 ----------|
 *           |------ read a -----|                             |-- read b --|   |---- read c ----|
 *
 * step 2: group reads by the shard they start in
 *
 *  |--------- shard 0 ----------|
 *            |------ read a -----|
 *
 *
 *  |--------- shard 2 ----------|
 *   |-- read b --|   |---- read c ----|
 *
 *  step 3: query the Google Genomics API for all bases needed for each shard
 *

 * |--- ref bases 1 ---|        |--------- ref bases 2 -----------|
 * |------ read a -----|        |-- read b --|   |---- read c ----|
 *
 *  step 4: pair the ref bases needed for each read with the read
 *
 * |------ read a -----|        |-- read b --|   |---- read c ----|
 * |-- ref bases 1a ---|        |ref bases 2b|   |- ref bases 2c -|
 *
 * or in code,
 *  KV<read a, ref bases 1a>
 *  KV<read b, ref bases 2b>
 *  KV<read c, ref bases 2c>
 *
 * The reference bases paired with each read can be customized by passing in a reference window function
 * inside the {@link ReferenceDataflowSource} argument to {@link #addBases}. See {@link RefWindowFunctions} for examples.
 */
public class RefBasesForReads {
    public static PCollection<KV<GATKRead, ReferenceBases>> addBases(PCollection<GATKRead> pReads, ReferenceDataflowSource referenceDataflowSource) {
        PCollection<KV<ReferenceShard, Iterable<GATKRead>>> shardAndRead = pReads.apply(new KeyReadsByRefShard(referenceDataflowSource.getReferenceWindowFunction()));
        PCollection<KV<ReferenceBases, Iterable<GATKRead>>> groupedReads =
                RefBasesFromAPI.getBasesForShard(shardAndRead, referenceDataflowSource);
        return groupedReads.apply(new PairReadWithRefBases(referenceDataflowSource.getReferenceWindowFunction()));
    }
}
