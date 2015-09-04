package org.broadinstitute.hellbender.engine.spark;

import com.google.cloud.dataflow.sdk.transforms.SerializableFunction;
import com.google.common.collect.Lists;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceDataflowSource;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceShard;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import scala.Tuple2;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

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
 * inside the {@link ReferenceDataflowSource} argument to {@link #addBases}. See {@link org.broadinstitute.hellbender.engine.dataflow.datasources.RefWindowFunctions} for examples.
 */
public class JoinReadsWithRefBases {
    public static JavaPairRDD<GATKRead, ReferenceBases> addBases(final ReferenceDataflowSource referenceDataflowSource,
                                                                 final JavaRDD<GATKRead> reads) {
        SerializableFunction<GATKRead, SimpleInterval> windowFunction = referenceDataflowSource.getReferenceWindowFunction();

        JavaPairRDD<ReferenceShard, GATKRead> shardRead = reads.mapToPair(gatkRead -> {
            ReferenceShard shard = ReferenceShard.getShardNumberFromInterval(windowFunction.apply(gatkRead));
            return new Tuple2<>(shard, gatkRead);
        });

        JavaPairRDD<ReferenceShard, Iterable<GATKRead>> shardiRead = shardRead.groupByKey();

        return shardiRead.flatMapToPair(in -> {
            List<Tuple2<GATKRead, ReferenceBases>> out = Lists.newArrayList();
            Iterable<GATKRead> iReads = in._2();

            // Apply the reference window function to each read to produce a set of intervals representing
            // the desired reference bases for each read.
            final List<SimpleInterval> readWindows = StreamSupport.stream(iReads.spliterator(), false).map(read -> windowFunction.apply(read)).collect(Collectors.toList());

            SimpleInterval interval = SimpleInterval.getSpanningInterval(readWindows);
            ReferenceBases bases = referenceDataflowSource.getReferenceBases(null, interval);
            for (GATKRead r : iReads) {
                final ReferenceBases subset = bases.getSubset(windowFunction.apply(r));
                out.add(new Tuple2<>(r, subset));
            }
            return out;
        });
    }
}