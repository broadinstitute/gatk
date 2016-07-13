package org.broadinstitute.hellbender.engine.spark;

import com.google.common.collect.Lists;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.VariantShard;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import scala.Tuple2;

import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

/**
 * PairReadsAndVariants takes two RDDs (GATKRead and Variant) and returns an RDD with a
 * (GATKRead,Iterable<Variant>) for every read and all the variants that overlap. We do this by first
 * Making Tuple2s of GATKRead,Variant, which means that a read or variant may be present
 * multiple times in the output. Also, there may be duplicate (GATKRead,Variant) pairs in the output.
 * Currently, all reads must be mapped. We then aggregrate by key to remove duplicate variants.
 *
 * The function works by creating a single RDD of Tuple2 where GATKRead is the key and Variant is the value.
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
 * Tuple2<read a, variant 1> // from shard 0
 * Tuple2<read a, variant 2> // from shard 1
 * Tuple2<read b, variant 3> // from shard 2
 * Tuple2<read b, variant 3> // from shard 3
 *
 * step 4: aggregate by key
 * Tuple2<read a, <variant 1, variant2>>
 * Tuple2<read b, <variant 3>>
 */
public class ShuffleJoinReadsWithVariants {
    public static JavaPairRDD<GATKRead, Iterable<GATKVariant>> join(
            final JavaRDD<GATKRead> reads, final JavaRDD<GATKVariant> variants) {

        JavaPairRDD<VariantShard, GATKRead> readsWShards = pairReadsWithVariantShards(reads);

        JavaPairRDD<VariantShard, GATKVariant> variantsWShards = pairVariantsWithVariantShards(variants);

        // generate read-variant pairs; however, the reads are replicated for each overlapping pair
        JavaPairRDD<GATKRead, GATKVariant> allPairs = pairReadsWithVariants(readsWShards, variantsWShards);

        // we group together all variants for each unique GATKRead.  As we combine through the Variants, they are added
        // to a HashSet that get continually merged together
        return allPairs.aggregateByKey(new LinkedHashSet<>(), (vs, v) -> {
            if (v != null) { // pairReadsWithVariants can produce null variant
                ((Set<GATKVariant>) vs).add(v);
            }
            return vs;
        }, (vs1, vs2) -> {
            ((Set<GATKVariant>) vs1).addAll((Set<GATKVariant>) vs2);
            return vs1;
        });
    }

    private static JavaPairRDD<VariantShard, GATKRead> pairReadsWithVariantShards(final JavaRDD<GATKRead> reads) {
        return reads.flatMapToPair(gatkRead -> {
            List<VariantShard> shards = VariantShard.getVariantShardsFromInterval(gatkRead);
            List<Tuple2<VariantShard, GATKRead>> out = Lists.newArrayList();
            for (VariantShard shard : shards) {
                out.add(new Tuple2<>(shard, gatkRead));
            }
            return out.iterator();
        });
    }

    private static JavaPairRDD<VariantShard, GATKVariant> pairVariantsWithVariantShards(final JavaRDD<GATKVariant> variants) {
        return  variants.flatMapToPair(variant -> {
            List<VariantShard> shards = VariantShard.getVariantShardsFromInterval(variant);
            List<Tuple2<VariantShard, GATKVariant>> out = Lists.newArrayList();
            for (VariantShard shard : shards) {
                out.add(new Tuple2<>(shard, variant));
            }
            return out.iterator();
        });
    }
    private static JavaPairRDD<GATKRead, GATKVariant> pairReadsWithVariants(final JavaPairRDD<VariantShard, GATKRead> readsWShards,
                                                                            final  JavaPairRDD<VariantShard, GATKVariant> variantsWShards) {
        JavaPairRDD<VariantShard, Tuple2<Iterable<GATKRead>, Iterable<GATKVariant>>> cogroup = readsWShards.cogroup(variantsWShards);

        return cogroup.flatMapToPair(cogroupValue -> {
            Iterable<GATKRead> iReads = cogroupValue._2()._1();
            Iterable<GATKVariant> iVariants = cogroupValue._2()._2();

            List<Tuple2<GATKRead, GATKVariant>> out = Lists.newArrayList();
            // For every read, find every overlapping variant.
            for (GATKRead r : iReads) {
                boolean foundVariants = false;
                SimpleInterval interval = new SimpleInterval(r);
                for (GATKVariant v : iVariants) {
                    if (interval.overlaps(v)) {
                        foundVariants = true;
                        out.add(new Tuple2<>(r, v));
                    }
                }
                // If no variants are found, we still want to output the read.
                if (!foundVariants) {
                    out.add(new Tuple2<>(r, null));
                }
            }
            return out.iterator();
        });
    }
}
