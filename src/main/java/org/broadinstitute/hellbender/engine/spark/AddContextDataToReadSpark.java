package org.broadinstitute.hellbender.engine.spark;

import com.google.common.base.Function;
import com.google.common.collect.Iterators;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction2;
import org.apache.spark.api.java.function.PairFlatMapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.ReadContextData;
import org.broadinstitute.hellbender.engine.Shard;
import org.broadinstitute.hellbender.engine.ShardBoundary;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipList;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import scala.Tuple2;

import javax.annotation.Nullable;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * AddContextDataToRead pairs reference bases and overlapping variants with each GATKRead in the RDD input.
 * The variants are obtained from a local file (later a GCS Bucket). The reference bases come from the Google Genomics API.
 *
 * This transform is intended for direct use in pipelines.
 *
 * This transform will filter out any unmapped reads.
 *
 * The reference bases paired with each read can be customized by passing in a reference window function
 * inside the {@link org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource} argument to {@link #add}. See
 * {@link org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions} for examples.
 */
public class AddContextDataToReadSpark {
    /**
     * Add context data ({@link ReadContextData}) to reads.
     * @param ctx the Spark context
     * @param reads the coordinate-sorted reads
     * @param referenceSource the reference source
     * @param variants the coordinate-sorted variants
     * @param joinStrategy the strategy to use to join context data to reads
     * @param sequenceDictionary the sequence dictionary for the reads (only used for OVERLAPS_PARTITIONER join strategy, use null otherwise)
     * @param shardSize the maximum size of each shard, in bases (only used for OVERLAPS_PARTITIONER join strategy, use 0 otherwise)
     * @param shardPadding amount of extra context around each shard, in bases (only used for OVERLAPS_PARTITIONER join strategy, use 0 otherwise)
     * @return a RDD of read-context pairs, in coordinate-sorted order
     */
    public static JavaPairRDD<GATKRead, ReadContextData> add(
            final JavaSparkContext ctx,
            final JavaRDD<GATKRead> reads, final ReferenceMultiSource referenceSource,
            final JavaRDD<GATKVariant> variants, final JoinStrategy joinStrategy,
            final SAMSequenceDictionary sequenceDictionary,
            final int shardSize, final int shardPadding) {
        // TODO: this static method should not be filtering the unmapped reads.  To be addressed in another issue.
        JavaRDD<GATKRead> mappedReads = reads.filter(read -> ReadFilterLibrary.MAPPED.test(read));
        JavaPairRDD<GATKRead, Tuple2<Iterable<GATKVariant>, ReferenceBases>> withVariantsWithRef;
        if (joinStrategy.equals(JoinStrategy.BROADCAST)) {
            // Join Reads and Variants
            JavaPairRDD<GATKRead, Iterable<GATKVariant>> withVariants = BroadcastJoinReadsWithVariants.join(mappedReads, variants);
            // Join Reads with ReferenceBases
            withVariantsWithRef = BroadcastJoinReadsWithRefBases.addBases(referenceSource, withVariants);
        } else if (joinStrategy.equals(JoinStrategy.SHUFFLE)) {
            // Join Reads and Variants
            JavaPairRDD<GATKRead, Iterable<GATKVariant>> withVariants = ShuffleJoinReadsWithVariants.join(mappedReads, variants);
            // Join Reads with ReferenceBases
            withVariantsWithRef = ShuffleJoinReadsWithRefBases.addBases(referenceSource, withVariants);
        } else if (joinStrategy.equals(JoinStrategy.OVERLAPS_PARTITIONER)) {
            return addUsingOverlapsPartitioning(ctx, reads, referenceSource, variants, sequenceDictionary, shardSize, shardPadding);
        } else {
            throw new UserException("Unknown JoinStrategy");
        }
        return withVariantsWithRef.mapToPair(in -> new Tuple2<>(in._1(), new ReadContextData(in._2()._2(), in._2()._1())));
    }

    /**
     * Add context data ({@link ReadContextData}) to reads, using overlaps partitioning to avoid a shuffle.
     * @param ctx the Spark context
     * @param mappedReads the coordinate-sorted reads
     * @param referenceSource the reference source
     * @param variants the coordinate-sorted variants
     * @param sequenceDictionary the sequence dictionary for the reads
     * @param shardSize the maximum size of each shard, in bases
     * @param shardPadding amount of extra context around each shard, in bases
     * @return a RDD of read-context pairs, in coordinate-sorted order
     */
    private static JavaPairRDD<GATKRead, ReadContextData> addUsingOverlapsPartitioning(
            final JavaSparkContext ctx,
            final JavaRDD<GATKRead> mappedReads, final ReferenceMultiSource referenceSource,
            final JavaRDD<GATKVariant> variants, final SAMSequenceDictionary sequenceDictionary,
            final int shardSize, final int shardPadding) {
        final List<SimpleInterval> intervals = IntervalUtils.getAllIntervalsForReference(sequenceDictionary);
        // use unpadded shards (padding is only needed for reference bases)
        final List<ShardBoundary> intervalShards = intervals.stream()
                .flatMap(interval -> Shard.divideIntervalIntoShards(interval, shardSize, 0, sequenceDictionary).stream())
                .collect(Collectors.toList());

        final Broadcast<ReferenceMultiSource> bReferenceSource = ctx.broadcast(referenceSource);

        int maxLocatableSize = Math.min(shardSize, shardPadding);
        JavaRDD<Shard<GATKRead>> shardedReads = SparkSharder.shard(ctx, mappedReads, GATKRead.class, sequenceDictionary, intervalShards, maxLocatableSize);
        JavaRDD<Shard<GATKVariant>> shardedVariants = SparkSharder.shard(ctx, variants, GATKVariant.class, sequenceDictionary, intervalShards, maxLocatableSize);
        JavaRDD<Shard<GATKVariant>> repartitionedShardedVariants = SparkSharder.repartition(ctx, shardedVariants, mappedReads, sequenceDictionary, maxLocatableSize); // partitions now match those of sharded reads
        return shardedReads.zipPartitions(repartitionedShardedVariants, (FlatMapFunction2<Iterator<Shard<GATKRead>>, Iterator<Shard<GATKVariant>>, Tuple2<GATKRead, ReadContextData>>) (shardIterator, shardIterator2) -> {
            List<GATKVariant> variants1 = new ArrayList<>();
            while (shardIterator2.hasNext()) {
                Shard<GATKVariant> variantShard = shardIterator2.next();
                for (GATKVariant variant : variantShard) {
                    variants1.add(variant);
                }
            }
            OverlapDetector<GATKVariant> overlapDetector = OverlapDetector.create(variants1);
            return Iterators.concat(Iterators.transform(shardIterator, new Function<Shard<GATKRead>, Iterator<Tuple2<GATKRead, ReadContextData>>>() {
                private static final long serialVersionUID = 1L;
                @Nullable
                @Override
                public Iterator<Tuple2<GATKRead, ReadContextData>> apply(@Nullable Shard<GATKRead> shard) {
                    SimpleInterval paddedInterval = shard.getInterval().expandWithinContig(shardPadding, sequenceDictionary);
                    final ReferenceBases referenceBases;
                    try {
                        referenceBases = bReferenceSource.getValue().getReferenceBases(null, paddedInterval);
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                    Iterator<Tuple2<GATKRead, ReadContextData>> transform = Iterators.transform(shard.iterator(), new Function<GATKRead, Tuple2<GATKRead, ReadContextData>>() {
                        private static final long serialVersionUID = 1L;
                        @Nullable
                        @Override
                        public Tuple2<GATKRead, ReadContextData> apply(@Nullable GATKRead r) {
                            List<GATKVariant> overlappingVariants;
                            if (SimpleInterval.isValid(r.getContig(), r.getStart(), r.getEnd())) {
                                overlappingVariants = new ArrayList<>(overlapDetector.getOverlaps(r));
                            } else {
                                //Sometimes we have reads that do not form valid intervals (reads that do not consume any ref bases, eg CIGAR 61S90I
                                //In those cases, we'll just say that nothing overlaps the read
                                overlappingVariants = Collections.emptyList();
                            }
                            return new Tuple2<>(r, new ReadContextData(referenceBases, overlappingVariants));
                        }
                    });
                    // only include reads that start in the shard
                    return Iterators.filter(transform, r -> r._1().getStart() >= shard.getStart()
                            && r._1().getStart() <= shard.getEnd());
                }
            }));
        }).mapToPair(t -> t);
    }
}
