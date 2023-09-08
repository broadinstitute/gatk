package org.broadinstitute.hellbender.engine.spark;

import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.SparkFiles;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.DownsampleableSparkReadShard;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileStateRange;
import org.broadinstitute.hellbender.utils.downsampling.PositionalDownsampler;
import org.broadinstitute.hellbender.utils.downsampling.ReadsDownsampler;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import javax.annotation.Nullable;
import java.util.Iterator;
import java.util.List;
import java.util.function.Supplier;

/**
 * Find assembly regions from reads in a distributed Spark setting. There are two algorithms available, <i>fast</i>,
 * which looks for assembly regions in each read shard in parallel, and <i>strict</i>, which looks for assembly regions
 * in each contig in parallel. Fast mode may produce read shard boundary artifacts for assembly regions compared to the
 * walker version. Strict mode should be identical to the walker version, at the cost of increased runtime compared to
 * the fast version.
 */
public class FindAssemblyRegionsSpark {

    /**
     * Get an RDD of assembly regions for the given reads and intervals using the <i>fast</i> algorithm (looks for
     * assembly regions in each read shard in parallel).
     * @param ctx the Spark context
     * @param reads the coordinate-sorted reads
     * @param header the header for the reads
     * @param sequenceDictionary the sequence dictionary for the reads
     * @param referenceFileName the file name for the reference
     * @param features source of arbitrary features (may be null)
     * @param intervalShards the sharded intervals to find assembly regions for
     * @param assemblyRegionEvaluatorSupplierBroadcast evaluator used to determine whether a locus is active
     * @param shardingArgs the arguments for sharding reads
     * @param assemblyRegionArgs the arguments for finding assembly regions
     * @param shuffle whether to use a shuffle or not when sharding reads
     * @return an RDD of assembly regions
     */
    public static JavaRDD<AssemblyRegionWalkerContext> getAssemblyRegionsFast(
            final JavaSparkContext ctx,
            final JavaRDD<GATKRead> reads,
            final SAMFileHeader header,
            final SAMSequenceDictionary sequenceDictionary,
            final String referenceFileName,
            final FeatureManager features,
            final List<ShardBoundary> intervalShards,
            final Broadcast<Supplier<AssemblyRegionEvaluator>> assemblyRegionEvaluatorSupplierBroadcast,
            final AssemblyRegionReadShardArgumentCollection shardingArgs,
            final AssemblyRegionArgumentCollection assemblyRegionArgs,
            final boolean shuffle,
            final boolean trackPileups) {
        JavaRDD<Shard<GATKRead>> shardedReads = SparkSharder.shard(ctx, reads, GATKRead.class, sequenceDictionary, intervalShards, shardingArgs.readShardSize, shuffle);
        Broadcast<FeatureManager> bFeatureManager = features == null ? null : ctx.broadcast(features);
        return shardedReads.mapPartitions(getAssemblyRegionsFunctionFast(referenceFileName, bFeatureManager, header,
                assemblyRegionEvaluatorSupplierBroadcast, assemblyRegionArgs, trackPileups));
    }

    private static FlatMapFunction<Iterator<Shard<GATKRead>>, AssemblyRegionWalkerContext> getAssemblyRegionsFunctionFast(
            final String referenceFileName,
            final Broadcast<FeatureManager> bFeatureManager,
            final SAMFileHeader header,
            final Broadcast<Supplier<AssemblyRegionEvaluator>> supplierBroadcast,
            final AssemblyRegionArgumentCollection assemblyRegionArgs,
            final boolean trackPileups) {
        return (FlatMapFunction<Iterator<Shard<GATKRead>>, AssemblyRegionWalkerContext>) shardedReadIterator -> {
            final ReferenceDataSource reference = referenceFileName == null ? null : new ReferenceFileSource(IOUtils.getPath(SparkFiles.get(referenceFileName)));
            final FeatureManager features = bFeatureManager == null ? null : bFeatureManager.getValue();
            final AssemblyRegionEvaluator assemblyRegionEvaluator = supplierBroadcast.getValue().get(); // one AssemblyRegionEvaluator instance per Spark partition
            final ReadsDownsampler readsDownsampler = assemblyRegionArgs.maxReadsPerAlignmentStart > 0 ?
                    new PositionalDownsampler(assemblyRegionArgs.maxReadsPerAlignmentStart, header) : null;

            Iterator<Iterator<AssemblyRegionWalkerContext>> iterators = Utils.stream(shardedReadIterator)
                    .map(shardedRead -> new ShardToMultiIntervalShardAdapter<>(
                            new DownsampleableSparkReadShard(
                                    new ShardBoundary(shardedRead.getInterval(), shardedRead.getPaddedInterval()), shardedRead, readsDownsampler)))
                    .map(downsampledShardedRead -> {
                        final Iterator<AssemblyRegion> assemblyRegionIter = new AssemblyRegionIterator(
                                new ShardToMultiIntervalShardAdapter<>(downsampledShardedRead),
                                header, reference, features, assemblyRegionEvaluator, assemblyRegionArgs, trackPileups);
                        return Utils.stream(assemblyRegionIter).map(assemblyRegion ->
                                new AssemblyRegionWalkerContext(assemblyRegion,
                                        new ReferenceContext(reference, assemblyRegion.getPaddedSpan()),
                                        new FeatureContext(features, assemblyRegion.getPaddedSpan()))).iterator();
                    }).iterator();
            return Iterators.concat(iterators);
        };
    }

    /**
     * Get an RDD of assembly regions for the given reads and intervals using the <i>strict</i> algorithm (looks for
     * assembly regions in each contig in parallel).
     * @param ctx the Spark context
     * @param reads the coordinate-sorted reads
     * @param header the header for the reads
     * @param sequenceDictionary the sequence dictionary for the reads
     * @param referenceFileName the file name for the reference
     * @param features source of arbitrary features (may be null)
     * @param intervalShards the sharded intervals to find assembly regions for
     * @param assemblyRegionEvaluatorSupplierBroadcast evaluator used to determine whether a locus is active
     * @param shardingArgs the arguments for sharding reads
     * @param assemblyRegionArgs the arguments for finding assembly regions
     * @param shuffle whether to use a shuffle or not when sharding reads
     * @return an RDD of assembly regions
     */
    public static JavaRDD<AssemblyRegionWalkerContext> getAssemblyRegionsStrict(
            final JavaSparkContext ctx,
            final JavaRDD<GATKRead> reads,
            final SAMFileHeader header,
            final SAMSequenceDictionary sequenceDictionary,
            final String referenceFileName,
            final FeatureManager features,
            final List<ShardBoundary> intervalShards,
            final Broadcast<Supplier<AssemblyRegionEvaluator>> assemblyRegionEvaluatorSupplierBroadcast,
            final AssemblyRegionReadShardArgumentCollection shardingArgs,
            final AssemblyRegionArgumentCollection assemblyRegionArgs,
            final boolean shuffle) {
        JavaRDD<Shard<GATKRead>> shardedReads = SparkSharder.shard(ctx, reads, GATKRead.class, sequenceDictionary, intervalShards, shardingArgs.readShardSize, shuffle);
        Broadcast<FeatureManager> bFeatureManager = features == null ? null : ctx.broadcast(features);

        // 1. Calculate activity for each locus in the desired intervals, in parallel.
        JavaRDD<ActivityProfileStateRange> activityProfileStates = shardedReads.mapPartitions(getActivityProfileStatesFunction(referenceFileName, bFeatureManager, header,
                assemblyRegionEvaluatorSupplierBroadcast, assemblyRegionArgs));

        // 2. Group by contig. We need to do this so we can perform the band pass filter over the whole contig, so we
        // produce assembly regions that are identical to those produced by AssemblyRegionWalker.
        // This step requires a shuffle, but the amount of data in the ActivityProfileStateRange should be small, so it
        // should not be prohibitive.
        JavaPairRDD<String, Iterable<ActivityProfileStateRange>> contigToGroupedStates = activityProfileStates
                .keyBy((Function<ActivityProfileStateRange, String>) range -> range.getContig())
                .groupByKey();

        // 3. Run the band pass filter to find AssemblyRegions. The filtering is fairly cheap, so should be fast
        // even though it has to scan a whole contig. Note that we *don't* fill in reads here, since after we have found
        // the assembly regions we want to do assembly using the full resources of the cluster. So if we have
        // very small assembly region objects, then we can repartition them for redistribution across the cluster,
        // at which points the reads can be filled in. (See next step.)
        JavaRDD<ReadlessAssemblyRegion> readlessAssemblyRegions = contigToGroupedStates
                .flatMap(getReadlessAssemblyRegionsFunction(header, assemblyRegionArgs));
        // repartition to distribute the data evenly across the cluster again
        readlessAssemblyRegions = readlessAssemblyRegions.repartition(readlessAssemblyRegions.getNumPartitions());

        // 4. Fill in the reads. Each shard is an assembly region, with its overlapping reads.
        JavaRDD<Shard<GATKRead>> assemblyRegionShardedReads = SparkSharder.shard(ctx, reads, GATKRead.class, header.getSequenceDictionary(), readlessAssemblyRegions, shardingArgs.readShardSize);

        // 5. Convert shards to assembly regions. Reads downsampling is done again here. Note it will only be
        // consistent with the downsampling done in step 1 when https://github.com/broadinstitute/gatk/issues/5437 is in.
        JavaRDD<AssemblyRegion> assemblyRegions = assemblyRegionShardedReads.mapPartitions((FlatMapFunction<Iterator<Shard<GATKRead>>, AssemblyRegion>) shardedReadIterator -> {
            final ReadsDownsampler readsDownsampler = assemblyRegionArgs.maxReadsPerAlignmentStart > 0 ?
                    new PositionalDownsampler(assemblyRegionArgs.maxReadsPerAlignmentStart, header) : null;
            return Utils.stream(shardedReadIterator)
                    .map(shardedRead -> toAssemblyRegion(shardedRead, header, readsDownsampler)).iterator();
        });

        // 6. Add reference and feature context.
        return assemblyRegions.mapPartitions(getAssemblyRegionWalkerContextFunction(referenceFileName, bFeatureManager));
    }

    private static FlatMapFunction<Iterator<Shard<GATKRead>>, ActivityProfileStateRange> getActivityProfileStatesFunction(
            final String referenceFileName,
            final Broadcast<FeatureManager> bFeatureManager,
            final SAMFileHeader header,
            final Broadcast<Supplier<AssemblyRegionEvaluator>> supplierBroadcast,
            final AssemblyRegionArgumentCollection assemblyRegionArgs) {
        return (FlatMapFunction<Iterator<Shard<GATKRead>>, ActivityProfileStateRange>) shardedReadIterator -> {
            final ReferenceDataSource reference = referenceFileName == null ? null : new ReferenceFileSource(IOUtils.getPath(SparkFiles.get(referenceFileName)));
            final FeatureManager features = bFeatureManager == null ? null : bFeatureManager.getValue();
            final AssemblyRegionEvaluator assemblyRegionEvaluator = supplierBroadcast.getValue().get(); // one AssemblyRegionEvaluator instance per Spark partition
            
            return Utils.stream(shardedReadIterator)
                    .map(shardedRead -> {
                        final ReadsDownsampler readsDownsampler = assemblyRegionArgs.maxReadsPerAlignmentStart > 0 ?
                                new PositionalDownsampler(assemblyRegionArgs.maxReadsPerAlignmentStart, header) : null;
                        return new ShardToMultiIntervalShardAdapter<>(
                                new DownsampleableSparkReadShard(
                                        new ShardBoundary(shardedRead.getInterval(), shardedRead.getPaddedInterval()), shardedRead, readsDownsampler));
                    })
                    .map(shardedRead -> {
                        final Iterator<ActivityProfileState> activityProfileStateIter = new ActivityProfileStateIterator(
                                new ShardToMultiIntervalShardAdapter<>(shardedRead),
                                header, reference, features, assemblyRegionEvaluator
                        );
                        return new ActivityProfileStateRange(shardedRead, activityProfileStateIter);
                    }).iterator();
        };
    }

    private static FlatMapFunction<Tuple2<String, Iterable<ActivityProfileStateRange>>, ReadlessAssemblyRegion> getReadlessAssemblyRegionsFunction(
            final SAMFileHeader header,
            final AssemblyRegionArgumentCollection assemblyRegionArgs) {
        return (FlatMapFunction<Tuple2<String, Iterable<ActivityProfileStateRange>>, ReadlessAssemblyRegion>) iter ->
                Iterators.transform(
                        new AssemblyRegionFromActivityProfileStateIterator(
                                ActivityProfileStateRange.toIteratorActivityProfileState(iter._2.iterator()),
                                header,
                                assemblyRegionArgs.minAssemblyRegionSize,
                                assemblyRegionArgs.maxAssemblyRegionSize,
                                assemblyRegionArgs.assemblyRegionPadding,
                                assemblyRegionArgs.activeProbThreshold,
                                assemblyRegionArgs.maxProbPropagationDistance), new com.google.common.base.Function<AssemblyRegion, ReadlessAssemblyRegion>() {
                            @Nullable
                            @Override
                            public ReadlessAssemblyRegion apply(@Nullable AssemblyRegion input) {
                                return new ReadlessAssemblyRegion(input);
                            }
                        });
    }

    private static AssemblyRegion toAssemblyRegion(Shard<GATKRead> shard, SAMFileHeader header, ReadsDownsampler readsDownsampler) {
        Shard<GATKRead> downsampledShardedRead =
                new DownsampleableSparkReadShard(
                        new ShardBoundary(shard.getInterval(), shard.getPaddedInterval()), shard, readsDownsampler);

        // TODO: interfaces could be improved to avoid casting
        ReadlessAssemblyRegion readlessAssemblyRegion = (ReadlessAssemblyRegion) ((ShardBoundaryShard<GATKRead>) shard).getShardBoundary();
        int extension = Math.max(shard.getInterval().getStart() - shard.getPaddedInterval().getStart(), shard.getPaddedInterval().getEnd() - shard.getInterval().getEnd());
        AssemblyRegion assemblyRegion = new AssemblyRegion(shard.getInterval(), readlessAssemblyRegion.isActive(), extension, header);
        assemblyRegion.addAll(Lists.newArrayList(downsampledShardedRead));
        return assemblyRegion;
    }

    private static FlatMapFunction<Iterator<AssemblyRegion>, AssemblyRegionWalkerContext> getAssemblyRegionWalkerContextFunction(
            final String referenceFileName,
            final Broadcast<FeatureManager> bFeatureManager) {

        return (FlatMapFunction<Iterator<AssemblyRegion>, AssemblyRegionWalkerContext>) assemblyRegionIter -> {
            final ReferenceDataSource reference = referenceFileName == null ? null : new ReferenceFileSource(IOUtils.getPath(SparkFiles.get(referenceFileName)));
            final FeatureManager features = bFeatureManager == null ? null : bFeatureManager.getValue();
            return Utils.stream(assemblyRegionIter).map(assemblyRegion ->
                    new AssemblyRegionWalkerContext(assemblyRegion,
                            new ReferenceContext(reference, assemblyRegion.getPaddedSpan()),
                            new FeatureContext(features, assemblyRegion.getPaddedSpan()))).iterator();
        };
    }
}
