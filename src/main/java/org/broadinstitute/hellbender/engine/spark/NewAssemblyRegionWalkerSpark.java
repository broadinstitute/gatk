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
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.tools.DownsampleableSparkReadShard;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.downsampling.PositionalDownsampler;
import org.broadinstitute.hellbender.utils.downsampling.ReadsDownsampler;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import javax.annotation.Nullable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.function.Supplier;
import java.util.stream.Collectors;

/**
 * A Spark version of {@link AssemblyRegionWalker} that produces identical assembly regions, unlike
 * {@link AssemblyRegionWalkerSpark} which can produce different assembly regions due to shard end artifacts.
 * Subclasses should implement {@link #processAssemblyRegions(JavaRDD, JavaSparkContext)}.
 */
public abstract class NewAssemblyRegionWalkerSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    public final AssemblyRegionReadShardArgumentCollection shardingArgs = new AssemblyRegionReadShardArgumentCollection();

    @ArgumentCollection
    public final AssemblyRegionArgumentCollection assemblyRegionArgs = getAssemblyRegionArgumentCollection();

    /**
     * @return a subclass of {@link AssemblyRegionArgumentCollection} with the default values filled in.
     */
    protected abstract AssemblyRegionArgumentCollection getAssemblyRegionArgumentCollection();

    /**
     * subclasses can override this to control if reads with deletions should be included in isActive pileups
     */
    protected abstract boolean includeReadsWithDeletionsInIsActivePileups();

    @Argument(doc = "whether to use the shuffle implementation or not", shortName = "shuffle", fullName = "shuffle", optional = true)
    public boolean shuffle = false;

    private String referenceFileName;

    @Override
    public final boolean requiresReads() { return true; }

    @Override
    public final boolean requiresReference() { return true; }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> defaultFilters = new ArrayList<>(2);
        defaultFilters.add(new WellformedReadFilter());
        defaultFilters.add(new ReadFilterLibrary.MappedReadFilter());
        return defaultFilters;
    }

    /**
     * @return The evaluator to be used to determine whether each locus is active or not. Must be implemented by tool authors.
     *         The results of this per-locus evaluator are used to determine the bounds of each active and inactive region.
     */
    public abstract AssemblyRegionEvaluator assemblyRegionEvaluator();

    /**
     * Tools that use an evaluator that is expensive to create, and/or that is not compatible with Spark broadcast, can
     * override this method to return a broadcast of a supplier of the evaluator. The supplier will be invoked once for
     * each Spark partition, thus each partition will have its own evaluator instance.
     */
    protected Broadcast<Supplier<AssemblyRegionEvaluator>> assemblyRegionEvaluatorSupplierBroadcast(final JavaSparkContext ctx) {
        return assemblyRegionEvaluatorSupplierBroadcastFunction(ctx, assemblyRegionEvaluator());
    }

    private static Broadcast<Supplier<AssemblyRegionEvaluator>> assemblyRegionEvaluatorSupplierBroadcastFunction(final JavaSparkContext ctx, final AssemblyRegionEvaluator assemblyRegionEvaluator) {
        Supplier<AssemblyRegionEvaluator> supplier = () -> assemblyRegionEvaluator;
        return ctx.broadcast(supplier);
    }

    private List<ShardBoundary> intervalShards;

    /**
     * Note that this sets {@code intervalShards} as a side effect, in order to add padding to the intervals.
     */
    @Override
    protected List<SimpleInterval> editIntervals(List<SimpleInterval> rawIntervals) {
        SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        List<SimpleInterval> intervals = rawIntervals == null ? IntervalUtils.getAllIntervalsForReference(sequenceDictionary) : rawIntervals;
        intervalShards = intervals.stream()
                .flatMap(interval -> Shard.divideIntervalIntoShards(interval, shardingArgs.readShardSize, shardingArgs.readShardPadding, sequenceDictionary).stream())
                .collect(Collectors.toList());
        List<SimpleInterval> paddedIntervalsForReads =
                intervals.stream().map(interval -> interval.expandWithinContig(shardingArgs.readShardPadding, sequenceDictionary)).collect(Collectors.toList());
        return paddedIntervalsForReads;
    }

    /**
     * Loads assembly regions and the corresponding reference and features into a {@link JavaRDD} for the intervals specified.
     *
     * If no intervals were specified, returns all the assembly regions.
     *
     * @return all assembly regions as a {@link JavaRDD}, bounded by intervals if specified.
     */
    protected JavaRDD<AssemblyRegionWalkerContext> getAssemblyRegions(JavaSparkContext ctx) {
        SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        return getAssemblyRegions(ctx, getReads(), getHeaderForReads(), sequenceDictionary, referenceFileName, features,
                intervalShards, assemblyRegionEvaluatorSupplierBroadcast(ctx), shardingArgs, assemblyRegionArgs,
                includeReadsWithDeletionsInIsActivePileups(), shuffle);
    }

    protected static JavaRDD<AssemblyRegionWalkerContext> getAssemblyRegions(
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
            final boolean includeReadsWithDeletionsInIsActivePileups,
            final boolean shuffle) {
        JavaRDD<Shard<GATKRead>> shardedReads = SparkSharder.shard(ctx, reads, GATKRead.class, sequenceDictionary, intervalShards, shardingArgs.readShardSize, shuffle);
        Broadcast<FeatureManager> bFeatureManager = features == null ? null : ctx.broadcast(features);

        // 1. Calculate activity for each locus in the desired intervals, in parallel.
        JavaRDD<ActivityProfileState> activityProfileStates = shardedReads.mapPartitions(getActivityProfileStatesFunction(referenceFileName, bFeatureManager, header,
                assemblyRegionEvaluatorSupplierBroadcast, assemblyRegionArgs, includeReadsWithDeletionsInIsActivePileups));

        // 2. Group by contig. We need to do this so we can perform the band pass filter over the whole contig, so we
        // produce assembly regions that are identical to those produced by AssemblyRegionWalker.
        // This step requires a shuffle, but the amount of data in the ActivityProfileState should be small, so it
        // should not be prohibitive (it might also be optimized).
        JavaPairRDD<String, Iterable<ActivityProfileState>> contigToGroupedStates = activityProfileStates
                .keyBy((Function<ActivityProfileState, String>) state -> state.getLoc().getContig())
                .groupByKey();

        // 3. Run the band pass filter to find AssemblyRegions. The filtering is fairly cheap, so should be fast
        // even though it has to scan a whole contig. Note that we *don't* fill in reads here, since after we have found
        // the assembly regions we want to do assembly using the full resources of the cluster. So if we have
        // very small assembly region objects, then we can easily collect them on the driver (or repartition them)
        // for redistribution across the cluster, at which points the reads can be filled in. (See next two steps.)
        // 3. Run the band pass filter to find AssemblyRegions.
        JavaRDD<ReadlessAssemblyRegion> readlessAssemblyRegions = contigToGroupedStates
                .flatMap(getReadlessAssemblyRegionsFunction(header, assemblyRegionArgs, includeReadsWithDeletionsInIsActivePileups));

        // 4. Pull the assembly region boundaries down to the driver, so we can fill in reads.
        List<ShardBoundary> assemblyRegionBoundaries = readlessAssemblyRegions
                .map((Function<ReadlessAssemblyRegion, ShardBoundary>) NewAssemblyRegionWalkerSpark::toShardBoundary)
                .collect();

        // 5. Fill in the reads. Each shard is an assembly region, with its overlapping reads.
        JavaRDD<Shard<GATKRead>> assemblyRegionShardedReads = SparkSharder.shard(ctx, reads, GATKRead.class, header.getSequenceDictionary(), assemblyRegionBoundaries, shardingArgs.readShardSize);

        // 6. Convert shards to assembly regions.
        JavaRDD<AssemblyRegion> assemblyRegions = assemblyRegionShardedReads.map((Function<Shard<GATKRead>, AssemblyRegion>) shard -> toAssemblyRegion(shard, header));

        // 7. Add reference and feature context.
        return assemblyRegions.mapPartitions(getAssemblyRegionWalkerContextFunction(referenceFileName, bFeatureManager));
    }

    private static FlatMapFunction<Iterator<Shard<GATKRead>>, ActivityProfileState> getActivityProfileStatesFunction(
            final String referenceFileName,
            final Broadcast<FeatureManager> bFeatureManager,
            final SAMFileHeader header,
            final Broadcast<Supplier<AssemblyRegionEvaluator>> supplierBroadcast,
            final AssemblyRegionArgumentCollection assemblyRegionArgs,
            final boolean includeReadsWithDeletionsInIsActivePileups) {
        return (FlatMapFunction<Iterator<Shard<GATKRead>>, ActivityProfileState>) shardedReadIterator -> {
            ReferenceDataSource reference = referenceFileName == null ? null : new ReferenceFileSource(IOUtils.getPath(SparkFiles.get(referenceFileName)));
            final FeatureManager features = bFeatureManager == null ? null : bFeatureManager.getValue();
            AssemblyRegionEvaluator assemblyRegionEvaluator = supplierBroadcast.getValue().get(); // one AssemblyRegionEvaluator instance per Spark partition
            final ReadsDownsampler readsDownsampler = assemblyRegionArgs.maxReadsPerAlignmentStart > 0 ?
                    new PositionalDownsampler(assemblyRegionArgs.maxReadsPerAlignmentStart, header) : null;

            Iterator<Iterator<ActivityProfileState>> iterators = Utils.stream(shardedReadIterator)
                    .map(shardedRead -> new ShardToMultiIntervalShardAdapter<>(shardedRead))
                    // TODO: reinstate downsampling (not yet working)
//                            new DownsampleableSparkReadShard(
//                                    new ShardBoundary(shardedRead.getInterval(), shardedRead.getPaddedInterval()), shardedRead, readsDownsampler)))
                    .map(shardedRead -> {
                final Iterator<ActivityProfileState> activityProfileStateIter = new ActivityProfileStateIterator(
                        new ShardToMultiIntervalShardAdapter<>(shardedRead),
                        header, reference, features, assemblyRegionEvaluator,
                        includeReadsWithDeletionsInIsActivePileups);
                return activityProfileStateIter;
            }).iterator();
            return Iterators.concat(iterators);
        };
    }

    private static FlatMapFunction<Tuple2<String, Iterable<ActivityProfileState>>, ReadlessAssemblyRegion> getReadlessAssemblyRegionsFunction(
            final SAMFileHeader header,
            final AssemblyRegionArgumentCollection assemblyRegionArgs,
            final boolean includeReadsWithDeletionsInIsActivePileups) {
        return (FlatMapFunction<Tuple2<String, Iterable<ActivityProfileState>>, ReadlessAssemblyRegion>) iter ->
                Iterators.transform(
                        new AssemblyRegionFromActivityProfileStateIterator(
                                iter._2.iterator(),
                                header,
                                assemblyRegionArgs.minAssemblyRegionSize,
                                assemblyRegionArgs.maxAssemblyRegionSize,
                                assemblyRegionArgs.assemblyRegionPadding,
                                assemblyRegionArgs.activeProbThreshold,
                                assemblyRegionArgs.maxProbPropagationDistance,
                                includeReadsWithDeletionsInIsActivePileups), new com.google.common.base.Function<AssemblyRegion, ReadlessAssemblyRegion>() {
                            @Nullable
                            @Override
                            public ReadlessAssemblyRegion apply(@Nullable AssemblyRegion input) {
                                return new ReadlessAssemblyRegion(input);
                            }
                        });
    }

    public static ShardBoundary toShardBoundary(ReadlessAssemblyRegion assemblyRegion) {
        return assemblyRegion;
    }

    private static AssemblyRegion toAssemblyRegion(Shard<GATKRead> shard, SAMFileHeader header) {
        // TODO: interfaces could be better designed to avoid casting
        ReadlessAssemblyRegion readlessAssemblyRegion = (ReadlessAssemblyRegion) ((ShardBoundaryShard<GATKRead>) shard).getShardBoundary();
        int extension = Math.max(shard.getInterval().getStart() - shard.getPaddedInterval().getStart(), shard.getPaddedInterval().getEnd() - shard.getInterval().getEnd());
        AssemblyRegion assemblyRegion = new AssemblyRegion(shard.getInterval(), Collections.emptyList(), readlessAssemblyRegion.isActive(), extension, header);
        assemblyRegion.addAll(Lists.newArrayList(shard));
        return assemblyRegion;
    }

    private static FlatMapFunction<Iterator<AssemblyRegion>, AssemblyRegionWalkerContext> getAssemblyRegionWalkerContextFunction(
            final String referenceFileName,
            final Broadcast<FeatureManager> bFeatureManager) {

        return (FlatMapFunction<Iterator<AssemblyRegion>, AssemblyRegionWalkerContext>) assemblyRegionIter -> {
            ReferenceDataSource reference = referenceFileName == null ? null : new ReferenceFileSource(IOUtils.getPath(SparkFiles.get(referenceFileName)));
            final FeatureManager features = bFeatureManager == null ? null : bFeatureManager.getValue();
            return Utils.stream(assemblyRegionIter).map(assemblyRegion ->
                    new AssemblyRegionWalkerContext(assemblyRegion,
                            new ReferenceContext(reference, assemblyRegion.getExtendedSpan()),
                            new FeatureContext(features, assemblyRegion.getExtendedSpan()))).iterator();
        };
    }

    @Override
    protected void runTool(JavaSparkContext ctx) {
        referenceFileName = addReferenceFilesForSpark(ctx, referenceArguments.getReferenceFileName());
        processAssemblyRegions(getAssemblyRegions(ctx), ctx);
    }

    /**
     * Process the assembly regions and write output. Must be implemented by subclasses.
     *
     * @param rdd a distributed collection of {@link AssemblyRegionWalkerContext}
     * @param ctx our Spark context
     */
    protected abstract void processAssemblyRegions(JavaRDD<AssemblyRegionWalkerContext> rdd, JavaSparkContext ctx);

}
