package org.broadinstitute.hellbender.engine.spark;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.AssemblyRegionWalkerShardingArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ShardingArgumentCollection;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * A Spark version of {@link AssemblyRegionWalker}. Subclasses should implement {@link #processAssemblyRegions(JavaRDD, JavaSparkContext)}
 * and operate on the passed in RDD.
 */
public abstract class AssemblyRegionWalkerSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    protected ShardingArgumentCollection shardingArgumentCollection = defaultShardingArgumentCollection();

    @Argument(fullName = "minAssemblyRegionSize", shortName = "minAssemblyRegionSize", doc = "Minimum size of an assembly region", optional = true)
    protected int minAssemblyRegionSize = defaultMinAssemblyRegionSize();

    @Argument(fullName = "maxAssemblyRegionSize", shortName = "maxAssemblyRegionSize", doc = "Maximum size of an assembly region", optional = true)
    protected int maxAssemblyRegionSize = defaultMaxAssemblyRegionSize();

    @Argument(fullName = "assemblyRegionPadding", shortName = "assemblyRegionPadding", doc = "Number of additional bases of context to include around each assembly region", optional = true)
    protected int assemblyRegionPadding = defaultAssemblyRegionPadding();

    @Argument(fullName = "maxReadsPerAlignmentStart", shortName = "maxReadsPerAlignmentStart", doc = "Maximum number of reads to retain per alignment start position. Reads above this threshold will be downsampled. Set to 0 to disable.", optional = true)
    protected int maxReadsPerAlignmentStart = defaultMaxReadsPerAlignmentStart();

    @Advanced
    @Argument(fullName = "activeProbabilityThreshold", shortName = "activeProbabilityThreshold", doc="Minimum probability for a locus to be considered active.", optional = true)
    protected double activeProbThreshold = defaultActiveProbThreshold();

    @Advanced
    @Argument(fullName = "maxProbPropagationDistance", shortName = "maxProbPropagationDistance", doc="Upper limit on how many bases away probability mass can be moved around when calculating the boundaries between active and inactive assembly regions", optional = true)
    protected int maxProbPropagationDistance = defaultMaxProbPropagationDistance();

    /**
     * Default implementation returns an {@link AssemblyRegionWalkerShardingArgumentCollection} with {@link #defaultReadShardPadding()} and {@link #defaultReadShardPadding()} default values.
     * @return
     */
    public ShardingArgumentCollection defaultShardingArgumentCollection() {
        return new AssemblyRegionWalkerShardingArgumentCollection(defaultReadShardSize(), defaultReadShardPadding());
    }

    /**
     * @return Default value for the {@link AssemblyRegionWalkerShardingArgumentCollection#readShardSize} parameter, if none is provided on the command line
     */
    protected abstract int defaultReadShardSize();

    /**
     * @return Default value for the {@link AssemblyRegionWalkerShardingArgumentCollection#readShardPadding} parameter, if none is provided on the command line
     */
    protected abstract int defaultReadShardPadding();

    /**
     * @return Default value for the {@link #minAssemblyRegionSize} parameter, if none is provided on the command line
     */
    protected abstract int defaultMinAssemblyRegionSize();

    /**
     * @return Default value for the {@link #maxAssemblyRegionSize} parameter, if none is provided on the command line
     */
    protected abstract int defaultMaxAssemblyRegionSize();

    /**
     * @return Default value for the {@link #assemblyRegionPadding} parameter, if none is provided on the command line
     */
    protected abstract int defaultAssemblyRegionPadding();

    /**
     * @return Default value for the {@link #maxReadsPerAlignmentStart} parameter, if none is provided on the command line
     */
    protected abstract int defaultMaxReadsPerAlignmentStart();

    /**
     * @return Default value for the {@link #activeProbThreshold} parameter, if none is provided on the command line
     */
    protected abstract double defaultActiveProbThreshold();

    /**
     * @return Default value for the {@link #maxProbPropagationDistance} parameter, if none is provided on the command line
     */
    protected abstract int defaultMaxProbPropagationDistance();

    @Argument(doc = "whether to use the shuffle implementation or not", shortName = "shuffle", fullName = "shuffle", optional = true)
    public boolean shuffle = false;

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

    private List<ShardBoundary> intervalShards;

    /**
     * Note that this sets {@code intervalShards} as a side effect, in order to add padding to the intervals.
     */
    @Override
    protected List<SimpleInterval> editIntervals(List<SimpleInterval> rawIntervals) {
        SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        List<SimpleInterval> intervals = rawIntervals == null ? IntervalUtils.getAllIntervalsForReference(sequenceDictionary) : rawIntervals;
        intervalShards = intervals.stream()
                .flatMap(interval -> Shard.divideIntervalIntoShards(interval, shardingArgumentCollection.getShardSize(), shardingArgumentCollection.getShardStep(), shardingArgumentCollection.getShardPadding(), sequenceDictionary).stream())
                .collect(Collectors.toList());
        List<SimpleInterval> paddedIntervalsForReads =
                intervals.stream().map(interval -> interval.expandWithinContig(shardingArgumentCollection.getShardPadding(), sequenceDictionary)).collect(Collectors.toList());
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
        JavaRDD<Shard<GATKRead>> shardedReads = SparkSharder.shard(ctx, getReads(), GATKRead.class, sequenceDictionary, intervalShards, shardingArgumentCollection.getShardSize(), shuffle);
        Broadcast<ReferenceMultiSource> bReferenceSource = hasReference() ? ctx.broadcast(getReference()) : null;
        Broadcast<FeatureManager> bFeatureManager = features == null ? null : ctx.broadcast(features);
        return shardedReads.flatMap(getAssemblyRegionsFunction(bReferenceSource, bFeatureManager, sequenceDictionary, getHeaderForReads(),
                assemblyRegionEvaluator(), minAssemblyRegionSize, maxAssemblyRegionSize, assemblyRegionPadding, activeProbThreshold, maxProbPropagationDistance));
    }

    private static FlatMapFunction<Shard<GATKRead>, AssemblyRegionWalkerContext> getAssemblyRegionsFunction(
            final Broadcast<ReferenceMultiSource> bReferenceSource,
            final Broadcast<FeatureManager> bFeatureManager,
            final SAMSequenceDictionary sequenceDictionary,
            final SAMFileHeader header,
            final AssemblyRegionEvaluator evaluator,
            final int minAssemblyRegionSize,
            final int maxAssemblyRegionSize,
            final int assemblyRegionPadding,
            final double activeProbThreshold,
            final int maxProbPropagationDistance) {
        return (FlatMapFunction<Shard<GATKRead>, AssemblyRegionWalkerContext>) shardedRead -> {
            SimpleInterval paddedInterval = shardedRead.getPaddedInterval();
            SimpleInterval assemblyRegionPaddedInterval = paddedInterval.expandWithinContig(assemblyRegionPadding, sequenceDictionary);

            ReferenceDataSource reference = bReferenceSource == null ? null :
                    new ReferenceMemorySource(bReferenceSource.getValue().getReferenceBases(null, assemblyRegionPaddedInterval), sequenceDictionary);
            FeatureManager features = bFeatureManager == null ? null : bFeatureManager.getValue();
            ReferenceContext referenceContext = new ReferenceContext(reference, paddedInterval);
            FeatureContext featureContext = new FeatureContext(features, paddedInterval);

            final Iterable<AssemblyRegion> assemblyRegions = AssemblyRegion.createFromReadShard(shardedRead,
                    header, referenceContext, featureContext, evaluator,
                    minAssemblyRegionSize, maxAssemblyRegionSize, assemblyRegionPadding, activeProbThreshold,
                    maxProbPropagationDistance);
            return StreamSupport.stream(assemblyRegions.spliterator(), false).map(assemblyRegion ->
                    new AssemblyRegionWalkerContext(assemblyRegion,
                        new ReferenceContext(reference, assemblyRegion.getExtendedSpan()),
                        new FeatureContext(features, assemblyRegion.getExtendedSpan()))).iterator();
        };
    }

    @Override
    protected void runTool(JavaSparkContext ctx) {
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
