package org.broadinstitute.hellbender.engine.spark;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A Spark version of {@link AssemblyRegionWalker}. Subclasses should implement {@link #processAssemblyRegions(JavaRDD, JavaSparkContext)}
 * and operate on the passed in RDD.
 */
public abstract class AssemblyRegionWalkerSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(fullName="readShardSize", shortName="readShardSize", doc = "Maximum size of each read shard, in bases. For good performance, this should be much larger than the maximum assembly region size.", optional = true)
    protected int readShardSize = defaultReadShardSize();

    @Argument(fullName="readShardPadding", shortName="readShardPadding", doc = "Each read shard has this many bases of extra context on each side. Read shards must have as much or more padding than assembly regions.", optional = true)
    protected int readShardPadding = defaultReadShardPadding();

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
     * @return Default value for the {@link #readShardSize} parameter, if none is provided on the command line
     */
    protected abstract int defaultReadShardSize();

    /**
     * @return Default value for the {@link #readShardPadding} parameter, if none is provided on the command line
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

    /**
     * subclasses can override this to control if reads with deletions should be included in isActive pileups
     */
    protected abstract boolean includeReadsWithDeletionsInIsActivePileups();

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
                .flatMap(interval -> Shard.divideIntervalIntoShards(interval, readShardSize, readShardPadding, sequenceDictionary).stream())
                .collect(Collectors.toList());
        List<SimpleInterval> paddedIntervalsForReads =
                intervals.stream().map(interval -> interval.expandWithinContig(readShardPadding, sequenceDictionary)).collect(Collectors.toList());
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
        JavaRDD<Shard<GATKRead>> shardedReads = SparkSharder.shard(ctx, getReads(), GATKRead.class, sequenceDictionary, intervalShards, readShardSize, shuffle);
        Broadcast<ReferenceMultiSparkSource> bReferenceSource = hasReference() ? ctx.broadcast(getReference()) : null;
        Broadcast<FeatureManager> bFeatureManager = features == null ? null : ctx.broadcast(features);
        return shardedReads.flatMap(getAssemblyRegionsFunction(bReferenceSource, bFeatureManager, sequenceDictionary, getHeaderForReads(),
                assemblyRegionEvaluator(), minAssemblyRegionSize, maxAssemblyRegionSize, assemblyRegionPadding, activeProbThreshold, maxProbPropagationDistance, includeReadsWithDeletionsInIsActivePileups()));
    }

    private static FlatMapFunction<Shard<GATKRead>, AssemblyRegionWalkerContext> getAssemblyRegionsFunction(
            final Broadcast<ReferenceMultiSparkSource> bReferenceSource,
            final Broadcast<FeatureManager> bFeatureManager,
            final SAMSequenceDictionary sequenceDictionary,
            final SAMFileHeader header,
            final AssemblyRegionEvaluator evaluator,
            final int minAssemblyRegionSize,
            final int maxAssemblyRegionSize,
            final int assemblyRegionPadding,
            final double activeProbThreshold,
            final int maxProbPropagationDistance,
            final boolean includeReadsWithDeletionsInIsActivePileups) {
        return (FlatMapFunction<Shard<GATKRead>, AssemblyRegionWalkerContext>) shardedRead -> {
            final SimpleInterval paddedInterval = shardedRead.getPaddedInterval();
            final SimpleInterval assemblyRegionPaddedInterval = paddedInterval.expandWithinContig(assemblyRegionPadding, sequenceDictionary);

            final ReferenceDataSource reference = bReferenceSource == null ? null :
                    new ReferenceMemorySource(bReferenceSource.getValue().getReferenceBases(assemblyRegionPaddedInterval), sequenceDictionary);
            final FeatureManager features = bFeatureManager == null ? null : bFeatureManager.getValue();

            final Iterator<AssemblyRegion> assemblyRegionIter = new AssemblyRegionIterator(
                    new ShardToMultiIntervalShardAdapter<>(shardedRead),
                    header, reference, features, evaluator,
                    minAssemblyRegionSize, maxAssemblyRegionSize, assemblyRegionPadding, activeProbThreshold,
                    maxProbPropagationDistance, includeReadsWithDeletionsInIsActivePileups);
            final Iterable<AssemblyRegion> assemblyRegions = () -> assemblyRegionIter;
            return Utils.stream(assemblyRegions).map(assemblyRegion ->
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
