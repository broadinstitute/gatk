package org.broadinstitute.hellbender.engine.spark;

import com.google.common.base.Function;
import com.google.common.collect.Iterators;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Advanced;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple3;

import javax.annotation.Nullable;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A Spark version of {@link AssemblyRegionWalker}.
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

    @Argument(doc = "whether to use the shuffle implementation or not", shortName = "shuffle", fullName = "shuffle", optional = true)
    public boolean shuffle = false;

    @Override
    public final boolean requiresReads() { return true; }

    @Override
    public final boolean requiresReference() { return true; }

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
    public JavaRDD<Tuple3<AssemblyRegion, ReferenceContext, FeatureContext>> getAssemblyRegions(JavaSparkContext ctx) {
        SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        JavaRDD<Shard<GATKRead>> shardedReads = SparkSharder.shard(ctx, getReads(), GATKRead.class, sequenceDictionary, intervalShards, readShardSize, shuffle);
        Broadcast<ReferenceMultiSource> bReferenceSource = hasReference() ? ctx.broadcast(getReference()) : null;
        Broadcast<FeatureManager> bFeatureManager = features == null ? null : ctx.broadcast(features);
        return shardedReads.flatMap(getAssemblyRegionsFunction(bReferenceSource, bFeatureManager, sequenceDictionary, getHeaderForReads(),
                assemblyRegionEvaluator(), minAssemblyRegionSize, maxAssemblyRegionSize, assemblyRegionPadding, activeProbThreshold, maxProbPropagationDistance));
    }

    private static FlatMapFunction<Shard<GATKRead>, Tuple3<AssemblyRegion, ReferenceContext, FeatureContext>> getAssemblyRegionsFunction(
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
        return (FlatMapFunction<Shard<GATKRead>, Tuple3<AssemblyRegion, ReferenceContext, FeatureContext>>) shardedRead -> {
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
            return Iterators.transform(assemblyRegions.iterator(), new Function<AssemblyRegion, Tuple3<AssemblyRegion, ReferenceContext, FeatureContext>>() {
                @Nullable
                @Override
                public Tuple3<AssemblyRegion, ReferenceContext, FeatureContext> apply(@Nullable AssemblyRegion assemblyRegion) {
                    return new Tuple3<>(assemblyRegion,
                            new ReferenceContext(reference, assemblyRegion.getExtendedSpan()),
                            new FeatureContext(features, assemblyRegion.getExtendedSpan()));
                }
            });
        };
    }

}
