package org.broadinstitute.hellbender.engine.spark;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * A Spark version of {@link ReadWalker}. Subclasses should implement {@link #processReads(JavaRDD, JavaSparkContext)}
 * and operate on the passed in RDD.
 */
public abstract class ReadWalkerSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() {
        return true;
    }

    /**
     * This number controls the size of the cache for our FeatureInputs
     * (specifically, the number of additional bases worth of overlapping records to cache when querying feature sources).
     */
    public static final int FEATURE_CACHE_LOOKAHEAD = 1_000;

    @Argument(fullName="readShardSize", shortName="readShardSize", doc = "Maximum size of each read shard, in bases.", optional = true)
    public int readShardSize = 10000;

    @Argument(fullName="readShardPadding", shortName="readShardPadding", doc = "Each read shard has this many bases of extra context on each side.", optional = true)
    public int readShardPadding = 1000;

    @Argument(doc = "whether to use the shuffle implementation or not", shortName = "shuffle", fullName = "shuffle", optional = true)
    public boolean shuffle = false;

    void initializeFeatures() {
        features = new FeatureManager(this, FEATURE_CACHE_LOOKAHEAD);
        if ( features.isEmpty() ) {  // No available sources of Features discovered for this tool
            features = null;
        }
    }

    /**
     * Loads reads and the corresponding reference and features into a {@link JavaRDD} for the intervals specified.
     *
     * If no intervals were specified, returns all the reads.
     *
     * @return all reads as a {@link JavaRDD}, bounded by intervals if specified.
     */
    public JavaRDD<ReadWalkerContext> getReads(JavaSparkContext ctx) {
        SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        List<SimpleInterval> intervals = hasUserSuppliedIntervals() ? getIntervals() : IntervalUtils.getAllIntervalsForReference(sequenceDictionary);
        // use unpadded shards (padding is only needed for reference bases)
        final List<ShardBoundary> intervalShards = intervals.stream()
                .flatMap(interval -> Shard.divideIntervalIntoShards(interval, readShardSize, 0, sequenceDictionary).stream())
                .collect(Collectors.toList());
        JavaRDD<Shard<GATKRead>> shardedReads = SparkSharder.shard(ctx, getReads(), GATKRead.class, sequenceDictionary, intervalShards, readShardSize, shuffle);
        Broadcast<ReferenceMultiSparkSource> bReferenceSource = hasReference() ? ctx.broadcast(getReference()) : null;
        Broadcast<FeatureManager> bFeatureManager = features == null ? null : ctx.broadcast(features);
        return shardedReads.flatMap(getReadsFunction(bReferenceSource, bFeatureManager, sequenceDictionary, readShardPadding));
    }

    private static FlatMapFunction<Shard<GATKRead>, ReadWalkerContext> getReadsFunction(
            Broadcast<ReferenceMultiSparkSource> bReferenceSource, Broadcast<FeatureManager> bFeatureManager,
            SAMSequenceDictionary sequenceDictionary, int readShardPadding) {
        return (FlatMapFunction<Shard<GATKRead>, ReadWalkerContext>) shard -> {
            // get reference bases for this shard (padded)
            SimpleInterval paddedInterval = shard.getInterval().expandWithinContig(readShardPadding, sequenceDictionary);
            ReferenceDataSource reference = bReferenceSource == null ? null :
                    new ReferenceMemorySource(bReferenceSource.getValue().getReferenceBases(paddedInterval), sequenceDictionary);
            FeatureManager features = bFeatureManager == null ? null : bFeatureManager.getValue();

            return StreamSupport.stream(shard.spliterator(), false)
                    .map(r -> {
                        final SimpleInterval readInterval = getReadInterval(r);
                        return new ReadWalkerContext(r, new ReferenceContext(reference, readInterval), new FeatureContext(features, readInterval));
                    }).iterator();
        };
    }

    /**
     * Returns an interval for the read.
     * Note: some walkers must be able to work on any read, including those whose coordinates do not form a valid SimpleInterval.
     * So here we check this condition and create null intervals for such reads.
     */
    static SimpleInterval getReadInterval(final GATKRead read) {
        return !read.isUnmapped() && SimpleInterval.isValid(read.getContig(), read.getStart(), read.getEnd()) ? new SimpleInterval(read) : null;
    }

    @Override
    protected void runTool(JavaSparkContext ctx) {
        processReads(getReads(ctx), ctx);
    }

    /**
     * Process the reads and write output. Must be implemented by subclasses.
     *
     * @param rdd a distributed collection of {@link ReadWalkerContext}
     * @param ctx our Spark context
     */
    protected abstract void processReads(JavaRDD<ReadWalkerContext> rdd, JavaSparkContext ctx);
}
