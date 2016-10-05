package org.broadinstitute.hellbender.engine.spark;

import com.google.common.base.Function;
import com.google.common.collect.Iterators;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple3;

import javax.annotation.Nullable;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A Spark version of {@link ReadWalker}.
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
    public JavaRDD<Tuple3<GATKRead, ReferenceContext, FeatureContext>> getReads(JavaSparkContext ctx) {
        SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        List<SimpleInterval> intervals = hasIntervals() ? getIntervals() : IntervalUtils.getAllIntervalsForReference(sequenceDictionary);
        // use unpadded shards (padding is only needed for reference bases)
        final List<ShardBoundary> intervalShards = intervals.stream()
                .flatMap(interval -> Shard.divideIntervalIntoShards(interval, readShardSize, 0, sequenceDictionary).stream())
                .collect(Collectors.toList());
        JavaRDD<Shard<GATKRead>> shardedReads = SparkSharder.shard(ctx, getReads(), GATKRead.class, sequenceDictionary, intervalShards, readShardSize, shuffle);
        Broadcast<ReferenceMultiSource> bReferenceSource = hasReference() ? ctx.broadcast(getReference()) : null;
        Broadcast<FeatureManager> bFeatureManager = features == null ? null : ctx.broadcast(features);
        return shardedReads.flatMap(getReadsFunction(bReferenceSource, bFeatureManager, sequenceDictionary, readShardPadding));
    }

    private static FlatMapFunction<Shard<GATKRead>, Tuple3<GATKRead, ReferenceContext, FeatureContext>> getReadsFunction(
            Broadcast<ReferenceMultiSource> bReferenceSource, Broadcast<FeatureManager> bFeatureManager,
            SAMSequenceDictionary sequenceDictionary, int readShardPadding) {
        return (FlatMapFunction<Shard<GATKRead>, Tuple3<GATKRead, ReferenceContext, FeatureContext>>) shard -> {
            // get reference bases for this shard (padded)
            SimpleInterval paddedInterval = shard.getInterval().expandWithinContig(readShardPadding, sequenceDictionary);
            ReferenceDataSource reference = bReferenceSource == null ? null :
                    new ReferenceMemorySource(bReferenceSource.getValue().getReferenceBases(null, paddedInterval), sequenceDictionary);
            FeatureManager features = bFeatureManager == null ? null : bFeatureManager.getValue();

            Iterator<Tuple3<GATKRead, ReferenceContext, FeatureContext>> transform = Iterators.transform(shard.iterator(), new Function<GATKRead, Tuple3<GATKRead, ReferenceContext, FeatureContext>>() {
                @Nullable
                @Override
                public Tuple3<GATKRead, ReferenceContext, FeatureContext> apply(@Nullable GATKRead read) {
                    final SimpleInterval readInterval = getReadInterval(read);
                    return new Tuple3<>(read, new ReferenceContext(reference, readInterval), new FeatureContext(features, readInterval));
                }
            });
            // only include reads that start in the shard
            return Iterators.filter(transform, r -> r._1().getStart() >= shard.getStart()
                    && r._1().getStart() <= shard.getEnd());
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
}
