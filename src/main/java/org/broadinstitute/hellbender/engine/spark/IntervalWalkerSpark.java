package org.broadinstitute.hellbender.engine.spark;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * A Spark version of {@link IntervalWalker}. Subclasses should implement {@link #processIntervals(JavaRDD, JavaSparkContext)}
 * and operate on the passed in RDD.
 */
public abstract class IntervalWalkerSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresIntervals() {
        return true;
    }

    @Argument(doc = "whether to use the shuffle implementation or not", shortName = "shuffle", fullName = "shuffle", optional = true)
    public boolean shuffle = false;

    @Argument(fullName="intervalShardPadding", shortName="intervalShardPadding", doc = "Each interval shard has this many bases of extra context on each side.", optional = true)
    public int intervalShardPadding = 1000;

    /**
     * Customize initialization of the Feature data source for this traversal type to disable query lookahead.
     */
    @Override
    void initializeFeatures() {
        // Disable query lookahead in our FeatureManager for this traversal type. Query lookahead helps
        // when our query intervals are overlapping and gradually increasing in position (as they are
        // with ReadWalkers, typically), but with IntervalWalkers our query intervals are guaranteed
        // to be non-overlapping, since our interval parsing code always merges overlapping intervals.
        features = new FeatureManager(this, 0);
        if ( features.isEmpty() ) {  // No available sources of Features for this tool
            features = null;
        }
    }

    /**
     * Loads intervals and the corresponding reads, reference and features into a {@link JavaRDD}.
     *
     * @return all intervals as a {@link JavaRDD}.
     */
    public JavaRDD<IntervalWalkerContext> getIntervals(JavaSparkContext ctx) {
        SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        // don't shard the intervals themselves, since we want each interval to be processed by a single task
        final List<ShardBoundary> intervalShardBoundaries = getIntervals().stream()
                .map(i -> new ShardBoundary(i, i)).collect(Collectors.toList());
        JavaRDD<Shard<GATKRead>> shardedReads = SparkSharder.shard(ctx, getReads(), GATKRead.class, sequenceDictionary, intervalShardBoundaries, Integer.MAX_VALUE, shuffle);
        Broadcast<ReferenceMultiSparkSource> bReferenceSource = hasReference() ? ctx.broadcast(getReference()) : null;
        Broadcast<FeatureManager> bFeatureManager = features == null ? null : ctx.broadcast(features);
        return shardedReads.map(getIntervalsFunction(bReferenceSource, bFeatureManager, sequenceDictionary, intervalShardPadding));
    }

    private static org.apache.spark.api.java.function.Function<Shard<GATKRead>, IntervalWalkerContext> getIntervalsFunction(
            Broadcast<ReferenceMultiSparkSource> bReferenceSource, Broadcast<FeatureManager> bFeatureManager,
            SAMSequenceDictionary sequenceDictionary, int intervalShardPadding) {
        return (org.apache.spark.api.java.function.Function<Shard<GATKRead>, IntervalWalkerContext>) shard -> {
            // get reference bases for this shard (padded)
            SimpleInterval interval = shard.getInterval();
            SimpleInterval paddedInterval = shard.getInterval().expandWithinContig(intervalShardPadding, sequenceDictionary);
            ReadsContext readsContext = new ReadsContext(new GATKDataSource<GATKRead>() {
                @Override
                public Iterator<GATKRead> iterator() {
                    return shard.iterator();
                }
                @Override
                public Iterator<GATKRead> query(SimpleInterval interval) {
                    return StreamSupport.stream(shard.spliterator(), false).filter(
                            r -> IntervalUtils.overlaps(r, interval)).iterator();
                }
            }, shard.getInterval());
            ReferenceDataSource reference = bReferenceSource == null ? null :
                    new ReferenceMemorySource(bReferenceSource.getValue().getReferenceBases(paddedInterval), sequenceDictionary);
            FeatureManager features = bFeatureManager == null ? null : bFeatureManager.getValue();
            return new IntervalWalkerContext(interval, readsContext, new ReferenceContext(reference, interval), new FeatureContext(features, interval));
        };
    }

    @Override
    protected void runTool(JavaSparkContext ctx) {
        processIntervals(getIntervals(ctx), ctx);
    }

    /**
     * Process the intervals and write output. Must be implemented by subclasses.
     *
     * @param rdd a distributed collection of {@link IntervalWalkerContext}
     * @param ctx our Spark context
     */
    protected abstract void processIntervals(JavaRDD<IntervalWalkerContext> rdd, JavaSparkContext ctx);
}
