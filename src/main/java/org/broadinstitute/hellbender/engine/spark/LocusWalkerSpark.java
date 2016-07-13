package org.broadinstitute.hellbender.engine.spark;

import com.google.common.base.Function;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.hellbender.utils.iterators.IntervalOverlappingIterator;
import org.broadinstitute.hellbender.utils.locusiterator.LIBSDownsamplingInfo;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple3;

import javax.annotation.Nullable;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public abstract class LocusWalkerSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = "maxDepthPerSample", shortName = "maxDepthPerSample", doc = "Maximum number of reads to retain per sample per locus. Reads above this threshold will be downsampled. Set to 0 to disable.", optional = true)
    protected int maxDepthPerSample = defaultMaxDepthPerSample();

    /**
     * Returns default value for the {@link #maxDepthPerSample} parameter, if none is provided on the command line.
     * Default implementation returns 0 (no downsampling by default).
     */
    protected int defaultMaxDepthPerSample() {
        return 0;
    }

    @Argument(fullName="readShardSize", shortName="readShardSize", doc = "Maximum size of each read shard, in bases.", optional = true)
    public int readShardSize = 10000;

    @Argument(fullName="readShardPadding", shortName="readShardPadding", doc = "Each read shard has this many bases of extra context on each side.", optional = true)
    public int readShardPadding = 1000;

    @Argument(doc = "whether to use the shuffle implementation or overlaps partitioning (the default)", shortName = "shuffle", fullName = "shuffle", optional = true)
    public boolean shuffle = false;

    private FeatureManager features; // TODO: move up to GATKSparkTool?

    @Override
    protected void runPipeline(JavaSparkContext sparkContext) {
        initializeFeatures();
        super.runPipeline(sparkContext);
    }

    void initializeFeatures() {
        features = new FeatureManager(this);
        if ( features.isEmpty() ) {  // No available sources of Features discovered for this tool
            features = null;
        }
    }

    /** Returns the downsampling info using {@link #maxDepthPerSample} as target coverage. */
    protected final LIBSDownsamplingInfo getDownsamplingInfo() {
        if (maxDepthPerSample < 0) {
            throw new UserException.BadArgumentValue("maxDepthPerSample",
                    String.valueOf(maxDepthPerSample),
                    "should be a positive number");
        }
        return (maxDepthPerSample == 0) ? LocusIteratorByState.NO_DOWNSAMPLING : new LIBSDownsamplingInfo(true, maxDepthPerSample);
    }

    /**
     * Loads alignments and the corresponding reference and features into a {@link JavaRDD} for the intervals specified.
     *
     * If no intervals were specified, returns all the alignments.
     *
     * @return all alignments from as a {@link JavaRDD}, bounded by intervals if specified.
     */
    public JavaRDD<Tuple3<AlignmentContext, ReferenceContext, FeatureContext>> getAlignments(JavaSparkContext ctx) {
        SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        List<SimpleInterval> intervals = hasIntervals() ? getIntervals() : IntervalUtils.getAllIntervalsForReference(sequenceDictionary);
        final List<ShardBoundary> intervalShards = intervals.stream()
                .flatMap(interval -> Shard.divideIntervalIntoShards(interval, readShardSize, readShardPadding, sequenceDictionary).stream())
                .collect(Collectors.toList());
        int maxLocatableSize = Math.min(readShardSize, readShardPadding);
        JavaRDD<Shard<GATKRead>> shardedReads = SparkSharder.shard(ctx, getReads(), GATKRead.class, sequenceDictionary, intervalShards, maxLocatableSize, shuffle);
        Broadcast<ReferenceMultiSource> bReferenceSource = hasReference() ? ctx.broadcast(getReference()) : null;
        Broadcast<FeatureManager> bFeatureManager = features == null ? null : ctx.broadcast(features);
        return shardedReads.flatMap(getAlignmentsFunction(bReferenceSource, bFeatureManager, sequenceDictionary, getHeaderForReads(), getDownsamplingInfo()));
    }

    /**
     * Return a function that maps a {@link Shard} of reads into a tuple of alignments and their corresponding reference and features.
     * @param bReferenceSource the reference source broadcast
     * @param bFeatureManager the feature manager broadcast
     * @param sequenceDictionary the sequence dictionary for the reads
     * @param header the reads header
     * @param downsamplingInfo the downsampling method for the reads
     * @return a function that maps a {@link Shard} of reads into a tuple of alignments and their corresponding reference and features.
     */
    private static FlatMapFunction<Shard<GATKRead>, Tuple3<AlignmentContext, ReferenceContext, FeatureContext>> getAlignmentsFunction(
            Broadcast<ReferenceMultiSource> bReferenceSource, Broadcast<FeatureManager> bFeatureManager,
            SAMSequenceDictionary sequenceDictionary, SAMFileHeader header, LIBSDownsamplingInfo downsamplingInfo) {
        return (FlatMapFunction<Shard<GATKRead>, Tuple3<AlignmentContext, ReferenceContext, FeatureContext>>) shardedRead -> {
            SimpleInterval interval = shardedRead.getInterval();
            SimpleInterval paddedInterval = shardedRead.getPaddedInterval();
            Iterator<GATKRead> readIterator = shardedRead.iterator();
            ReferenceDataSource reference = bReferenceSource == null ? null :
                    new ReferenceMemorySource(bReferenceSource.getValue().getReferenceBases(null, paddedInterval), sequenceDictionary);
            FeatureManager fm = bFeatureManager == null ? null : bFeatureManager.getValue();

            final Set<String> samples = header.getReadGroups().stream()
                    .map(SAMReadGroupRecord::getSample)
                    .collect(Collectors.toSet());
            LocusIteratorByState libs = new LocusIteratorByState(readIterator, downsamplingInfo, false, samples, header, true, false);
            IntervalOverlappingIterator<AlignmentContext> alignmentContexts = new IntervalOverlappingIterator<>(libs, ImmutableList.of(interval), sequenceDictionary);
            Iterator<Tuple3<AlignmentContext, ReferenceContext, FeatureContext>> transform = Iterators.transform(alignmentContexts, new Function<AlignmentContext, Tuple3<AlignmentContext, ReferenceContext, FeatureContext>>() {
                @Nullable
                @Override
                public Tuple3<AlignmentContext, ReferenceContext, FeatureContext> apply(@Nullable AlignmentContext alignmentContext) {
                    final SimpleInterval alignmentInterval = new SimpleInterval(alignmentContext);
                    return new Tuple3<>(alignmentContext, new ReferenceContext(reference, alignmentInterval), new FeatureContext(fm, alignmentInterval));
                }
            });
            return transform;
        };
    }
}
