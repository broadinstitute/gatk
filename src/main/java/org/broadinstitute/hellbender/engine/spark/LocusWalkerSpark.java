package org.broadinstitute.hellbender.engine.spark;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.locusiterator.AlignmentContextIteratorBuilder;
import org.broadinstitute.hellbender.utils.locusiterator.LIBSDownsamplingInfo;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * A Spark version of {@link LocusWalker}. Subclasses should implement {@link #processAlignments(JavaRDD, JavaSparkContext)}
 * and operate on the passed in RDD.
 */
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

    @Override
    public boolean requiresReads() {
        return true;
    }

    /** Returns the downsampling info using {@link #maxDepthPerSample} as target coverage. */
    protected final LIBSDownsamplingInfo getDownsamplingInfo() {
        if (maxDepthPerSample < 0) {
            throw new CommandLineException.BadArgumentValue("maxDepthPerSample",
                    String.valueOf(maxDepthPerSample),
                    "should be a positive number");
        }
        return (maxDepthPerSample == 0) ? LocusIteratorByState.NO_DOWNSAMPLING : new LIBSDownsamplingInfo(true, maxDepthPerSample);
    }

    /**
     * Does this tool emit information for uncovered loci? Tools that do should override to return {@code true}.
     *
     * NOTE:  Typically, this should only be used when intervals are specified.
     * NOTE:  If MappedReadFilter is removed, then emitting empty loci will fail.
     * NOTE:  If there is no available sequence dictionary and this is set to true, there should be a failure.  Please
     *  consider requiring reads and/or references for all tools that wish to set this to {@code true}.
     *
     * @return {@code true} if this tool requires uncovered loci information to be emitted, {@code false} otherwise
     */
    public boolean emitEmptyLoci() {
        return false;
    }

    /**
     * Loads alignments and the corresponding reference and features into a {@link JavaRDD} for the intervals specified.
     *
     * If no intervals were specified, returns all the alignments.
     *
     * @return all alignments as a {@link JavaRDD}, bounded by intervals if specified.
     */
    public JavaRDD<LocusWalkerContext> getAlignments(JavaSparkContext ctx) {
        SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        List<SimpleInterval> intervals = hasUserSuppliedIntervals() ? getIntervals() : IntervalUtils.getAllIntervalsForReference(sequenceDictionary);
        final List<ShardBoundary> intervalShards = intervals.stream()
                .flatMap(interval -> Shard.divideIntervalIntoShards(interval, readShardSize, readShardPadding, sequenceDictionary).stream())
                .collect(Collectors.toList());
        int maxLocatableSize = Math.min(readShardSize, readShardPadding);
        JavaRDD<Shard<GATKRead>> shardedReads = SparkSharder.shard(ctx, getReads(), GATKRead.class, sequenceDictionary, intervalShards, maxLocatableSize, shuffle);
        Broadcast<ReferenceMultiSparkSource> bReferenceSource = hasReference() ? ctx.broadcast(getReference()) : null;
        Broadcast<FeatureManager> bFeatureManager = features == null ? null : ctx.broadcast(features);
        return shardedReads.flatMap(getAlignmentsFunction(bReferenceSource, bFeatureManager, sequenceDictionary, getHeaderForReads(), getDownsamplingInfo(), emitEmptyLoci()));
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
    private static FlatMapFunction<Shard<GATKRead>, LocusWalkerContext> getAlignmentsFunction(
            Broadcast<ReferenceMultiSparkSource> bReferenceSource, Broadcast<FeatureManager> bFeatureManager,
            SAMSequenceDictionary sequenceDictionary, SAMFileHeader header, LIBSDownsamplingInfo downsamplingInfo, boolean isEmitEmptyLoci) {
        return (FlatMapFunction<Shard<GATKRead>, LocusWalkerContext>) shardedRead -> {
            SimpleInterval interval = shardedRead.getInterval();
            SimpleInterval paddedInterval = shardedRead.getPaddedInterval();
            Iterator<GATKRead> readIterator = shardedRead.iterator();
            ReferenceDataSource reference = bReferenceSource == null ? null :
                    new ReferenceMemorySource(bReferenceSource.getValue().getReferenceBases(paddedInterval), sequenceDictionary);
            FeatureManager fm = bFeatureManager == null ? null : bFeatureManager.getValue();

            final SAMSequenceDictionary referenceDictionary = reference == null? null : reference.getSequenceDictionary();

            final AlignmentContextIteratorBuilder alignmentContextIteratorBuilder = new AlignmentContextIteratorBuilder();
            alignmentContextIteratorBuilder.setDownsamplingInfo(downsamplingInfo);
            alignmentContextIteratorBuilder.setEmitEmptyLoci(isEmitEmptyLoci);
            alignmentContextIteratorBuilder.setIncludeDeletions(true);
            alignmentContextIteratorBuilder.setKeepUniqueReadListInLibs(false);
            alignmentContextIteratorBuilder.setIncludeNs(false);

            final Iterator<AlignmentContext> alignmentContextIterator = alignmentContextIteratorBuilder.build(
                    readIterator, header, Collections.singletonList(interval), sequenceDictionary,
                    reference != null);

            return StreamSupport.stream(Spliterators.spliteratorUnknownSize(alignmentContextIterator, 0), false).map(alignmentContext -> {
                final SimpleInterval alignmentInterval = new SimpleInterval(alignmentContext);
                return new LocusWalkerContext(alignmentContext, new ReferenceContext(reference, alignmentInterval), new FeatureContext(fm, alignmentInterval));
            }).iterator();
        };
    }

    @Override
    protected void runTool(JavaSparkContext ctx) {
        processAlignments(getAlignments(ctx), ctx);
    }

    /**
     * Process the alignments and write output. Must be implemented by subclasses.
     *
     * @param rdd a distributed collection of {@link LocusWalkerContext}
     * @param ctx our Spark context
     */
    protected abstract void processAlignments(JavaRDD<LocusWalkerContext> rdd, JavaSparkContext ctx);
}
