package org.broadinstitute.hellbender.engine;

import com.google.common.base.Function;
import com.google.common.collect.Iterators;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKCommandLinePluginDescriptor;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptor;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.engine.spark.SparkSharder;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple3;

import javax.annotation.Nullable;
import java.io.IOException;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * A ReadWalker is a tool that processes a single read at a time from one or multiple sources of reads, with
 * optional contextual information from a reference and/or sets of variants/Features.
 *
 * If multiple sources of reads are specified, they are merged together into a single sorted stream of reads.
 *
 * ReadWalker authors must implement the apply() method to process each read, and may optionally implement
 * onTraversalStart() and/or onTraversalSuccess(). See the PrintReadsWithReference walker for an example.
 */
public abstract class ReadWalker extends GATKTool {

    @Override
    public boolean requiresReads() {
        return true;
    }

    /**
     * This number controls the size of the cache for our FeatureInputs
     * (specifically, the number of additional bases worth of overlapping records to cache when querying feature sources).
     */
    public static final int FEATURE_CACHE_LOOKAHEAD = 1_000;

    /**
     * Return the list of GATKCommandLinePluginDescriptors to be used for this tool.
     * Uses the read filter plugin.
     */
    @Override
    protected List<? extends GATKCommandLinePluginDescriptor<?>> getPluginDescriptors() {
        return Collections.singletonList(new GATKReadFilterPluginDescriptor(getDefaultReadFilters()));
    }

    /**
     * Initialize data sources for traversal.
     *
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalStart() instead.
     */
    @Override
    protected final void onStartup() {
        super.onStartup();

        setReadTraversalBounds();
    }

    /**
     * Initialize traversal bounds if intervals are specified
     */
    void setReadTraversalBounds() {
        if ( hasIntervals() ) {
            reads.setTraversalBounds(intervalArgumentCollection.getTraversalParameters(getHeaderForReads().getSequenceDictionary()));
        }
    }

    @Override
    void initializeFeatures() {
        //We override this method to change lookahead of the cache
        features = new FeatureManager(this, FEATURE_CACHE_LOOKAHEAD);
        if ( features.isEmpty() ) {  // No available sources of Features discovered for this tool
            features = null;
        }
    }

    /**
     * Implementation of read-based traversal.
     * Subclasses can override to provide own behavior but default implementation should be suitable for most uses.
     *
     * The default implementation creates filters using {@link #makeReadFilter}
     * and then iterates over all reads, applies the filter and hands the resulting reads to the {@link #apply}
     * function of the walker (along with additional contextual information, if present, such as reference bases).
     */
    @Override
    public void traverse() {
        // Process each read in the input stream.
        // Supply reference bases spanning each read, if a reference is available.
        final CountingReadFilter countedFilter = makeReadFilter();

        StreamSupport.stream(reads.spliterator(), false)
                .filter(countedFilter)
                .forEach(read -> {
                    final SimpleInterval readInterval = getReadInterval(read);
                    apply(read,
                          new ReferenceContext(reference, readInterval), // Will create an empty ReferenceContext if reference or readInterval == null
                          new FeatureContext(features, readInterval));   // Will create an empty FeatureContext if features or readInterval == null

                    progressMeter.update(readInterval);
                });

        logger.info(countedFilter.getSummaryLine());
    }

    /**
     * Get a {@link Stream} of reads, for custom processing using the Java Streams API.
     *
     * Subclasses should override {@link #traverse()} and call this method to provide their own custom processing
     * of the stream.
     *
     * @return all reads from in a {@link Stream}, bounded by intervals if specified.
     */
    protected Stream<GATKRead> getReadsStream() {
        // Process each read in the input stream.
        // Supply reference bases spanning each read, if a reference is available.
        final CountingReadFilter countedFilter = makeReadFilter();

        return StreamSupport.stream(reads.spliterator(), false)
                .filter(countedFilter);
    }

    /**
     * Get a {@link Stream} of reads, (along with additional contextual information, if present, such as reference
     * bases), for custom processing using the Java Streams API.
     *
     * Subclasses should override {@link #traverse()} and call this method to provide their own custom processing
     * of the stream.
     *
     * @return all reads from in a {@link Stream}, bounded by intervals if specified.
     */
    protected Stream<Tuple3<GATKRead, ReferenceContext, FeatureContext>> getReadsWithContextStream() {
        // Process each read in the input stream.
        // Supply reference bases spanning each read, if a reference is available.
        final CountingReadFilter countedFilter = makeReadFilter();

        return StreamSupport.stream(reads.spliterator(), false)
                .filter(countedFilter)
                .map(read -> {
                    final SimpleInterval readInterval = getReadInterval(read);
                    progressMeter.update(readInterval);
                    return new Tuple3<>(read,
                            new ReferenceContext(reference, readInterval), // Will create an empty ReferenceContext if reference or readInterval == null
                            new FeatureContext(features, readInterval));   // Will create an empty FeatureContext if features or readInterval == null)
                });
    }

    /**
     * Returns an interval for the read.
     * Note: some walkers must be able to work on any read, including those whose coordinates do not form a valid SimpleInterval.
     * So here we check this condition and create null intervals for such reads.
     */
    static SimpleInterval getReadInterval(final GATKRead read) {
        return !read.isUnmapped() && SimpleInterval.isValid(read.getContig(), read.getStart(), read.getEnd()) ? new SimpleInterval(read) : null;
    }

    /**
     * Returns the read filter (simple or composite) that will be applied to the reads before calling {@link #apply}.
     * The default implementation combines the default read filters for this tool (returned by
     * {@link org.broadinstitute.hellbender.engine.ReadWalker#getDefaultReadFilters} with any read filter command
     * line arguments specified by the user; wraps each filter in the resulting list with a CountingReadFilter;
     * and returns a single composite filter resulting from the list by and'ing them together.
     *
     * Default tool implementation of {@link #traverse()} calls this method once before iterating
     * over the reads and reuses the filter object to avoid object allocation. Nevertheless, keeping state in filter
     * objects is strongly discouraged.
     *
     * Multiple filters can be composed by using {@link org.broadinstitute.hellbender.engine.filters.ReadFilter}
     * composition methods.
     */
    public CountingReadFilter makeReadFilter(){
        final GATKReadFilterPluginDescriptor readFilterPlugin =
                commandLineParser.getPluginDescriptor(GATKReadFilterPluginDescriptor.class);
        return readFilterPlugin.getMergedCountingReadFilter(getHeaderForReads());
    }

    /**
     * Returns the default list of CommandLineReadFilters that are used for this tool. The filters returned
     * by this method are subject to selective enabling/disabling by the user via the command line. The
     * default implementation uses the {@link WellformedReadFilter} filter with all default options. Subclasses
     * can override to provide alternative filters.
     *
     * Note: this method is called before command line parsing begins, and thus before a SAMFileHeader is
     * available through {link #getHeaderForReads}.
     *
     * @return List of individual filters to be applied for this tool.
     */
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(new WellformedReadFilter());
    }

    /**
     * Process an individual read (with optional contextual information). Must be implemented by tool authors.
     * In general, tool authors should simply stream their output from apply(), and maintain as little internal state
     * as possible.
     *
     * TODO: Determine whether and to what degree the GATK engine should provide a reduce operation
     * TODO: to complement this operation. At a minimum, we should make apply() return a value to
     * TODO: discourage statefulness in walkers, but how this value should be handled is TBD.
     * @param read current read
     * @param referenceContext Reference bases spanning the current read. Will be an empty, but non-null, context object
     *                         if there is no backing source of reference data (in which case all queries on it will return
     *                         an empty array/iterator). Can request extra bases of context around the current read's interval
     *                         by invoking {@link ReferenceContext#setWindow} on this object before calling {@link ReferenceContext#getBases}
     * @param featureContext Features spanning the current read. Will be an empty, but non-null, context object
     *                       if there is no backing source of Feature data (in which case all queries on it will return an
     *                       empty List).
     */
    public abstract void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext );

    /**
     * Shutdown data sources.
     *
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalSuccess() instead.
     */
    @Override
    protected final void onShutdown() {
        // Overridden only to make final so that concrete tool implementations don't override
        super.onShutdown();
    }

    @Argument(fullName="readShardSize", shortName="readShardSize", doc = "Maximum size of each read shard, in bases.", optional = true)
    public int readShardSize = 10000;

    @Argument(fullName="readShardPadding", shortName="readShardPadding", doc = "Each read shard has this many bases of extra context on each side.", optional = true)
    public int readShardPadding = 1000;

    @Argument(doc = "whether to use the shuffle implementation or not", shortName = "shuffle", fullName = "shuffle", optional = true)
    public boolean shuffle = false;

    /**
     * Loads reads into a {@link JavaRDD} for the intervals specified.
     * If no intervals were specified, returns all the reads.
     *
     * Subclasses should override {@link #traverse()} and call this method to provide their own custom processing
     * of the RDD.
     *
     * @return all reads from as a {@link JavaRDD}, bounded by intervals if specified.
     */
    protected JavaRDD<GATKRead> getReadsRdd() {
        final ReadFilter filter = makeReadFilter();
        return getUnfilteredReads().filter(read -> filter.test(read));
    }

    /**
     * Loads reads and the corresponding reference and features into a {@link JavaRDD} for the intervals specified.
     * If no intervals were specified, returns all the reads.
     *
     * Subclasses should override {@link #traverse()} and call this method to provide their own custom processing
     * of the RDD.
     *
     * @return all reads from as a {@link JavaRDD}, bounded by intervals if specified.
     */
    protected JavaRDD<Tuple3<GATKRead, ReferenceContext, FeatureContext>> getReadsWithContextRdd() {
        SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        List<SimpleInterval> intervals = hasIntervals() ? intervalsForTraversal : IntervalUtils.getAllIntervalsForReference(sequenceDictionary);
        // use unpadded shards (padding is only needed for reference bases)
        final List<ShardBoundary> intervalShards = intervals.stream()
                .flatMap(interval -> Shard.divideIntervalIntoShards(interval, readShardSize, 0, sequenceDictionary).stream())
                .collect(Collectors.toList());
        final ReadFilter filter = makeReadFilter();
        JavaRDD<GATKRead> reads = getUnfilteredReads().filter(read -> filter.test(read));
        JavaRDD<Shard<GATKRead>> shardedReads = SparkSharder.shard(ctx, reads, GATKRead.class, sequenceDictionary, intervalShards, readShardSize, shuffle);
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
            return () -> Iterators.filter(transform, r -> r._1().getStart() >= shard.getStart()
                    && r._1().getStart() <= shard.getEnd());
        };
    }

    /**
     * Loads the reads into a {@link JavaRDD} using the intervals specified, and returns them
     * without applying any filtering.
     *
     * If no intervals were specified, returns all the reads (both mapped and unmapped).
     *
     * @return all reads from our reads input(s) as a {@link JavaRDD}, bounded by intervals if specified, and unfiltered.
     */
    JavaRDD<GATKRead> getUnfilteredReads() {
        // TODO: This if statement is a temporary hack until #959 gets resolved.
        if (readInput.endsWith(".adam")) {
            try {
                return readsSource.getADAMReads(readInput, intervalsForTraversal, getHeaderForReads());
            } catch (IOException e) {
                throw new UserException("Failed to read ADAM file " + readInput, e);
            }

        } else {
            if (hasCramInput() && !hasReference()){
                throw new UserException.MissingReference("A reference file is required when using CRAM files.");
            }
            final String refPath = hasReference() ?  referenceArguments.getReferenceFile().getAbsolutePath() : null;
            // If no intervals were specified (intervals == null), this will return all reads (mapped and unmapped)
            return readsSource.getParallelReads(readInput, refPath, intervalsForTraversal, 0 /* TODO bamPartitionSplitSize */);
        }
    }

    /**
     * @return our reference source, or null if no reference is present
     */
    public ReferenceMultiSource getReference() {
        return referenceSource;
    }
}
