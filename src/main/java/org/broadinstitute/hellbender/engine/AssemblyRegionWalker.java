package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.cmdline.Advanced;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.ArrayList;
import java.util.List;

/**
 * An AssemblyRegionWalker is a tool that processes an entire region of reads at a time, each marked as either "active"
 * (containing possible variation) or "inactive" (not likely to contain actual variation). Tool authors must implement
 * {@link #assemblyRegionEvaluator} to provide a means of determining whether a site is active or not, as well as
 * {@link #apply} to process each region. Authors must also implement methods providing default values for the
 * various traversal parameters.
 *
 * Each region passed to {@link #apply} will come pre-marked as either "active" or "inactive" using the results of the
 * configured {@link #assemblyRegionEvaluator}. The {@link #assemblyRegionEvaluator} is used to evaluate the probability
 * that each individual locus is active, and these probabilities are used to determine the bounds of each active and
 * inactive region.
 *
 * {@link #apply} will be called once for each active AND inactive region, and it is up to the implementation how to
 * handle/process active vs. inactive regions.
 *
 * Internally, the reads are loaded in chunks called read shards, which are then subdivided into active/inactive regions
 * for processing by the tool implementation. Read shards should typically be much larger than the maximum assembly
 * region size to achieve good performance, and should have sufficient padding on either end to avoid boundary artifacts
 * for events near the shard boundaries.
 *
 * Since the assembly regions cover the padded regions around each read shard in addition to the shard's main span,
 * tool implementations may need to query {@link #getCurrentReadShardBounds} to determine whether an event is within
 * the main span or the padded regions, and to implement a scheme for avoiding reporting events more than once for events
 * that span shard boundaries.
 *
 * Read shards exist mainly as a proof-of-concept that we can shard the reads without introducing calling artifacts,
 * which will be important for the Spark equivalent of this traversal.
 */
public abstract class AssemblyRegionWalker extends GATKTool {

    @Argument(fullName="readShardSize", shortName="readShardSize", doc = "Maximum size of each read shard, in bases. For good performance, this should be much larger than the maximum assembly region size.", optional = true)
    protected int readShardSize = defaultReadShardSize();

    @Argument(fullName="readShardPadding", shortName="readShardPadding", doc = "Each read shard has this many bases of extra context on each side. Read shards must have as much or more padding than assembly regions.", optional = true)
    protected int readShardPadding = defaultReadShardPadding();

    @Argument(fullName = "minAssemblyRegionSize", shortName = "minAssemblyRegionSize", doc = "Minimum size of an assembly region", optional = true)
    protected int minAssemblyRegionSize = defaultMinAssemblyRegionSize();

    @Argument(fullName = "maxAssemblyRegionSize", shortName = "maxAssemblyRegionSize", doc = "Maximum size of an assembly region", optional = true)
    protected int maxAssemblyRegionSize = defaultMaxAssemblyRegionSize();

    @Argument(fullName = "assemblyRegionPadding", shortName = "assemblyRegionPadding", doc = "Amount of additional bases of context to include around each assembly region", optional = true)
    protected int assemblyRegionPadding = defaultAssemblyRegionPadding();

    @Advanced
    @Argument(fullName = "activeProbabilityThreshold", shortName = "activeProbabilityThreshold", doc="Minimum probability for a locus to be considered active.", optional = true)
    protected double activeProbThreshold = defaultActiveProbThreshold();

    @Advanced
    @Argument(fullName = "maxProbPropagationDistance", shortName = "maxProbPropagationDistance", doc="Upper limit on how many bases away probability mass can be moved around when calculating the boundaries between active and inactive assembly regions", optional = true)
    protected int maxProbPropagationDistance = defaultMaxProbPropagationDistance();

    @Argument(fullName = "disable_all_read_filters", shortName = "f", doc = "Disable all read filters", common = false, optional = true)
    public boolean disableAllReadFilters = false;

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
     * @return Default value for the {@link #activeProbThreshold} parameter, if none is provided on the command line
     */
    protected abstract double defaultActiveProbThreshold();

    /**
     * @return Default value for the {@link #maxProbPropagationDistance} parameter, if none is provided on the command line
     */
    protected abstract int defaultMaxProbPropagationDistance();

    @Override
    public final boolean requiresReads() { return true; }

    @Override
    public final boolean requiresReference() { return true; }

    private List<ReadShard> readShards;
    private ReadShard currentReadShard;

    /**
     * Initialize data sources for traversal.
     *
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalStart() instead.
     */
    @Override
    protected final void onStartup() {
        super.onStartup();

        if ( readShardSize <= 0 ) {
            throw new UserException.BadArgumentValue("read shard size must be > 0");
        }

        if ( readShardPadding < 0 ) {
            throw new UserException.BadArgumentValue("read shard padding must be >= 0");
        }

        if ( minAssemblyRegionSize > maxAssemblyRegionSize ) {
            throw new UserException.BadArgumentValue("minAssemblyRegionSize must be <= maxAssemblyRegionSize");
        }

        if ( maxAssemblyRegionSize > readShardSize ) {
            throw new UserException.BadArgumentValue("maxAssemblyRegionSize must be <= readShardSize");
        }

        if ( assemblyRegionPadding > readShardPadding ) {
            throw new UserException.BadArgumentValue("assemblyRegionPadding must be <= readShardPadding");
        }

        final List<SimpleInterval> intervals = hasIntervals() ? intervalsForTraversal : IntervalUtils.getAllIntervalsForReference(getHeaderForReads().getSequenceDictionary());
        readShards = makeReadShards(intervals);
    }

    /**
     * Shard our intervals for traversal into ReadShards using the {@link #readShardSize} and {@link #readShardPadding} arguments
     *
     * @param intervals unmodified intervals for traversal
     * @return List of {@link ReadShard} objects, sharded and padded as necessary
     */
    private List<ReadShard> makeReadShards( final List<SimpleInterval> intervals ) {
        final List<ReadShard> shards = new ArrayList<>();

        for ( final SimpleInterval interval : intervals ) {
            shards.addAll(ReadShard.divideIntervalIntoShards(interval, readShardSize, readShardPadding, reads, getHeaderForReads().getSequenceDictionary()));
        }

        return shards;
    }

    /**
     * Returns the read filter (simple or composite) that will be applied to the reads.
     *
     * The default implementation uses the {@link org.broadinstitute.hellbender.engine.filters.WellformedReadFilter} filter with all default options,
     * as well as the {@link org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary#MAPPED} filter.
     *
     * Default implementation of {@link #traverse()} calls this method once before iterating
     * over the reads and reuses the filter object to avoid object allocation. Nevertheless, keeping state in filter objects is strongly discouraged.
     *
     * Subclasses can override to provide their own filters
     * Multiple filters can be composed by using {@link org.broadinstitute.hellbender.engine.filters.ReadFilter} composition methods.
     */
    public CountingReadFilter makeReadFilter(){
        return new CountingReadFilter("Wellformed", new WellformedReadFilter(getHeaderForReads()))
                .and(new CountingReadFilter("Mapped", ReadFilterLibrary.MAPPED));
    }

    /**
     * @return The boundaries of the read shard we're currently operating within (ignoring any padding).
     */
    public SimpleInterval getCurrentReadShardBounds() {
        return currentReadShard.getInterval();
    }

    @Override
    public final void traverse() {
        CountingReadFilter countedFilter = disableAllReadFilters ?
                new CountingReadFilter("Allow all", ReadFilterLibrary.ALLOW_ALL_READS ) :
                makeReadFilter();

        // Since we're processing regions rather than individual reads, tell the progress
        // meter to check the time more frequently (every 10 regions instead of every 1000 regions).
        progressMeter.setRecordsBetweenTimeChecks(10L);

        for ( final ReadShard readShard : readShards ) {
            // Since reads in each shard are lazily fetched, we need to pass the filter to the window
            // instead of filtering the reads directly here
            readShard.setReadFilter(countedFilter);
            currentReadShard = readShard;

            processReadShard(readShard,
                    new ReferenceContext(reference, readShard.getPaddedInterval()), // use the fully-padded window to fetch overlapping data
                    new FeatureContext(features, readShard.getPaddedInterval()));
        }

        logger.info(countedFilter.getSummaryLine());
    }

    /**
     * Divide the given ReadShard up into active/inactive AssemblyRegions using the {@link #assemblyRegionEvaluator},
     * and send each region to the tool implementation for processing.
     *
     * @param shard ReadShard to process
     * @param referenceContext Reference bases spanning the fully-padded interval of the shard
     * @param featureContext Features spanning the fully-padded interval of the shard
     */
    private void processReadShard( ReadShard shard, ReferenceContext referenceContext, FeatureContext featureContext ) {
        // Divide each shard into one or more assembly regions using our AssemblyRegionEvaluator:
        final Iterable<AssemblyRegion> assemblyRegions = AssemblyRegion.createFromReadShard(shard,
                getHeaderForReads(), referenceContext, featureContext, assemblyRegionEvaluator(),
                minAssemblyRegionSize, maxAssemblyRegionSize, assemblyRegionPadding, activeProbThreshold,
                maxProbPropagationDistance);

        // Call into the tool implementation to process each assembly region from this shard.
        for ( final AssemblyRegion assemblyRegion : assemblyRegions ) {
            logger.debug("Processing assembly region at " + assemblyRegion.getSpan() + " isActive: " + assemblyRegion.isActive() + " numReads: " + assemblyRegion.getReads().size() + " in read shard " + shard.getInterval());

            apply(assemblyRegion,
                    new ReferenceContext(reference, assemblyRegion.getExtendedSpan()),
                    new FeatureContext(features, assemblyRegion.getExtendedSpan()));

            // For this traversal, the progress meter unit is the assembly region rather than the read shard
            progressMeter.update(assemblyRegion.getSpan());
        }
    }

    /**
     * Shutdown data sources.
     *
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalDone() instead.
     */
    @Override
    protected final void onShutdown() {
        // Overridden only to make final so that concrete tool implementations don't override
        super.onShutdown();
    }

    /**
     * @return The evaluator to be used to determine whether each locus is active or not. Must be implemented by tool authors.
     *         The results of this per-locus evaluator are used to determine the bounds of each active and inactive region.
     */
    public abstract AssemblyRegionEvaluator assemblyRegionEvaluator();

    /**
     * Process an individual AssemblyRegion. Must be implemented by tool authors.
     *
     * Each region will come pre-marked as either "active" or "inactive" using the results of the configured
     * {@link #assemblyRegionEvaluator}. This method will be called once for each active AND inactive region,
     * and it is up to the implementation how to handle/process active vs. inactive regions.
     *
     * @param region region to process (pre-marked as either active or inactive)
     * @param referenceContext reference data overlapping the full extended span of the assembly region
     * @param featureContext features overlapping the full extended span of the assembly region
     */
    public abstract void apply( final AssemblyRegion region, final ReferenceContext referenceContext, final FeatureContext featureContext );
}
