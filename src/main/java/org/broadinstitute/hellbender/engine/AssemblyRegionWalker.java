package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IGVUtils;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.downsampling.PositionalDownsampler;
import org.broadinstitute.hellbender.utils.downsampling.ReadsDownsampler;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
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

    /**
     * If readShardSize is set to this value, we will not shard the user's intervals. Instead,
     * we'll create one shard per interval (or one shard per contig, if intervals are not explicitly specified)
     */
    public static final int NO_INTERVAL_SHARDING = -1;

    @Advanced
    @Argument(fullName="readShardSize", shortName="readShardSize", doc = "Maximum size of each read shard, in bases. Set to " + NO_INTERVAL_SHARDING + " for one shard per interval (or one shard per contig, if intervals are not explicitly specified). For good performance, this should typically be much larger than the maximum assembly region size.", optional = true)
    protected int readShardSize = defaultReadShardSize();

    @Advanced
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
     * If provided, this walker will write out its activity profile (per bp probabilities of being active)
     * to this file in the IGV formatted TAB deliminated output:
     *
     * http://www.broadinstitute.org/software/igv/IGV
     *
     * Intended to make debugging the activity profile calculations easier
     */
    @Argument(fullName="activityProfileOut", shortName="APO", doc="Output the raw activity profile results in IGV format", optional = true)
    protected String activityProfileOut = null;

    private PrintStream activityProfileOutStream;

    /**
     * If provided, this walker will write out its assembly regions
     * to this file in the IGV formatted TAB-delimited output:
     *
     * http://www.broadinstitute.org/software/igv/IGV
     *
     * Intended to make debugging the active region calculations easier
     */
    @Argument(fullName="assemblyRegionOut", shortName="ARO", doc="Output the assembly region to this IGV formatted file", optional = true)
    protected String assemblyRegionOut = null;

    private PrintStream assemblyRegionOutStream;

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
     * @return If true, include reads with deletions at the current locus in the pileups passed to the AssemblyRegionEvaluator.
     */
    protected abstract boolean includeReadsWithDeletionsInIsActivePileups();

    @Override
    public final boolean requiresReads() { return true; }

    @Override
    public final boolean requiresReference() { return true; }

    @Override
    public String getProgressMeterRecordLabel() { return "regions"; }
    
    private List<LocalReadShard> readShards;
    private Shard<GATKRead> currentReadShard;

    /**
     * Initialize data sources for traversal.
     *
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalStart() instead.
     */
    @Override
    protected final void onStartup() {
        super.onStartup();

        if ( readShardSize <= 0 && readShardSize != NO_INTERVAL_SHARDING ) {
            throw new CommandLineException.BadArgumentValue("read shard size must be > 0 or " + NO_INTERVAL_SHARDING + " for no sharding");
        }

        if ( readShardPadding < 0 ) {
            throw new CommandLineException.BadArgumentValue("read shard padding must be >= 0");
        }

        if ( minAssemblyRegionSize > maxAssemblyRegionSize ) {
            throw new CommandLineException.BadArgumentValue("minAssemblyRegionSize must be <= maxAssemblyRegionSize");
        }

        if ( readShardSize != NO_INTERVAL_SHARDING && maxAssemblyRegionSize > readShardSize ) {
            throw new CommandLineException.BadArgumentValue("maxAssemblyRegionSize must be <= readShardSize");
        }

        if ( assemblyRegionPadding > readShardPadding ) {
            throw new CommandLineException.BadArgumentValue("assemblyRegionPadding must be <= readShardPadding");
        }

        final List<SimpleInterval> intervals = hasIntervals() ? intervalsForTraversal : IntervalUtils.getAllIntervalsForReference(getHeaderForReads().getSequenceDictionary());
        readShards = makeReadShards(intervals);

        initializeAssemblyRegionOutputStreams();
    }

    /**
     * Shard our intervals for traversal into ReadShards using the {@link #readShardSize} and {@link #readShardPadding} arguments
     *
     * @param intervals unmodified intervals for traversal
     * @return List of {@link LocalReadShard} objects, sharded and padded as necessary
     */
    private List<LocalReadShard> makeReadShards(final List<SimpleInterval> intervals ) {
        final List<LocalReadShard> shards = new ArrayList<>();

        for ( final SimpleInterval interval : intervals ) {
            if ( readShardSize != NO_INTERVAL_SHARDING ) {
                shards.addAll(LocalReadShard.divideIntervalIntoShards(interval, readShardSize, readShardPadding, reads, getHeaderForReads().getSequenceDictionary()));
            }
            else {
                shards.add(new LocalReadShard(interval, interval.expandWithinContig(readShardPadding, getHeaderForReads().getSequenceDictionary()), reads));
            }
        }

        return shards;
    }

    private void initializeAssemblyRegionOutputStreams() {
        if ( activityProfileOut != null ) {
            try {
                activityProfileOutStream = new PrintStream(activityProfileOut);
            }
            catch ( IOException e ) {
                throw new UserException.CouldNotCreateOutputFile(activityProfileOut, "Error writing activity profile to output file", e);
            }

            logger.info("Writing activity profile to " + activityProfileOut);
            IGVUtils.printIGVFormatHeader(activityProfileOutStream, "line", "ActivityProfile");
        }

        if ( assemblyRegionOut != null ) {
            try {
                assemblyRegionOutStream = new PrintStream(assemblyRegionOut);
            }
            catch ( IOException e ) {
                throw new UserException.CouldNotCreateOutputFile(assemblyRegionOut, "Error writing assembly regions to output file", e);
            }

            logger.info("Writing assembly regions to " + assemblyRegionOut);
            IGVUtils.printIGVFormatHeader(assemblyRegionOutStream, "line", "AssemblyRegions");
        }
    }

    /**
     * @return The boundaries of the read shard we're currently operating within (ignoring any padding).
     */
    public SimpleInterval getCurrentReadShardBounds() {
        return currentReadShard.getInterval();
    }

    /**
     * Returns the default list of CommandLineReadFilters that are used for this tool. The filters
     * returned by this method are subject to selective enabling/disabling and customization by the
     * user via the command line. The default implementation uses the {@link WellformedReadFilter}
     * filter with all default options, as well as the {@link ReadFilterLibrary.MappedReadFilter}.
     * Subclasses can override to provide alternative filters.
     *
     * Note: this method is called before command line parsing begins, and thus before a SAMFileHeader is
     * available through {link #getHeaderForReads}.
     *
     * @return List of default filter instances to be applied for this tool.
     */
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> defaultFilters = new ArrayList<>(2);
        defaultFilters.add(new WellformedReadFilter());
        defaultFilters.add(new ReadFilterLibrary.MappedReadFilter());
        return defaultFilters;
    }

    protected ReadsDownsampler createDownsampler() {
        return maxReadsPerAlignmentStart > 0 ? new PositionalDownsampler(maxReadsPerAlignmentStart, getHeaderForReads()) : null;
    }

    @Override
    public final void traverse() {

        CountingReadFilter countedFilter = makeReadFilter();

        // Since we're processing regions rather than individual reads, tell the progress
        // meter to check the time more frequently (every 10 regions instead of every 1000 regions).
        progressMeter.setRecordsBetweenTimeChecks(10L);

        for ( final LocalReadShard readShard : readShards ) {
            // Since reads in each shard are lazily fetched, we need to pass the filter and transformers to the window
            // instead of filtering the reads directly here
            readShard.setPreReadFilterTransformer(makePreReadFilterTransformer());
            readShard.setReadFilter(countedFilter);
            readShard.setDownsampler(createDownsampler());
            readShard.setPostReadFilterTransformer(makePostReadFilterTransformer());
            currentReadShard = readShard;

            processReadShard(readShard, reference, features);
        }

        logger.info(countedFilter.getSummaryLine());
    }

    /**
     * Divide the given Shard up into active/inactive AssemblyRegions using the {@link #assemblyRegionEvaluator},
     * and send each region to the tool implementation for processing.
     *
     * @param shard Shard to process
     * @param reference Reference data source
     * @param features FeatureManager
     */
    private void processReadShard(Shard<GATKRead> shard, ReferenceDataSource reference, FeatureManager features ) {
        final Iterator<AssemblyRegion> assemblyRegionIter = new AssemblyRegionIterator(shard, getHeaderForReads(), reference, features, assemblyRegionEvaluator(), minAssemblyRegionSize, maxAssemblyRegionSize, assemblyRegionPadding, activeProbThreshold, maxProbPropagationDistance, includeReadsWithDeletionsInIsActivePileups());

        // Call into the tool implementation to process each assembly region from this shard.
        while ( assemblyRegionIter.hasNext() ) {
            final AssemblyRegion assemblyRegion = assemblyRegionIter.next();
            
            logger.debug("Processing assembly region at " + assemblyRegion.getSpan() + " isActive: " + assemblyRegion.isActive() + " numReads: " + assemblyRegion.getReads().size() + " in read shard " + shard.getInterval());
            writeAssemblyRegion(assemblyRegion);

            apply(assemblyRegion,
                    new ReferenceContext(reference, assemblyRegion.getExtendedSpan()),
                    new FeatureContext(features, assemblyRegion.getExtendedSpan()));

            // For this traversal, the progress meter unit is the assembly region rather than the read shard
            progressMeter.update(assemblyRegion.getSpan());
        }
    }

    private void writeAssemblyRegion(final AssemblyRegion region) {
        writeActivityProfile(region.getSupportingStates());

        if ( assemblyRegionOutStream != null ) {
            IGVUtils.printIGVFormatRow(assemblyRegionOutStream, new SimpleInterval(region.getContig(), region.getStart(), region.getStart()),
                    "end-marker", 0.0);
            IGVUtils.printIGVFormatRow(assemblyRegionOutStream, region,
                    "size=" + new SimpleInterval(region).size(), region.isActive() ? 1.0 : -1.0);
        }
    }

    private void writeActivityProfile(final List<ActivityProfileState> states) {
        if ( activityProfileOutStream != null ) {
            for ( final ActivityProfileState state : states ) {
                IGVUtils.printIGVFormatRow(activityProfileOutStream, state.getLoc(), "state", Math.min(state.isActiveProb(), 1.0));
            }
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

        if ( assemblyRegionOutStream != null ) {
            assemblyRegionOutStream.close();
        }

        if ( activityProfileOutStream != null ) {
            activityProfileOutStream.close();
        }
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
