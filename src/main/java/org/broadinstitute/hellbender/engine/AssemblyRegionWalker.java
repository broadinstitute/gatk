package org.broadinstitute.hellbender.engine;

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
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.downsampling.PositionalDownsampler;
import org.broadinstitute.hellbender.utils.downsampling.ReadsDownsampler;

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
 * for processing by the tool implementation. One read shard is created per contig.
 */
public abstract class AssemblyRegionWalker extends GATKTool {

    //NOTE: these argument names are referenced by HaplotypeCallerSpark
    public static final String MIN_ASSEMBLY_LONG_NAME = "min-assembly-region-size";
    public static final String MAX_ASSEMBLY_LONG_NAME = "max-assembly-region-size";
    public static final String ASSEMBLY_PADDING_LONG_NAME = "assembly-region-padding";
    public static final String MAX_STARTS_LONG_NAME = "max-reads-per-alignment-start";
    public static final String THRESHOLD_LONG_NAME = "active-probability-threshold";
    public static final String PROPAGATION_LONG_NAME = "max-prob-propagation-distance";
    public static final String PROFILE_OUT_LONG_NAME = "activity-profile-out";
    public static final String ASSEMBLY_REGION_OUT_LONG_NAME = "assembly-region-out";

    @Advanced
    @Argument(fullName = MIN_ASSEMBLY_LONG_NAME, doc = "Minimum size of an assembly region", optional = true)
    protected int minAssemblyRegionSize = defaultMinAssemblyRegionSize();

    @Advanced
    @Argument(fullName = MAX_ASSEMBLY_LONG_NAME, doc = "Maximum size of an assembly region", optional = true)
    protected int maxAssemblyRegionSize = defaultMaxAssemblyRegionSize();

    @Advanced
    @Argument(fullName = ASSEMBLY_PADDING_LONG_NAME, doc = "Number of additional bases of context to include around each assembly region", optional = true)
    protected int assemblyRegionPadding = defaultAssemblyRegionPadding();

    @Argument(fullName = MAX_STARTS_LONG_NAME, doc = "Maximum number of reads to retain per alignment start position. Reads above this threshold will be downsampled. Set to 0 to disable.", optional = true)
    protected int maxReadsPerAlignmentStart = defaultMaxReadsPerAlignmentStart();

    @Advanced
    @Argument(fullName = THRESHOLD_LONG_NAME, doc="Minimum probability for a locus to be considered active.", optional = true)
    protected double activeProbThreshold = defaultActiveProbThreshold();

    @Advanced
    @Argument(fullName = PROPAGATION_LONG_NAME, doc="Upper limit on how many bases away probability mass can be moved around when calculating the boundaries between active and inactive assembly regions", optional = true)
    protected int maxProbPropagationDistance = defaultMaxProbPropagationDistance();

    /**
     * If provided, this walker will write out its activity profile (per bp probabilities of being active)
     * to this file in the IGV formatted TAB deliminated output:
     *
     * http://www.broadinstitute.org/software/igv/IGV
     *
     * Intended to make debugging the activity profile calculations easier
     */
    @Argument(fullName = PROFILE_OUT_LONG_NAME, doc="Output the raw activity profile results in IGV format", optional = true)
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
    @Argument(fullName = ASSEMBLY_REGION_OUT_LONG_NAME, doc="Output the assembly region to this IGV formatted file", optional = true)
    protected String assemblyRegionOut = null;

    private PrintStream assemblyRegionOutStream;

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

    private List<MultiIntervalLocalReadShard> readShards;

    /**
     * Initialize data sources for traversal.
     *
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalStart() instead.
     */
    @Override
    protected final void onStartup() {
        super.onStartup();

        if ( minAssemblyRegionSize <= 0 || maxAssemblyRegionSize <= 0 ) {
            throw new CommandLineException.BadArgumentValue("min/max assembly region size must be > 0");
        }

        if ( minAssemblyRegionSize > maxAssemblyRegionSize ) {
            throw new CommandLineException.BadArgumentValue("minAssemblyRegionSize must be <= maxAssemblyRegionSize");
        }

        if ( assemblyRegionPadding < 0 ) {
            throw new CommandLineException.BadArgumentValue("assemblyRegionPadding must be >= 0");
        }

        if ( maxReadsPerAlignmentStart < 0 ) {
            throw new CommandLineException.BadArgumentValue("maxReadsPerAlignmentStart must be >= 0");
        }

        final List<SimpleInterval> intervals = hasUserSuppliedIntervals() ? userIntervals : IntervalUtils.getAllIntervalsForReference(getHeaderForReads().getSequenceDictionary());
        readShards = makeReadShards(intervals);

        initializeAssemblyRegionOutputStreams();
    }

    /**
     * Shard our intervals for traversal into ReadShards, each shard containing all of the
     * intervals for one contig.
     *
     * We pad the intervals within each shard by the same amount as the assembly region padding
     * to avoid boundary artifacts.
     *
     * @param intervals unmodified intervals for traversal
     * @return List of {@link MultiIntervalLocalReadShard} objects, sharded and padded as necessary
     */
    private List<MultiIntervalLocalReadShard> makeReadShards(final List<SimpleInterval> intervals ) {
        final List<MultiIntervalLocalReadShard> shards = new ArrayList<>();
        final List<List<SimpleInterval>> intervalsGroupedByContig = IntervalUtils.groupIntervalsByContig(intervals);

        for ( final List<SimpleInterval> allIntervalsOnContig : intervalsGroupedByContig ) {
            shards.add(new MultiIntervalLocalReadShard(allIntervalsOnContig, assemblyRegionPadding, reads));
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

        for ( final MultiIntervalLocalReadShard readShard : readShards ) {
            // Since reads in each shard are lazily fetched, we need to pass the filter and transformers to the window
            // instead of filtering the reads directly here
            readShard.setPreReadFilterTransformer(makePreReadFilterTransformer());
            readShard.setReadFilter(countedFilter);
            readShard.setDownsampler(createDownsampler());
            readShard.setPostReadFilterTransformer(makePostReadFilterTransformer());

            processReadShard(readShard, reference, features);
        }

        logger.info(countedFilter.getSummaryLine());
    }

    /**
     * Divide the given Shard up into active/inactive AssemblyRegions using the {@link #assemblyRegionEvaluator},
     * and send each region to the tool implementation for processing.
     *
     * @param shard MultiIntervalLocalReadShard to process
     * @param reference Reference data source
     * @param features FeatureManager
     */
    private void processReadShard(MultiIntervalLocalReadShard shard, ReferenceDataSource reference, FeatureManager features ) {
        final Iterator<AssemblyRegion> assemblyRegionIter = new AssemblyRegionIterator(shard, getHeaderForReads(), reference, features, assemblyRegionEvaluator(), minAssemblyRegionSize, maxAssemblyRegionSize, assemblyRegionPadding, activeProbThreshold, maxProbPropagationDistance, includeReadsWithDeletionsInIsActivePileups());

        // Call into the tool implementation to process each assembly region from this shard.
        while ( assemblyRegionIter.hasNext() ) {
            final AssemblyRegion assemblyRegion = assemblyRegionIter.next();
            
            logger.debug("Processing assembly region at " + assemblyRegion.getSpan() + " isActive: " + assemblyRegion.isActive() + " numReads: " + assemblyRegion.getReads().size());
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
