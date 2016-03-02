package org.broadinstitute.hellbender.engine;


import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
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
 * A ReadWindowWalker is a tool that processes entire regions/windows of reads at a time. These windows may or may not be overlapping.
 *
 * Reads in each window are lazily fetched rather than loaded all at once (though tools may choose to load all the reads in a window
 * at once).
 *
 * {@link #windowSize} controls the maximum size of each window. Large intervals are sharded into regions of up to this size,
 * while intervals smaller than this size are not subdivided.
 *
 * {@link #windowPadding} controls the number of bases of extra context to fetch on either side of the window. Overlapping reads
 * are fetched for the padded regions, but they are marked in the {@link ReadWindow} objects as not being within the window proper.
 *
 * If multiple sources of reads are specified, they are merged together into a single sorted stream of reads.
 *
 * ReadWindowWalker authors must implement the {@link #apply}, {@link #defaultWindowSize}, and {@link #defaultWindowPadding} methods,
 * and may optionally implement {@link #onTraversalStart} and/or {@link @onTraversalDone}.
 */
public abstract class ReadWindowWalker extends GATKTool {

    @Argument(fullName="windowSize", shortName="windowSize", doc = "Maximum size of each window, in bases", optional = true)
    protected int windowSize = defaultWindowSize();

    @Argument(fullName="windowPadding", shortName="windowPadding", doc = "Each window has this many bases of extra context on each side", optional = true)
    protected int windowPadding = defaultWindowPadding();

    @Argument(fullName = "disable_all_read_filters", shortName = "f", doc = "Disable all read filters", common = false, optional = true)
    public boolean disableAllReadFilters = false;

    private List<ReadWindow> windows;

    @Override
    public boolean requiresReads() {
        return true;
    }

    /**
     * @return Default value for the {@link #windowSize} parameter, if none is provided on the command line
     */
    protected abstract int defaultWindowSize();

    /**
     * @return Default value for the {@link #windowPadding} parameter, if none is provided on the command line
     */
    protected abstract int defaultWindowPadding();

    /**
     * Initialize data sources.
     *
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalDone() instead.
     */
    @Override
    protected final void onStartup() {
        super.onStartup();

        if ( windowSize <= 0 ) {
            throw new UserException.BadArgumentValue("--windowSize", Integer.toString(windowSize), "window size must be > 0");
        }
        if ( windowPadding < 0 ) {
            throw new UserException.BadArgumentValue("--windowPadding", Integer.toString(windowPadding), "window padding must be >= 0");
        }

        final SAMSequenceDictionary readsDictionary = getHeaderForReads().getSequenceDictionary();
        final List<SimpleInterval> intervals = hasIntervals() ? intervalsForTraversal : IntervalUtils.getAllIntervalsForReference(readsDictionary);
        windows = makeWindows(intervals, windowSize, windowPadding, reads, readsDictionary);
    }

    /**
     * Shard our intervals for traversal into windows using the {@link #windowSize} and {@link #windowPadding} arguments
     *
     * @param intervals unmodified intervals for traversal
     * @param windowSize desired window size; intervals larger than this will be sharded into windows up to this size
     * @param windowPadding desired window padding; each windowed interval will be padded on both sides by this number of bases
     * @param readsSource data source for reads
     * @param dictionary sequence dictionary for reads
     * @return List of {@link ReadWindow} objects, sharded and padded as necessary
     */
    @VisibleForTesting
    protected static final List<ReadWindow> makeWindows( final List<SimpleInterval> intervals, final int windowSize, final int windowPadding, final ReadsDataSource readsSource, final SAMSequenceDictionary dictionary ) {
        final List<ReadWindow> windows = new ArrayList<>();

        for ( final SimpleInterval interval : intervals ) {
            windows.addAll(shardInterval(interval, windowSize, windowPadding, readsSource, dictionary));
        }

        return windows;
    }

    /**
     * Shard a single interval into windows using the {@link #windowSize} and {@link #windowPadding} arguments
     *
     * @param interval interval to shard; must be on the contig according to the provided dictionary
     * @param windowSize desired window size; intervals larger than this will be sharded into windows up to this size
     * @param windowPadding desired window padding; each windowed interval will be padded on both sides by this number of bases
     * @param readsSource data source for reads
     * @param dictionary sequence dictionary for reads
     * @return List of {@link ReadWindow} objects, sharded and padded as necessary
     */
    @VisibleForTesting
    protected static List<ReadWindow> shardInterval( final SimpleInterval interval, final int windowSize, final int windowPadding, final ReadsDataSource readsSource, final SAMSequenceDictionary dictionary ) {
        if ( ! IntervalUtils.intervalIsOnDictionaryContig(interval, dictionary) ) {
            throw new IllegalArgumentException("Interval " + interval + " not within the bounds of a contig in the provided dictionary");
        }

        final List<ReadWindow> shards = new ArrayList<>();

        int start = interval.getStart();

        while ( start <= interval.getEnd() ) {
            int end = Math.min(start + windowSize - 1, interval.getEnd());

            final SimpleInterval nextWindow = new SimpleInterval(interval.getContig(), start, end);
            final SimpleInterval nextWindowPadded = nextWindow.expandWithinContig(windowPadding, dictionary);
            shards.add(new ReadWindow(nextWindow, nextWindowPadded, readsSource));

            start += windowSize;
        }

        return shards;
    }

    /**
     * Returns the read filter (simple or composite) that will be applied to the reads in each window.
     *
     * The default implementation uses the {@link org.broadinstitute.hellbender.engine.filters.WellformedReadFilter} filter with all default options,
     * as well as the {@link ReadFilterLibrary#MAPPED} filter.
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

    @Override
    public void traverse() {
        CountingReadFilter countedFilter = disableAllReadFilters ?
                new CountingReadFilter("Allow all", ReadFilterLibrary.ALLOW_ALL_READS ) :
                makeReadFilter();

        // Since we're processing regions rather than individual reads, tell the progress
        // meter to check the time more frequently (every 10 regions instead of every 1000 regions).
        progressMeter.setRecordsBetweenTimeChecks(10L);

        for ( final ReadWindow window : windows ) {
            // Since reads in each window are lazily fetched, we need to pass the filter to the window
            // instead of filtering the reads directly here
            window.setReadFilter(countedFilter);

            apply(window,
                  new ReferenceContext(reference, window.getPaddedInterval()), // use the fully-padded window to fetch overlapping data
                  new FeatureContext(features, window.getPaddedInterval()));

            progressMeter.update(window.getInterval());
        }

        logger.info(countedFilter.getSummaryLine());
    }

    /**
     * Process a window of reads, potentially with overlapping reference/Feature data
     *
     * @param window window of reads to process; reads in this window are lazily loaded rather than pre-loaded
     * @param referenceContext reference bases overlapping the window (including padded regions)
     * @param featureContext Features overlapping the window (including padded regions)
     */
    public abstract void apply( final ReadWindow window, final ReferenceContext referenceContext, final FeatureContext featureContext );

    /**
     * Close data sources.
     *
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalDone() instead.
     */
    @Override
    protected final void onShutdown() {
        // Overridden only to make final so that concrete tool implementations don't override
        super.onShutdown();
    }
}
