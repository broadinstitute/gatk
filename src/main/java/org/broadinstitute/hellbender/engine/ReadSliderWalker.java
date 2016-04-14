package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKCommandLinePluginDescriptor;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptor;
import org.broadinstitute.hellbender.cmdline.argumentcollections.SliderArgumentCollection;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * A ReadSliderWalker is a tool that processes a single {@link LocalReadShard} window over the genome/intervals at a time,
 * with the ability to query overlapping optional sources reference data, and/or variants/features.
 *
 * The windows are constructed over each interval with a concrete window and step size, and padding overlapping
 * sources if requested.
 *
 * SlidingWindow authors must implement the {@link #apply}, {@link #defaultWindowSize}, {@link #defaultWindowStep} and {@link
 * #defaultWindowPadding} methods, and may optionally implement {@link #onTraversalStart} and/or {@link #onTraversalSuccess}.
 *
 * For control the window arguments (window-size, window-step and window-padding) exposed to the user, authors may implement
 * {@link #provideWindowSizeArgument()}, {@link #provideWindowStepArgument()} or {@link #provideWindowPaddingArgument()}.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public abstract class ReadSliderWalker extends GATKTool {

    @Argument(fullName = "disable_all_read_filters", shortName = "f", doc = "Disable all read filters", common = false, optional = true)
    public boolean disableAllReadFilters = false;

    @ArgumentCollection
    protected SliderArgumentCollection sliderArgumentCollection = new SliderArgumentCollection(
            defaultWindowSize(), provideWindowSizeArgument(), defaultWindowStep(), provideWindowStepArgument(), defaultWindowPadding(), provideWindowPaddingArgument()
    );

    @Override
    public final boolean requiresReads() {
        return true;
    }

    /**
     * Should the tool provide a window-size argument? Tools that do should override to return {@code true}.
     */
    protected boolean provideWindowSizeArgument() {
        return false;
    }

    /**
     * Should the tool provide a window-step argument? Tools that do should override to return {@code true}.
     */
    protected boolean provideWindowStepArgument() {
        return false;
    }

    /**
     * Should the tool provide a window-padding argument? Tools that do should override to return {@code true}.
     */
    protected boolean provideWindowPaddingArgument() {
        return false;
    }

    /**
     * Get the default window size for making the windows
     *
     * @return window-size (must be a positive integer)
     */
    protected abstract int defaultWindowSize();

    /**
     * Get the window step for making the windows
     *
     * @return window-step (must be a positive integer)
     */
    protected abstract int defaultWindowStep();

    /**
     * Get the length in each direction of the window to pad
     *
     * @return number of bases to pad (must be a positive integer or 0)
     */
    protected abstract int defaultWindowPadding();

    /**
     * Return the list of GATKCommandLinePluginDescriptors to be used for this tool.
     * Uses the read filter plugin.
     */
    protected List<? extends GATKCommandLinePluginDescriptor<?>> getPluginDescriptors() {
        return Collections.singletonList(new GATKReadFilterPluginDescriptor(getDefaultReadFilters()));
    }


    /**
     * Returns the read filter (simple or composite) that will be applied to the reads before calling {@link #apply}.
     * The default implementation combines the default read filters for this tool (returned by
     * {@link org.broadinstitute.hellbender.engine.ReadWalker#getDefaultReadFilters} with any read filter command
     * line arguments specified by the user; wraps each filter in the resulting list with a CountingReadFilter;
     * and returns a single composite filter resulting from the list by and'ing them together.
     * <p>
     * Default tool implementation of {@link #traverse()} calls this method once before iterating
     * over the reads and reuses the filter object to avoid object allocation. Nevertheless, keeping state in filter
     * objects is strongly discouraged.
     * <p>
     * Multiple filters can be composed by using {@link org.broadinstitute.hellbender.engine.filters.ReadFilter}
     * composition methods.
     */
    public CountingReadFilter makeReadFilter() {
        final GATKReadFilterPluginDescriptor readFilterPlugin =
                commandLineParser.getPluginDescriptor(GATKReadFilterPluginDescriptor.class);
        return readFilterPlugin.getMergedCountingReadFilter(getHeaderForReads());
    }

    /**
     * Returns the default list of CommandLineReadFilters that are used for this tool. The filters returned
     * by this method are subject to selective enabling/disabling by the user via the command line. The
     * default implementation uses the {@link WellformedReadFilter} filter with all default options. Subclasses
     * can override to provide alternative filters.
     * <p>
     * Note: this method is called before command line parsing begins, and thus before a SAMFileHeader is
     * available through {link #getHeaderForReads}.
     *
     * @return List of individual filters to be applied for this tool.
     */
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(new WellformedReadFilter());
    }

    private List<LocalReadShard> windows;

    /**
     * Initialize data sources for traversal.
     *
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalStart() instead.
     */
    @Override
    protected final void onStartup() {
        super.onStartup();
        SAMSequenceDictionary dictionary = getBestAvailableSequenceDictionary();
        // the tool needs a sequence dictionary to get the windows
        if (dictionary == null) {
            throw new UserException("Tool " + getClass().getSimpleName() + " requires some source for sequence dictionary, but none were provided");
        }
        final List<SimpleInterval> intervals = hasIntervals() ? intervalsForTraversal : IntervalUtils.getAllIntervalsForReference(dictionary);
        windows = makeReadShards(intervals,
                sliderArgumentCollection.getWindowSize(),
                sliderArgumentCollection.getWindowStep(),
                sliderArgumentCollection.getWindowPadding(),
                dictionary);
    }

    /**
     * Shard our intervals for traversal into LocalReadShard
     *
     * @param intervals     unmodified intervals for traversal
     * @param windowSize    the size for each window
     * @param windowStep    the step for each window
     * @param windowPadding the padding around the window
     * @return List of {@link LocalReadShard} objects, sharded and padded as necessary
     */
    private List<LocalReadShard> makeReadShards(final List<SimpleInterval> intervals, final int windowSize, final int windowStep, final int windowPadding, final SAMSequenceDictionary dictionary) {
        final List<LocalReadShard> shards = new ArrayList<>();
        for (final SimpleInterval interval : intervals) {
            shards.addAll(LocalReadShard.divideIntervalIntoShards(interval, windowSize, windowStep, windowPadding, reads, dictionary));
        }
        return shards;
    }

    /**
     * Implementation of window-based traversal. Subclasses can override to provide their own behavior but default
     * implementation should be suitable for most uses.
     *
     * The default implementation iterates over every LocalReadShard stored in {@link #windows} and pass all the available
     * information to {@link #apply}
     */
    @Override
    public void traverse() {
        final CountingReadFilter countedFilter = makeReadFilter();
        // Since we're processing regions rather than individual reads, tell the progress meter to check the time more frequently (every 10 regions instead of every 1000 regions).
        progressMeter.setRecordsBetweenTimeChecks(10L);
        // iterate over the windows
        for (final LocalReadShard window : windows) {
            // Since reads in each window are lazily fetched, we need to pass the filter to the window instead of filtering the reads directly here
            window.setReadFilter(countedFilter);
            apply(window,
                    new ReferenceContext(reference, window.getPaddedInterval()), // use the fully-padded window to fetch overlapping data
                    new FeatureContext(features, window.getPaddedInterval()));
            progressMeter.update(window.getInterval());
        }
        logger.info(countedFilter.getSummaryLine());
    }

    /**
     * Process an individual window. Must be implemented by tool authors. In general, tool authors should simply stream
     * their output from apply(), and maintain as little internal state as possible.
     *
     * @param readShard        Reads overlapping the current window (including padded regions). Will be an empty, but
     *                         non-null, context object if there is no backing source of reads data (in which case all
     *                         queries on it will return an empty array/iterator).
     * @param referenceContext Reference bases spanning the current window (including padded regions). Will be an empty,
     *                         but non-null, context object if there is no backing source of reference data (in which
     *                         case all queries on it will return an empty array/iterator). Can request extra bases of
     *                         context around the current read's interval by invoking {@link ReferenceContext#setWindow}
     *                         on this object before calling {@link ReferenceContext#getBases}
     * @param featureContext   Features spanning the current window (including padded regions). Will be an empty, but
     *                         non-null, context object if there is no backing source of Feature data (in which case all
     *                         queries on it will return an empty List).
     */
    public abstract void apply(Shard<GATKRead> readShard, ReferenceContext referenceContext, FeatureContext featureContext);

    /**
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalDone() instead.
     */
    @Override
    protected final void onShutdown() {
        // Overridden only to make final so that concrete tool implementations don't override
        super.onShutdown();
    }
}
