package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ShardingArgumentCollection;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.iterators.ShardingIterator;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.List;

/**
 * A {@link ReadSliderWalker} is a tool that processes a single {@link Shard<GATKRead>} window over
 * the genome/intervals at a time, with the ability to query overlapping optional sources (reference
 * data and/or variants/features).
 *
 * <p>Windows are constructed over each interval, using the window-size, window-step and window-padding
 * provided by the {@link #getShardingArgs()} arguments.
 *
 * <p>{@link ReadSliderWalker} authors must implement the {@link #apply(Shard, ReferenceContext, FeatureContext)}
 * and {@link #getShardingArgs()} methods, and optionally implement {@link #onTraversalStart()}
 * and/or {@link #onTraversalSuccess()}. For closing other sources, such as output writers,
 * {@link #closeTool()} might be implemented.
 *
 * <p>For an example walker, see {@link org.broadinstitute.hellbender.tools.examples.ExampleReadSliderWalker}.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public abstract class ReadSliderWalker extends GATKTool {

    /**
     * Argument collection for sliding window traversal.
     */
    @ArgumentCollection
    protected final ShardingArgumentCollection shardingArgs = getShardingArgs();

    /**
     * Returns the arguments for sliding over the genome/intervals.
     */
    public abstract ShardingArgumentCollection getShardingArgs();

    /**
     * Default label is "windows".
     */
    @Override
    public String getProgressMeterRecordLabel() {
        return "windows";
    }

    /** Always require reads. */
    @Override
    public final boolean requiresReads() {
        return true;
    }

    // list of shard boundaries to iterate over
    private List<ShardBoundary> windows;

    /**
     * Initializes windows and data sources for traversal.
     *
     * <p>Marked as final so that tool authors don't override it. Tool authors should override
     * {@link #onTraversalStart()} instead.
     */
    @Override
    protected final void onStartup() {
        super.onStartup();

        // validate first the sharding arguments
        shardingArgs.validate();

        // get the dictionary and check that it is provided
        final SAMSequenceDictionary dictionary = getBestAvailableSequenceDictionary();
        if (dictionary == null) {
            throw new UserException(String.format("Could not find sequence dictionary in data sources.A dictionary file is necessary for {}",
                    this.getToolName()));
        }

        // set traversal bounds if necessary
        if (hasIntervals() && hasReads()) {
            reads.setTraversalBounds(intervalArgumentCollection.getTraversalParameters(dictionary));
        }

        // generate the windows
        windows = makeWindows((hasIntervals()) ? intervalsForTraversal : IntervalUtils.getAllIntervalsForReference(dictionary),
                dictionary,
                shardingArgs.getWindowSize(), shardingArgs.getWindowStep(), shardingArgs.getWindowPad());
    }

    // generate the windows
    private static List<ShardBoundary> makeWindows(final List<SimpleInterval> intervals,
            final SAMSequenceDictionary dictionary,
            final int windowSize, final int windowStep, final int windowPad) {
        final List<ShardBoundary> windows = new ArrayList<>(intervals.size());
        for (final SimpleInterval i: intervals) {
            windows.addAll(Shard.divideIntervalIntoShards(i, windowSize, windowStep, windowPad, dictionary));
        }
        return windows;
    }

    /**
     * Implementation of window-based traversal.
     *
     * <p>Subclasses can override to provide their own behavior but default implementation should be suitable for most uses.
     *
     * <p>Default implementation iterates over the read source (with filters/transformers already applied)
     * by shards and pass all the available information to {@link #apply(Shard, ReferenceContext, FeatureContext)}.
     */
    @Override
    public void traverse() {
        final CountingReadFilter readFilter = makeReadFilter();
        final ShardingIterator<GATKRead> it = new ShardingIterator<>(getTransformedReadStream(readFilter).iterator(), windows, getBestAvailableSequenceDictionary());
        Utils.stream(it).forEach(shard -> {
            apply(shard,
                    new ReferenceContext(reference, shard.getPaddedInterval()),
                    new FeatureContext(features, shard.getPaddedInterval()));
            progressMeter.update(shard);
        });
        logger.info(readFilter.getSummaryLine());
    }

    /**
     * Process an individual window.
     *
     * <p>Must be implemented by tool authors. In general, tool authors should simply stream
     * their output from apply(), and maintain as little internal state as possible.
     *
     * <p>Warning: if overlapping windows are processed, modifying in-place reads on the shard
     * will affect the next iteration. In some cases this is desirable (e.g., cached per-read
     * statistic); otherwise, developers should perform a copy of the object before modification.
     *
     * @param reads            Non-null shard object with the window information (interval and
     *                         padded interval), containing reads overlapping it (might be empty if
     *                         no reads overlaps).
     * @param referenceContext Reference bases spanning the current window (including padded bases).
     *                         Will be an empty, but non-null, context object if there is no backing
     *                         source of reference data (in which case all queries on it will return
     *                         an empty array/iterator). Can request extra bases of context around
     *                         the current read's interval by invoking {@link ReferenceContext#setWindow}
     *                         on this object before calling {@link ReferenceContext#getBases}.
     * @param featureContext   Features spanning the current window (including padded bases). Will
     *                         be an empty, but non-null, context object if there is no backing
     *                         source of Feature data (in which case all queries on it will return
     *                         an empty List).
     */
    public abstract void apply(final Shard<GATKRead> reads, final ReferenceContext referenceContext, final FeatureContext featureContext);

    /**
     * Tool authors should override onTraversalDone() instead.
     *
     * <p>Marked final so that tool authors don't override it.
     */
    @Override
    protected final void onShutdown() {
        // Overridden only to make final so that concrete tool implementations don't override
        super.onShutdown();
    }
}
