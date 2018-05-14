package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ShardingArgumentCollection;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.examples.ExampleSlidingWindowReadWalker;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.iterators.ShardingIterator;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * A {@link SlidingWindowReadWalker} is a tool that processes a single {@link Shard<GATKRead>} window over
 * the genome/intervals at a time, with the ability to query overlapping optional sources (reference
 * data and/or variants/features).
 *
 * <p>Windows are constructed using the window size, step and padding provided by the
 * {@link #getShardingArgs()} arguments. Each interval provided by the user with the <code>-L</code>
 * option (or per-chromosome from the sequence dictionary if not provided) is divided into windows
 * up to <b>window size</b> base pairs if larger than that size. Each window will overlap
 * <b>window step</b> base pairs with the contiguous one (without window padding), and is extended
 * in both sides by <b>window padding</b> base pairs.
 *
 * <p>{@link SlidingWindowReadWalker} authors must implement the {@link #apply(Shard, ReferenceContext, FeatureContext)}
 * and {@link #getShardingArgs()} methods, and optionally implement {@link #onTraversalStart()}
 * and/or {@link #onTraversalSuccess()}. For closing other sources, such as output writers,
 * {@link #closeTool()} might be implemented.
 *
 * <p>For an example walker, see {@link ExampleSlidingWindowReadWalker}.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public abstract class SlidingWindowReadWalker extends GATKTool {

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

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        // include the MappedReadFilter
        return Arrays.asList(new WellformedReadFilter(), new ReadFilterLibrary.MappedReadFilter());
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
            throw new UserException(String.format("Could not find sequence dictionary in data sources. A dictionary file is necessary for %s",
                    this.getToolName()));
        }

        // set traversal bounds if necessary
        if (hasIntervals()) {
            reads.setTraversalBounds(intervalArgumentCollection.getTraversalParameters(dictionary));
        }

        // generate the windows
        windows = makeWindows((hasIntervals()) ? intervalsForTraversal : IntervalUtils.getAllIntervalsForReference(dictionary),
                dictionary,
                shardingArgs.getWindowSize(), shardingArgs.getWindowStep(), shardingArgs.getWindowPadding());
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
     * <p>Warning: if overlapping windows are processed, modifying reads in-place on the shard
     * will affect the next iteration. In some cases this is desirable (e.g., cached per-read
     * statistic); otherwise, developers should perform a copy of the object before modification.
     *
     * @param reads            Non-null shard object with the window information (interval and
     *                         padded interval), containing reads overlapping the padded interal
     *                         (might be empty if no reads overlaps).
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
