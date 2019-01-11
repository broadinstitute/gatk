package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.downsampling.ReadsDownsampler;
import org.broadinstitute.hellbender.utils.downsampling.ReadsDownsamplingIterator;
import org.broadinstitute.hellbender.utils.iterators.ReadFilteringIterator;
import org.broadinstitute.hellbender.utils.iterators.ReadTransformingIterator;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;

/**
 * A class to represent shards of read data spanning multiple intervals.
 *
 * The reads are lazily loaded by default (when calling {@link #iterator}).
 * The reads returned will overlap the intervals, including any padding requested.
 *
 * Shards may be configured to apply custom filters, transformers and downsamplers to
 * the stream of reads returned by {@link #iterator}.
 *
 * IMPORTANT: For efficiency, all intervals within each shard are queried simultaneously.
 * This avoids the problem of decompressing the same file regions multiple times for
 * intervals that are close together, and is critical for performance!
 */
public final class MultiIntervalLocalReadShard implements MultiIntervalShard<GATKRead> {

    private final List<SimpleInterval> intervals;
    private final List<SimpleInterval> paddedIntervals;
    private final ReadsDataSource readsSource;

    private ReadTransformer preReadFilterTransformer;
    private ReadFilter readFilter;
    private ReadTransformer postReadFilterTransformer;
    private ReadsDownsampler downsampler;

    /**
     * Create a new MultiIntervalLocalReadShard spanning the given intervals, with each interval expanded
     * on both sides by the specified number of padding bases.
     *
     * @param intervals The intervals that this shard spans
     * @param intervalPadding Number of bases to pad each of the shard's intervals (on both sides)
     * @param readsSource Source of reads
     */
    public MultiIntervalLocalReadShard(final List<SimpleInterval> intervals, final int intervalPadding, final ReadsDataSource readsSource) {
        Utils.nonNull(intervals);
        Utils.nonNull(readsSource);
        Utils.validateArg(intervalPadding >= 0, "intervalPadding must be >= 0");

        // Feed intervals through IntervalUtils.getIntervalsWithFlanks() to ensure they get sorted using
        // the same comparator as the paddedIntervals below.
        this.intervals = Collections.unmodifiableList(IntervalUtils.getIntervalsWithFlanks(intervals, 0, readsSource.getSequenceDictionary()));

        // This will both pad each interval and merge any intervals that are overlapping or adjacent after padding,
        // in addition to sorting the intervals
        this.paddedIntervals = Collections.unmodifiableList(IntervalUtils.getIntervalsWithFlanks(intervals, intervalPadding, readsSource.getSequenceDictionary()));

        this.readsSource = readsSource;
    }

    /**
     * Create a new MultiIntervalLocalReadShard spanning the given intervals, with no padding around
     * intervals.
     *
     * @param intervals The intervals that this shard spans
     * @param readsSource Source of reads
     */
    public MultiIntervalLocalReadShard(final List<SimpleInterval> intervals, final ReadsDataSource readsSource) {
        this(intervals, 0, readsSource);
    }

    @Override
    public List<SimpleInterval> getIntervals() {
        return intervals;
    }

    @Override
    public List<SimpleInterval> getPaddedIntervals() {
        return paddedIntervals;
    }

    /**
     * Reads in this shard will be transformed before filtering
     *
     * @param transformer read transformer to apply before read filtering
     */
    public void setPreReadFilterTransformer(final ReadTransformer transformer) {
        preReadFilterTransformer = transformer;
    }

    /**
     * Reads in this shard will be filtered using this filter before being returned.
     * Read filtering will be performed before any requested downsampling.
     *
     * @param filter filter to use (may be null, which signifies that no filtering is to be performed)
     */
    public void setReadFilter(final ReadFilter filter) {
        this.readFilter = filter;
    }

    /**
     * Reads in this shard will be downsampled using this downsampler before being returned.
     * Downsampling will be performed after any requested read filtering.
     *
     * @param downsampler downsampler to use (may be null, which signifies that no downsampling is to be performed)
     */
    public void setDownsampler(final ReadsDownsampler downsampler) {
        this.downsampler = downsampler;
    }

    /**
     * Reads in this shard will be transformed after filtering and before downsampling
     *
     * @param transformer read transformer to apply after read filtering and before downsampling
     */
    public void setPostReadFilterTransformer(final ReadTransformer transformer) {
        postReadFilterTransformer = transformer;
    }

    /**
     * @return an iterator over reads in this shard, as filtered using the configured read filter
     *         and downsampled using the configured downsampler; reads are lazily loaded rather than pre-loaded
     *
     * Note that any read filtering is always performed before any downsampling.
     */
    @Override
    public Iterator<GATKRead> iterator() {
        // Query all intervals in this shard at once. This is critical for performance, to avoid
        // decompressing the same blocks multiple times for intervals that are close together!
        readsSource.setTraversalBounds(paddedIntervals);
        Iterator<GATKRead> readsIterator = readsSource.iterator();

        if (preReadFilterTransformer != null) {
            readsIterator = new ReadTransformingIterator(readsIterator, preReadFilterTransformer);
        }

        if ( readFilter != null ) {
            readsIterator = new ReadFilteringIterator(readsIterator, readFilter);
        }

        if (postReadFilterTransformer != null) {
            readsIterator = new ReadTransformingIterator(readsIterator, postReadFilterTransformer);
        }

        if ( downsampler != null ) {
            readsIterator = new ReadsDownsamplingIterator(readsIterator, downsampler);
        }

        return readsIterator;
    }
}

