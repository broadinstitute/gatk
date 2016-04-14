package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.downsampling.ReadsDownsampler;
import org.broadinstitute.hellbender.utils.downsampling.ReadsDownsamplingIterator;
import org.broadinstitute.hellbender.utils.iterators.ReadFilteringIterator;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A class to represent a shard of reads data, optionally expanded by a configurable amount of padded data.
 *
 * The reads are lazily loaded by default (when accessing the reads via {@link #iterator}. Loading all the
 * reads in the shard at once is possible via {@link #loadAllReads}.
 *
 * The reads returned will overlap the expanded padded interval. It's possible to query whether they are within
 * the main part of the shard via {@link #contains} and {@link #containsStartPosition}.
 *
 * The reads in the shard can be filtered via {@link #setReadFilter} (no filtering is performed by default).
 */
public final class LocalReadShard implements Shard<GATKRead> {

    private final SimpleInterval interval;
    private final SimpleInterval paddedInterval;
    private final ReadsDataSource readsSource;
    private ReadFilter readFilter;
    private ReadsDownsampler downsampler;

    /**
     * Create a new Shard spanning the specified interval, with the specified amount of padding.
     *
     * @param interval the genomic span covered by this shard
     * @param paddedInterval the span covered by this shard, plus any additional padding on each side (must contain the un-padded interval)
     * @param readsSource source of reads from which to populate this shard
     */
    public LocalReadShard(final SimpleInterval interval, final SimpleInterval paddedInterval, final ReadsDataSource readsSource) {
        Utils.nonNull(interval);
        Utils.nonNull(paddedInterval);
        Utils.nonNull(readsSource);
        Utils.validateArg(paddedInterval.contains(interval), "The padded interval must contain the un-padded interval");

        this.interval = interval;
        this.paddedInterval = paddedInterval;
        this.readsSource = readsSource;
    }

    /**
     * Create a new Shard spanning the specified interval, with no additional padding
     *
     * @param interval the genomic span covered by this shard
     * @param readsSource source of reads from which to populate this shard
     */
    public LocalReadShard(final SimpleInterval interval, final ReadsDataSource readsSource) {
        this(interval, interval, readsSource);
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
     * @return the interval this shard spans
     */
    @Override
    public SimpleInterval getInterval() {
        return interval;
    }

    /**
     * @return the interval this shard spans, potentially with additional padding on each side
     */
    @Override
    public SimpleInterval getPaddedInterval() {
        return paddedInterval;
    }

    /**
     * @return an iterator over reads in this shard, as filtered using the configured read filter
     *         and downsampled using the configured downsampler; reads are lazily loaded rather than pre-loaded
     *
     * Note that any read filtering is always performed before any downsampling.
     */
    @Override
    public Iterator<GATKRead> iterator() {
        Iterator<GATKRead> readsIterator = readsSource.query(paddedInterval);

        if ( readFilter != null ) {
            readsIterator = new ReadFilteringIterator(readsIterator, readFilter);
        }

        if ( downsampler != null ) {
            readsIterator = new ReadsDownsamplingIterator(readsIterator, downsampler);
        }

        return readsIterator;
    }

    /**
     * @return a List containing all reads in this shard, pre-loaded, filtered using the configured read filter,
     *         and downsampled using the configured downsampler
     *
     * Call {@link #iterator} instead to avoid pre-loading all reads at once.
     *
     * Note that any read filtering is always performed before any downsampling.
     */
    public List<GATKRead> loadAllReads() {
        return loadAllRecords();
    }

    /**
     * Divide an interval into LocalReadShards. Each shard will cover up to shardSize bases, include shardPadding
     * bases of extra padding on either side, and begin shardSize bases after the previous shard (ie., shards will
     * not overlap except potentially in the padded regions).
     *
     * @param interval interval to shard; must be on the contig according to the provided dictionary
     * @param shardSize desired shard size; intervals larger than this will be divided into shards of up to this size
     * @param shardPadding desired shard padding; each shard's interval will be padded on both sides by this number of bases (may be 0)
     * @param readsSource data source for reads
     * @param dictionary sequence dictionary for reads
     * @return List of {@link LocalReadShard} objects spanning the interval
     */
    public static List<LocalReadShard> divideIntervalIntoShards(final SimpleInterval interval, final int shardSize, final int shardPadding, final ReadsDataSource readsSource, final SAMSequenceDictionary dictionary) {
        return divideIntervalIntoShards(interval, shardSize, shardSize, shardPadding, readsSource, dictionary);
    }

    /**
     * Divide an interval into LocalReadShards. Each shard will cover up to shardSize bases, include shardPadding
     * bases of extra padding on either side, and begin shardStep bases after the previous shard.
     *
     * @param interval interval to shard; must be on the contig according to the provided dictionary
     * @param shardSize desired shard size; intervals larger than this will be divided into shards of up to this size
     * @param shardStep each shard will begin this many bases away from the previous shard
     * @param shardPadding desired shard padding; each shard's interval will be padded on both sides by this number of bases (may be 0)
     * @param readsSource data source for reads
     * @param dictionary sequence dictionary for reads
     * @return List of {@link LocalReadShard} objects spanning the interval
     */
    public static List<LocalReadShard> divideIntervalIntoShards(final SimpleInterval interval, final int shardSize, final int shardStep, final int shardPadding, final ReadsDataSource readsSource, final SAMSequenceDictionary dictionary) {
        Utils.nonNull(readsSource);
        return Shard.divideIntervalIntoShards(interval, shardSize, shardStep, shardPadding, dictionary)
                    .stream().map(shardBoundary -> new LocalReadShard(shardBoundary.getInterval(), shardBoundary.getPaddedInterval(), readsSource))
                    .collect(Collectors.toList());
    }
}
