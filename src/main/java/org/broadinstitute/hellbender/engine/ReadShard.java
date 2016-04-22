package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.iterators.ReadFilteringIterator;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

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
public final class ReadShard implements Iterable<GATKRead>, Locatable {

    private final SimpleInterval interval;
    private final SimpleInterval paddedInterval;
    private final ReadsDataSource readsSource;
    private ReadFilter readFilter;

    /**
     * Create a new ReadShard spanning the specified interval, with the specified amount of padding.
     *
     * @param interval the genomic span covered by this shard
     * @param paddedInterval the span covered by this shard, plus any additional padding on each side (must contain the un-padded interval)
     * @param readsSource source of reads from which to populate this shard
     */
    public ReadShard( final SimpleInterval interval, final SimpleInterval paddedInterval, final ReadsDataSource readsSource ) {
        Utils.nonNull(interval);
        Utils.nonNull(paddedInterval);
        Utils.nonNull(readsSource);

        if ( ! paddedInterval.contains(interval) ) {
            throw new IllegalArgumentException("The padded interval must contain the un-padded interval");
        }

        this.interval = interval;
        this.paddedInterval = paddedInterval;
        this.readsSource = readsSource;
    }

    /**
     * Create a new ReadShard spanning the specified interval, with no additional padding
     *
     * @param interval the genomic span covered by this shard
     * @param readsSource source of reads from which to populate this shard
     */
    public ReadShard( final SimpleInterval interval, final ReadsDataSource readsSource ) {
        this(interval, interval, readsSource);
    }

    /**
     * Reads in this shard will be filtered using this filter before being returned
     *
     * @param filter filter to use (may be null)
     */
    public void setReadFilter( final ReadFilter filter ) {
        this.readFilter = filter;
    }

    /**
     * @return Contig this shard belongs to
     */
    @Override
    public String getContig() {
        return interval.getContig();
    }

    /**
     * @return Start position of this shard
     */
    @Override
    public int getStart() {
        return interval.getStart();
    }

    /**
     * @return End position of this shard
     */
    @Override
    public int getEnd() {
        return interval.getEnd();
    }

    /**
     * @return the interval this shard spans
     */
    public SimpleInterval getInterval() {
        return interval;
    }

    /**
     * @return the interval this shard spans, potentially with additional padding on each side
     */
    public SimpleInterval getPaddedInterval() {
        return paddedInterval;
    }

    /**
     * @return number of bases of padding to the left of our interval
     */
    public int numLeftPaddingBases() {
        return interval.getStart() - paddedInterval.getStart();
    }

    /**
     * @return number of bases of padding to the right of our interval
     */
    public int numRightPaddingBases() {
        return paddedInterval.getEnd() - interval.getEnd();
    }

    /**
     * @param loc Locatable to test
     * @return true if loc is completely contained within this shard's interval, otherwise false
     */
    public boolean contains( final Locatable loc ) {
        Utils.nonNull(loc);
        return interval.contains(loc);
    }

    /**
     * @param loc Locatable to test
     * @return true if loc starts within this shard's interval, otherwise false
     */
    public boolean containsStartPosition( final Locatable loc ) {
        Utils.nonNull(loc);
        return interval.contains(new SimpleInterval(loc.getContig(), loc.getStart(), loc.getStart()));
    }

    /**
     * @return an iterator over reads in this shard, as filtered using the configured read filter;
     *         reads are lazily loaded rather than pre-loaded
     */
    @Override
    public Iterator<GATKRead> iterator() {
        final Iterator<GATKRead> readsIterator = readsSource.query(paddedInterval);

        return readFilter != null ? new ReadFilteringIterator(readsIterator, readFilter) : readsIterator;
    }

    /**
     * @return a List containing all reads in this shard, pre-loaded, and filtered using the configured read filter
     *
     * note: call {@link #iterator} instead to avoid pre-loading all reads at once
     */
    public List<GATKRead> loadAllReads() {
        return StreamSupport.stream(Spliterators.spliteratorUnknownSize(iterator(), Spliterator.ORDERED), false).collect(Collectors.toList());
    }

    /**
     * Divide an interval into ReadShards. Each shard will cover up to shardSize bases, include shardPadding
     * bases of extra padding on either side, and begin shardSize bases after the previous shard (ie., shards will
     * not overlap except potentially in the padded regions).
     *
     * @param interval interval to shard; must be on the contig according to the provided dictionary
     * @param shardSize desired shard size; intervals larger than this will be divided into shards of up to this size
     * @param shardPadding desired shard padding; each shard's interval will be padded on both sides by this number of bases (may be 0)
     * @param readsSource data source for reads
     * @param dictionary sequence dictionary for reads
     * @return List of {@link ReadShard} objects spanning the interval
     */
    public static List<ReadShard> divideIntervalIntoShards( final SimpleInterval interval, final int shardSize, final int shardPadding, final ReadsDataSource readsSource, final SAMSequenceDictionary dictionary ) {
        return divideIntervalIntoShards(interval, shardSize, shardSize, shardPadding, readsSource, dictionary);
    }

    /**
     * Divide an interval into ReadShards. Each shard will cover up to shardSize bases, include shardPadding
     * bases of extra padding on either side, and begin shardStep bases after the previous shard.
     *
     * @param interval interval to shard; must be on the contig according to the provided dictionary
     * @param shardSize desired shard size; intervals larger than this will be divided into shards of up to this size
     * @param shardStep each shard will begin this many bases away from the previous shard
     * @param shardPadding desired shard padding; each shard's interval will be padded on both sides by this number of bases (may be 0)
     * @param readsSource data source for reads
     * @param dictionary sequence dictionary for reads
     * @return List of {@link ReadShard} objects spanning the interval
     */
    public static List<ReadShard> divideIntervalIntoShards( final SimpleInterval interval, final int shardSize, final int shardStep, final int shardPadding, final ReadsDataSource readsSource, final SAMSequenceDictionary dictionary ) {
        Utils.nonNull(interval);
        Utils.nonNull(readsSource);
        Utils.nonNull(dictionary);
        Utils.validateArg(shardSize >= 1, "shardSize must be >= 1");
        Utils.validateArg(shardStep >= 1, "shardStep must be >= 1");
        Utils.validateArg(shardPadding >= 0, "shardPadding must be >= 0");

        if ( ! IntervalUtils.intervalIsOnDictionaryContig(interval, dictionary) ) {
            throw new IllegalArgumentException("Interval " + interval + " not within the bounds of a contig in the provided dictionary");
        }

        final List<ReadShard> shards = new ArrayList<>();
        int start = interval.getStart();

        while ( start <= interval.getEnd() ) {
            int end = Math.min(start + shardSize - 1, interval.getEnd());

            final SimpleInterval nextShardInterval = new SimpleInterval(interval.getContig(), start, end);
            final SimpleInterval nextShardIntervalPadded = nextShardInterval.expandWithinContig(shardPadding, dictionary);
            shards.add(new ReadShard(nextShardInterval, nextShardIntervalPadded, readsSource));

            start += shardStep;
        }

        return shards;
    }
}
