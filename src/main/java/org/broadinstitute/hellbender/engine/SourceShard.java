package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.iterators.FilteringIterator;

import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * A class to represent a shard of data, optionally expanded by a configurable amount of padded data.
 *
 * The data is lazily loaded by default from a {@link GATKDataSource<T>} (when accessing via {@link #iterator}.
 * Loading all the data in the shard at once is possible via {@link #loadAllData}.
 *
 * The data returned will overlap the expanded padded interval. It's possible to query whether they are within
 * the main part of the shard via {@link #contains} and {@link #containsStartPosition}.
 */
public class SourceShard<T> implements Iterable<T>, Locatable{

    private final SimpleInterval interval;
    private final SimpleInterval paddedInterval;
    private final GATKDataSource<T> dataSource;
    private Predicate<T> filter;

    /**
     * Create a new SourceShard spanning the specified interval, with the specified amount of padding.
     *
     * @param interval the genomic span covered by this shard
     * @param paddedInterval the span covered by this shard, plus any additional padding on each side (must contain the un-padded interval)
     * @param dataSource source of data from which to populate this shard
     */
    protected SourceShard(final SimpleInterval interval, final SimpleInterval paddedInterval, final GATKDataSource<T> dataSource ) {
        Utils.nonNull(interval);
        Utils.nonNull(paddedInterval);
        Utils.nonNull(dataSource);
        Utils.validateArg(paddedInterval.contains(interval), "The padded interval must contain the un-padded interval");

        this.interval = interval;
        this.paddedInterval = paddedInterval;
        this.dataSource = dataSource;
    }

    /**
     * Create a new SourceShard spanning the specified interval, with the specified amount of padding.
     *
     * @param interval the genomic span covered by this shard
     * @param paddedInterval the span covered by this shard, plus any additional padding on each side (must contain the un-padded interval)
     * @param dataSource source of data from which to populate this shard
     */
    public static <T> SourceShard<T> create(final SimpleInterval interval, final SimpleInterval paddedInterval, final GATKDataSource<T> dataSource) {
        return new SourceShard<>(interval, paddedInterval, dataSource);
    }

    /**
     * Create a new SourceShard spanning the specified interval, with no additional padding
     *
     * @param interval the genomic span covered by this shard
     * @param dataSource source of data from which to populate this shard
     */
    public static <T> SourceShard<T> create(final SimpleInterval interval, final GATKDataSource<T> dataSource) {
        return new SourceShard<>(interval, interval, dataSource);
    }

    /**
     * Source data in this shard will be filtered using this filter before being returned
     *
     * @param filter filter to use (may be null)
     */
    public void setFilter( final Predicate<T> filter ) {
        this.filter = filter;
    }

    /**
     * @return an iterator over objects in this shard, as filtered using the configured filter;
     *         objects are lazily loaded rather than pre-loaded
     */
    @Override
    public Iterator<T> iterator() {
        final Iterator<T> it = dataSource.query(paddedInterval);
        return filter != null ? new FilteringIterator<>(it, filter) : it;
    }

    @Override
    public String getContig() {
        return interval.getContig();
    }

    @Override
    public int getStart() {
        return interval.getStart();
    }

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
     * @return a List containing all data in this shard, pre-loaded
     *
     * note: call {@link #iterator} instead to avoid pre-loading all data at once
     */
    public List<T> loadAllData() {
        return StreamSupport.stream(Spliterators.spliteratorUnknownSize(iterator(), Spliterator.ORDERED), false).collect(Collectors.toList());
    }

    /**
     * Divide an interval into SourceShards. Each shard will cover up to shardSize bases, include shardPadding
     * bases of extra padding on either side, and begin shardSize bases after the previous shard (ie., shards will
     * not overlap except potentially in the padded regions).
     *
     * @param interval interval to shard; must be on the contig according to the provided dictionary
     * @param shardSize desired shard size; intervals larger than this will be divided into shards of up to this size
     * @param shardPadding desired shard padding; each shard's interval will be padded on both sides by this number of bases (may be 0)
     * @param dataSource source of data
     * @param dictionary sequence dictionary for data
     * @return List of {@link SourceShard<T>} objects spanning the interval
     */
    public static <T> List<SourceShard<T>> divideIntervalIntoShards( final SimpleInterval interval, final int shardSize, final int shardPadding, final GATKDataSource<T> dataSource, final SAMSequenceDictionary dictionary ) {
        return divideIntervalIntoShards(interval, shardSize, shardSize, shardPadding, dataSource, dictionary);
    }

    /**
     * Divide an interval into SourceShards. Each shard will cover up to shardSize bases, include shardPadding
     * bases of extra padding on either side, and begin shardStep bases after the previous shard.
     *
     * @param interval interval to shard; must be on the contig according to the provided dictionary
     * @param shardSize desired shard size; intervals larger than this will be divided into shards of up to this size
     * @param shardStep each shard will begin this many bases away from the previous shard
     * @param shardPadding desired shard padding; each shard's interval will be padded on both sides by this number of bases (may be 0)
     * @param dataSource source of data
     * @param dictionary sequence dictionary for data
     * @return List of {@link SourceShard<T>} objects spanning the interval
     */
    public static <T> List<SourceShard<T>> divideIntervalIntoShards( final SimpleInterval interval, final int shardSize, final int shardStep, final int shardPadding, final GATKDataSource<T> dataSource, final SAMSequenceDictionary dictionary ) {
        Utils.nonNull(interval);
        Utils.nonNull(dataSource);
        Utils.nonNull(dictionary);
        Utils.validateArg(shardSize >= 1, "shardSize must be >= 1");
        Utils.validateArg(shardStep >= 1, "shardStep must be >= 1");
        Utils.validateArg(shardPadding >= 0, "shardPadding must be >= 0");

        if ( ! IntervalUtils.intervalIsOnDictionaryContig(interval, dictionary) ) {
            throw new IllegalArgumentException("Interval " + interval + " not within the bounds of a contig in the provided dictionary");
        }

        final List<SourceShard<T>> shards = new ArrayList<>();
        int start = interval.getStart();

        while ( start <= interval.getEnd() ) {
            int end = Math.min(start + shardSize - 1, interval.getEnd());

            final SimpleInterval nextShardInterval = new SimpleInterval(interval.getContig(), start, end);
            final SimpleInterval nextShardIntervalPadded = nextShardInterval.expandWithinContig(shardPadding, dictionary);
            shards.add(create(nextShardInterval, nextShardIntervalPadded, dataSource));

            start += shardStep;
        }

        return shards;
    }

}
