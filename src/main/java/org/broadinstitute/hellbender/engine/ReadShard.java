package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.downsampling.ReadsDownsampler;
import org.broadinstitute.hellbender.utils.downsampling.ReadsDownsamplingIterator;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * A class to represent a shard of reads data, optionally expanded by a configurable amount of padded data.
 *
 * The reads are lazily loaded by default (when accessing the reads via {@link #iterator}. Loading all the
 * reads in the shard at once is possible via {@link #loadAllReads} or {@link #loadAllData()}
 *
 * The reads returned will overlap the expanded padded interval. It's possible to query whether they are within
 * the main part of the shard via {@link #contains} and {@link #containsStartPosition}.
 *
 * The reads in the shard can be filtered via {@link #setReadFilter} (no filtering is performed by default).
 */
public final class ReadShard extends SourceShard<GATKRead> {

    private ReadsDownsampler downsampler;

    public ReadShard(SimpleInterval interval, SimpleInterval paddedInterval, ReadsDataSource readsSource) {
        super(interval, paddedInterval, readsSource);
    }

    public ReadShard(SimpleInterval interval, ReadsDataSource readsSource) {
        super(interval, interval, readsSource);
    }

    /**
     * Reads in this shard will be filtered using this filter before being returned.
     * Read filtering will be performed before any requested downsampling.
     *
     * @param filter filter to use (may be null, which signifies that no filtering is to be performed)
     */
    public void setReadFilter( final ReadFilter filter ) {
        setFilter(filter);;
    }

    /**
     * Reads in this shard will be downsampled using this downsampler before being returned.
     * Downsampling will be performed after any requested read filtering.
     *
     * @param downsampler downsampler to use (may be null, which signifies that no downsampling is to be performed)
     */
    public void setDownsampler( final ReadsDownsampler downsampler ) {
        this.downsampler = downsampler;
    }

    /**
     * @return an iterator over reads in this shard, as filtered using the configured read filter
     *         and downsampled using the configured downsampler; reads are lazily loaded rather than pre-loaded
     *
     * Note that any read filtering is always performed before any downsampling.
     */
    @Override
    public Iterator<GATKRead> iterator() {
        return ( downsampler == null ) ? super.iterator() : new ReadsDownsamplingIterator(super.iterator(), downsampler);
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
        return loadAllData();
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

        Utils.validateArg(IntervalUtils.intervalIsOnDictionaryContig(interval, dictionary), () ->
                "Interval " + interval + " not within the bounds of a contig in the provided dictionary");

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
