package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.List;


/**
 * A Shard of records of type T covering a specific genomic interval, optionally expanded by a configurable
 * amount of padded data, that provides the ability to iterate over its records.
 */
public interface Shard<T> extends Iterable<T>, Locatable {

    /**
     * @return the interval this shard spans
     */
    SimpleInterval getInterval();

    /**
     * @return the interval this shard spans, potentially with additional padding on each side
     * it must be the case that for a given Shard getPaddedInterval().contains(getInterval())
     */
    SimpleInterval getPaddedInterval();

    /**
     * @return the start of the non-padded interval this shard covers
     */
    @Override
    default int getStart() {
        return getInterval().getStart();
    }

    /**
     * @return the end of the non-padded interval this shard covers
     */
    @Override
    default int getEnd() {
        return getInterval().getEnd();
    }

    /**
     * @return contig this shard belongs to
     */
    @Override
    default String getContig() {
        return getInterval().getContig();
    }

    /**
     * Divide an interval into ShardBoundaries. Each shard will cover up to shardSize bases, include shardPadding
     * bases of extra padding on either side, and begin shardSize bases after the previous shard (ie., shards will
     * not overlap except potentially in the padded regions).
     *
     * @param interval interval to shard; must be on the contig according to the provided dictionary
     * @param shardSize desired shard size; intervals larger than this will be divided into shards of up to this size
     * @param shardPadding desired shard padding; each shard's interval will be padded on both sides by this number of bases (may be 0)
     * @param dictionary sequence dictionary for reads
     * @return List of {@link ShardBoundary} objects spanning the interval
     */
    static List<ShardBoundary> divideIntervalIntoShards(final SimpleInterval interval, final int shardSize, final int shardPadding, final SAMSequenceDictionary dictionary) {
        return  divideIntervalIntoShards(interval, shardSize, shardSize, shardPadding, dictionary);
    }

    /**
     * Divide an interval into ShardBoundaries. Each shard will cover up to shardSize bases, include shardPadding
     * bases of extra padding on either side, and begin shardStep bases after the previous shard.
     *
     * @param interval interval to shard; must be on the contig according to the provided dictionary
     * @param shardSize desired shard size; intervals larger than this will be divided into shards of up to this size
     * @param shardStep each shard will begin this many bases away from the previous shard
     * @param shardPadding desired shard padding; each shard's interval will be padded on both sides by this number of bases (may be 0)
     * @param dictionary sequence dictionary for reads
     * @return List of {@link ShardBoundary} objects spanning the interval
     */
    static List<ShardBoundary> divideIntervalIntoShards(final SimpleInterval interval, final int shardSize, final int shardStep, final int shardPadding, final SAMSequenceDictionary dictionary) {
        Utils.nonNull(interval);
        Utils.nonNull(dictionary);
        Utils.validateArg(shardSize >= 1, "shardSize must be >= 1");
        Utils.validateArg(shardStep >= 1, "shardStep must be >= 1");
        Utils.validateArg(shardPadding >= 0, "shardPadding must be >= 0");

        Utils.validateArg(IntervalUtils.intervalIsOnDictionaryContig(interval, dictionary), () ->
                "Interval " + interval + " not within the bounds of a contig in the provided dictionary");

        final List<ShardBoundary> shards = new ArrayList<>();
        int start = interval.getStart();

        while ( start <= interval.getEnd() ) {
            final int end = Math.min(start + shardSize - 1, interval.getEnd());
            final SimpleInterval nextShardInterval = new SimpleInterval(interval.getContig(), start, end);
            final SimpleInterval nextShardIntervalPadded = nextShardInterval.expandWithinContig(shardPadding, dictionary);
            shards.add(new ShardBoundary(nextShardInterval, nextShardIntervalPadded));
            start += shardStep;
        }

        return shards;
    }
}
