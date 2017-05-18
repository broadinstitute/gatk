package org.broadinstitute.hellbender.utils.iterators;

import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Iterator that will break up each input interval into shards.  The advantage of this iterator is that it is not RAM
 *  intensive.
 *
 *  Empty iterator of intervals is supported.
 *
 * Developer note:  This class queues up the next interval to return and then has (minimal) code to detect when
 *  an empty iterator was passed in for initialization.
 */
public class ShardedIntervalIterator implements Iterator<SimpleInterval> {

    private Iterator<SimpleInterval> intervals;
    private int shardSize;

    /**
     * The input interval that we are currently sharding.
     */
    private SimpleInterval currentInterval;

    private int currentOffsetInCurrentInterval;

    private int lastOffsetInCurrentInterval;
    private SimpleInterval shardedInterval;

    public ShardedIntervalIterator(Iterator<SimpleInterval> intervals, int shardSizeInBases) {
        Utils.validate(shardSizeInBases > 0, "Invalid shard size.  Must be greater than zero.");
        this.intervals = intervals;
        this.shardSize = shardSizeInBases;
        this.currentOffsetInCurrentInterval = 0;
        this.lastOffsetInCurrentInterval = 0;
        this.currentInterval = null;
        this.shardedInterval = null;


        // Queue up the next interval to shard or keep it as null if the interval iterator is exhausted.
        advanceInterval();
    }

    /**
     * This should only be called when the shards within an interval have been exhausted.  I.e. when starting a new
     *   interval to be sharded.
     */
    private void advanceInterval() {

        if (this.intervals.hasNext()) {
            currentInterval = this.intervals.next();
            currentOffsetInCurrentInterval = IntervalUtils.shardIndex(1, shardSize);
            lastOffsetInCurrentInterval = IntervalUtils.shardIndex(this.currentInterval.size(), shardSize);
            shardedInterval = calculateShardedInterval();

        } else {
            lastOffsetInCurrentInterval = 0;
            currentInterval = null;
            shardedInterval = null;
        }
    }

    private SimpleInterval calculateShardedInterval() {
        return new SimpleInterval(this.currentInterval.getContig(), currentInterval.getStart() + IntervalUtils.beginOfShard(currentOffsetInCurrentInterval, shardSize) - 1,
                Integer.min(currentInterval.getStart() + IntervalUtils.endOfShard(currentOffsetInCurrentInterval, shardSize) - 1, currentInterval.getEnd()));
    }

    @Override
    public boolean hasNext() {
        // if the current shard has exhausted the current interval AND there is no next interval, return false

        if (shardedInterval != null) {
            return true;
        }
        return false;
    }

    @Override
    public SimpleInterval next() {

        if (shardedInterval == null) {
            throw new NoSuchElementException();
        }

        final SimpleInterval result = shardedInterval;

        // Advance the shard index and (if necessary) set it back to zero and get the next interval.
        advanceShardInInterval();

        return result;
    }

    /**
     * This can also advance the interval if needed.
     */
    private void advanceShardInInterval() {
        currentOffsetInCurrentInterval ++;
        if (currentOffsetInCurrentInterval > lastOffsetInCurrentInterval) {

            // Advance the interval.
            advanceInterval();
        }
        if (currentInterval != null) {
            shardedInterval = calculateShardedInterval();
        } else {
            shardedInterval = null;
        }
    }
}
