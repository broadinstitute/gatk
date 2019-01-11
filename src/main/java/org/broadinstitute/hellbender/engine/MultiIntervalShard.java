package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.List;

/**
 * An interface to represent shards of arbitrary data spanning multiple intervals.
 *
 * @param <T> Type of data in this shard
 */
public interface MultiIntervalShard<T> extends Iterable<T> {

    /**
     * @return A List of this shard's intervals
     */
    List<SimpleInterval> getIntervals();

    /**
     * @return A List of this shard's intervals, with padding added to each interval on both sides
     */
    List<SimpleInterval> getPaddedIntervals();
}
