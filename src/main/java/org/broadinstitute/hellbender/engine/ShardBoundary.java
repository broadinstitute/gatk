package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;
import java.util.Objects;

/**
 * Holds the bounds of a {@link Shard}, both with and without padding
 */
public class ShardBoundary implements Locatable, Serializable {
    private static final long serialVersionUID = 1L;

    protected final SimpleInterval interval;
    protected final SimpleInterval paddedInterval;
    protected final boolean padded;

    /**
     * Create a new ShardBoundary from the given intervals
     *
     * @param interval the interval covered by the shard
     * @param paddedInterval the interval covered by the shard's padding, must contain the shard interval
     */
    public ShardBoundary(final SimpleInterval interval, final SimpleInterval paddedInterval) {
        this(interval, paddedInterval, false);
    }

    protected ShardBoundary(final SimpleInterval interval, final SimpleInterval paddedInterval, final boolean padded) {
        Utils.nonNull(interval);
        Utils.nonNull(paddedInterval);
        Utils.validateArg(paddedInterval.contains(interval), "interval must be contained within paddedInterval");
        this.interval = interval;
        this.paddedInterval = paddedInterval;
        this.padded = padded;
    }

    @Override
    public String getContig() {
        return (padded ? paddedInterval : interval).getContig();
    }

    /**
     * @return start of the shard boundary without padding
     */
    @Override
    public int getStart() {
        return (padded ? paddedInterval : interval).getStart();
    }

    /**
     * @return end of the shard boundary without padding
     */
    @Override
    public int getEnd() {
        return (padded ? paddedInterval : interval).getEnd();
    }

    /**
     * @return the interval this boundary covers without padding
     */
    public SimpleInterval getInterval() {
        return interval;
    }

    /**
     * @return the interval this boundary covers including padding
     */
    public SimpleInterval getPaddedInterval() {
        return paddedInterval;
    }

    public ShardBoundary paddedShardBoundary() {
        return padded ? this : new ShardBoundary(interval, paddedInterval, true);
    }

    public <T> Shard<T> createShard(Iterable<T> locatables) {
        return new ShardBoundaryShard<>(this, locatables);
    }

    @Override
    public boolean equals(Object o) {
        if ( this == o ) {
            return true;
        }
        if ( o == null || getClass() != o.getClass() ) {
            return false;
        }

        final ShardBoundary key = (ShardBoundary) o;

        if ( !Objects.equals(interval, key.interval) ) {
            return false;
        }
        return Objects.equals(paddedInterval, key.paddedInterval);

    }

    @Override
    public int hashCode() {
        return Objects.hash(interval, paddedInterval);
    }

    @Override
    public String toString() {
        return "ShardBoundary{" +
                "interval=" + interval +
                ", paddedInterval=" + paddedInterval +
                ", padded=" + padded +
                '}';
    }
}
