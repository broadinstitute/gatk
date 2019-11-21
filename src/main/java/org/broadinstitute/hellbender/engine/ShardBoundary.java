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
    protected final SimpleInterval paddedSpan;
    protected final boolean padded;

    /**
     * Create a new ShardBoundary from the given intervals
     *
     * @param interval the interval covered by the shard
     * @param paddedSpan the interval covered by the shard's padding, must contain the shard interval
     */
    public ShardBoundary(final SimpleInterval interval, final SimpleInterval paddedSpan) {
        this(interval, paddedSpan, false);
    }

    protected ShardBoundary(final SimpleInterval interval, final SimpleInterval paddedSpan, final boolean padded) {
        Utils.nonNull(interval);
        Utils.nonNull(paddedSpan);
        Utils.validateArg(paddedSpan.contains(interval), "interval must be contained within paddedSpan");
        this.interval = interval;
        this.paddedSpan = paddedSpan;
        this.padded = padded;
    }

    @Override
    public String getContig() {
        return (padded ? paddedSpan : interval).getContig();
    }

    /**
     * @return start of the shard boundary without padding
     */
    @Override
    public int getStart() {
        return (padded ? paddedSpan : interval).getStart();
    }

    /**
     * @return end of the shard boundary without padding
     */
    @Override
    public int getEnd() {
        return (padded ? paddedSpan : interval).getEnd();
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
    public SimpleInterval getPaddedSpan() {
        return paddedSpan;
    }

    public ShardBoundary paddedShardBoundary() {
        return padded ? this : new ShardBoundary(interval, paddedSpan, true);
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
        return Objects.equals(paddedSpan, key.paddedSpan);

    }

    @Override
    public int hashCode() {
        return Objects.hash(interval, paddedSpan);
    }

    @Override
    public String toString() {
        return "ShardBoundary{" +
                "interval=" + interval +
                ", paddedSpan=" + paddedSpan +
                ", padded=" + padded +
                '}';
    }
}
